/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 					  jlvilas@cnb.csic.es
 * 			   Carlos Oscar S. Sorzano            coss@cnb.csic.es (2016)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "resolution_monogenic_signal.h"

#include <data/args.h>
#include <data/projection.h>
#include <data/xmipp_fftw.h>
//#include <data/xmipp_fft.h>
#include <data/metadata_extension.h>
#include "reconstruct_art.h"
#include "fourier_projection.h"
#include <vector>
#include <math.h>
#include <limits>
#include <complex>


void ProgMonogenicSignalRes::readParams()
{
	fnVol = getParam("--vol");
	fnTilt = getParam("--tiltparticles");
	fnDir = getParam("--odir");
	ang_dist = getDoubleParam("--angular_distance");
	ang_acc = getDoubleParam("--angular_accuracy");
	N_cls = getIntParam("--classes");
}


void ProgMonogenicSignalRes::defineParams()
{
	//usage
	addUsageLine("This function determines the local resolution of a map");
	//params
	addParamsLine("  [--vol <md_file=\"\">]    : Input volume");
	addParamsLine("  [--tiltparticles <md_file=\"\">]    : Tilt particles stack");
	addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
	addParamsLine("  [--angular_distance <s=5>]   : Angular distance");
	addParamsLine("  [--angular_accuracy <s=5>]   : Angular distance");
	addParamsLine("  [--classes <s=2>]   : Number of classes");
}


void ProgMonogenicSignalRes::RieszTransform3Dreal(const MultidimArray<double> &inputVol,
													std::vector<MultidimArray< std::complex<double> > > &RieszVector)
{
	MultidimArray< std::complex<double> > RieszVector0, RieszVector1, RieszVector2;

	MultidimArray< std::complex<double> > fftV;

	FourierTransformer transformer, transformer_inv;
	MultidimArray<double> I, ifftV, inputVol_aux = inputVol;

	transformer.FourierTransform(inputVol_aux,fftV, false);

	RieszVector0.initZeros(fftV);
	RieszVector1.initZeros(fftV);
	RieszVector2.initZeros(fftV);

//	std::cout << std::numeric_limits<double>::epsilon() << std::endl;
	double uz, uy, ux, uz2, u2, uz2y2, u;
	for(size_t k=0; k<ZSIZE(fftV); ++k)
	{
		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol),uz);
		uz2=uz*uz;
		for(size_t i=0; i<YSIZE(fftV); ++i)
		{
			FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
			uz2y2=uz2+uy*uy;
			for(size_t j=0; j<XSIZE(fftV); ++j)
			{
				FFT_IDX2DIGFREQ(j,XSIZE(inputVol),ux);
				u2=uz2y2+ux*ux;
				u = sqrt(u2);
				if ((k == 0) && (i == 0) && (j == 0))
				{
					u = std::numeric_limits<double>::epsilon();
				}
				A3D_ELEM(RieszVector0, k ,i, j) = (ux/u)*A3D_ELEM(fftV, k ,i, j);
				A3D_ELEM(RieszVector1, k ,i, j) = (uy/u)*A3D_ELEM(fftV, k ,i, j);
				A3D_ELEM(RieszVector2, k ,i, j) = (uz/u)*A3D_ELEM(fftV, k ,i, j);
			}
		}
	}

	MultidimArray<std::complex <double> >aux0, aux1, aux2;

	aux0.initZeros(fftV);
	aux1.initZeros(fftV);
	aux2.initZeros(fftV);

	std::cout << aux0 << std::endl;

	std::cout << "---------------------------" << std::endl;

	transformer_inv.inverseFourierTransform(RieszVector0, aux0);
	transformer_inv.inverseFourierTransform(RieszVector1, aux1);
	transformer_inv.inverseFourierTransform(RieszVector2, aux2);

	std::cout << aux0 << std::endl;

	RieszVector[0] = aux0;
	RieszVector[1] = aux1;
	RieszVector[2] = aux2;
}


void ProgMonogenicSignalRes::amplitudeMonogenicSignal(const MultidimArray<double> &inputVol,
													  const std::vector<MultidimArray< std::complex<double> > > &RieszVector,
															MultidimArray<double> &amplitude)
{
	const MultidimArray< std::complex<double> > &v0=RieszVector[0], &v1=RieszVector[1], &v2=RieszVector[2];
	amplitude.initZeros(inputVol);
	FOR_ALL_ELEMENTS_IN_ARRAY3D(inputVol)
	{
		A3D_ELEM(amplitude, k, i, j) = sqrt((double) std::norm(A3D_ELEM(v0, k, i, j)) +
											(double) std::norm(A3D_ELEM(v1, k, i, j)) +
											(double) std::norm(A3D_ELEM(v2, k, i, j)) +
											A3D_ELEM(inputVol, k, i, j)*A3D_ELEM(inputVol, k, i, j));
	}
}


void ProgMonogenicSignalRes::passbandfiltervol(const MultidimArray<double> &inputVol, double freq, MultidimArray<double> &filteredVol, double FWHM)
{
	filteredVol = inputVol;

	const MultidimArray<double> &v_flt=filteredVol;
	int idx;

	//Full Width at Half Maximum FWHM = s*sqrt(2*log(2)) *sigma = 2.35482*sigma;
	double sigma = FWHM/2.35482;
	double sigma2 = sigma * sigma;
	double K = 1. / pow(sqrt(2.*PI)*sigma,(double)inputVol.getDim());
	double r2;

	DIGFREQ2FFT_IDX(freq, ZSIZE(inputVol), idx);

	FOR_ALL_ELEMENTS_IN_ARRAY3D(inputVol)
	{
		r2 = (k - idx) * (k - idx) + (i - idx) * (i - idx) + (j - idx) * (j - idx);
		A3D_ELEM(filteredVol, k, i, j) = K * exp(-0.5 * r2 / sigma2);
	}
}


void ProgMonogenicSignalRes::run()
{
	std::cout << "Starting ... " << std::endl;
	std::vector< MultidimArray< std::complex <double> > > RieszVector(3);

	Image<double> imgVol, amplitudeMS;
	imgVol.read(fnVol);

	std::cout << "Riesz Transform " << std::endl;
	std::cout << "Riesz vector size = " << RieszVector.size() << std::endl;
	RieszTransform3Dreal(imgVol(), RieszVector);

	std::cout << "amplitude Monogenic Signal " << std::endl;
	amplitudeMonogenicSignal(imgVol(), RieszVector, amplitudeMS());

	FileName amplitude_out = "amplitude.vol";

	std::cout << "before saving " << std::endl;
	amplitudeMS.write(amplitude_out);
	std::cout << "saved! " << std::endl;

	double step = 0.05;

	int count = 0;
	size_t N_elems = (size_t) (0.5/step);
	std::cout <<  "Number of elements =  " << N_elems << std::endl;
	std::vector< MultidimArray <double> > filteredVol(N_elems);


	for (size_t idx_freq = 0; idx_freq<N_elems; idx_freq++)
	{
		filteredVol[idx_freq].initZeros(amplitudeMS());
		passbandfiltervol(amplitudeMS(), idx_freq*step, filteredVol[idx_freq]);
		count++;
		std::cout <<  "idx_freq =  " << idx_freq << std::endl;
	}

	std::cout <<  "Finished!" << std::endl;

	//exit(0);


}


