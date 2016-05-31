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
#include <data/morphology.h>
#include <data/args.h>
#include <data/projection.h>
#include <data/xmipp_fftw.h>
#include <data/xmipp_fft.h>
#include <data/metadata_extension.h>
#include "reconstruct_art.h"
#include "fourier_projection.h"
#include <vector>
#include <math.h>
#include <limits>
#include <complex>
//#define DEBUG
//#define METHOD1
//#define METHOD2
//#define METHOD3
//#define METHOD4
#define TESTING
//#define MAJORIZATION


void ProgMonogenicSignalRes::readParams()
{
	fnVol = getParam("--vol");
	fnDir = getParam("--odir");
	fnMask = getParam("--mask");
	fType = getDoubleParam("--filter_type");
	smpr = getDoubleParam("--sampling_rate");
	condpremask = checkParam("--premask");
}


void ProgMonogenicSignalRes::defineParams()
{
	//usage
	addUsageLine("This function determines the local resolution of a map");
	//params
	addParamsLine("  [--vol <md_file=\"\">]    : Input volume");
	addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
	addParamsLine("  [--mask <md_file=\"\">]   : Mask for locate the volume");
	addParamsLine("  [--filter_type <s=0>]   : pass band or pass high");
	addParamsLine("  [--sampling_rate <s=0>]   : Sampling rate (A/px)");
	addParamsLine("  [--premask]   : Only if the introduced map has already masked");
}


void ProgMonogenicSignalRes::RieszTransform3Dreal(const MultidimArray<double> &inputVol,
													std::vector<MultidimArray<double> > &RieszVector)
{
	MultidimArray< std::complex<double> > RieszVector0, RieszVector1, RieszVector2;

	MultidimArray< std::complex<double> > fftV;

	FourierTransformer transformer, transformer_inv;
	MultidimArray<double> I, ifftV, inputVol_aux = inputVol;

	transformer.FourierTransform(inputVol_aux,fftV, false);

	RieszVector0.initZeros(fftV);
	RieszVector1.initZeros(fftV);
	RieszVector2.initZeros(fftV);

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

	MultidimArray<double> aux0, aux1, aux2;

	aux0.initZeros(inputVol);
	aux1.initZeros(inputVol);
	aux2.initZeros(inputVol);

	transformer_inv.inverseFourierTransform(RieszVector0, aux0);
	transformer_inv.inverseFourierTransform(RieszVector1, aux1);
	transformer_inv.inverseFourierTransform(RieszVector2, aux2);

	RieszVector[0] = aux0;
	RieszVector[1] = aux1;
	RieszVector[2] = aux2;
}


void ProgMonogenicSignalRes::amplitudeMonogenicSignal3D(const MultidimArray<double> &inputVol,
													  const std::vector<MultidimArray<double> > &RieszVector,
															MultidimArray<double> &amplitude)
{
	const MultidimArray<double> &v0=RieszVector[0], &v1=RieszVector[1], &v2=RieszVector[2];
	//amplitude.initZeros(inputVol);
	amplitude.initZeros(RieszVector[0]);

	FOR_ALL_ELEMENTS_IN_ARRAY3D(amplitude)
	{
		A3D_ELEM(amplitude, k, i, j) = sqrt((double) A3D_ELEM(v0, k, i, j)*A3D_ELEM(v0, k, i, j) +
											(double) A3D_ELEM(v1, k, i, j)*A3D_ELEM(v1, k, i, j) +
											(double) A3D_ELEM(v2, k, i, j)*A3D_ELEM(v2, k, i, j) +
											A3D_ELEM(inputVol, k, i, j)*A3D_ELEM(inputVol, k, i, j));
	}
}

void ProgMonogenicSignalRes::medianFilter3x3x3(const MultidimArray<double> &inputVol, MultidimArray<double> &FilteredVol)
{
	FilteredVol.initZeros(inputVol);
	double median;

	const MultidimArray<double> &p2inputVol = inputVol;

	std::vector<double> aux(6);

	std::cout << "median filter ..." << std::endl;

	for (size_t k = 1; k<(ZSIZE(FilteredVol)-1); k++)
	{
		for (size_t i = 1; i<(YSIZE(FilteredVol)-1); i++)
		{
			for (size_t j = 1; j<(XSIZE(FilteredVol)-1); j++)
			{
				aux[0] = A3D_ELEM(inputVol, k+1, i, j);
				aux[1] = A3D_ELEM(inputVol, k-1, i, j);
				aux[2] = A3D_ELEM(inputVol, k, i+1, j);
				aux[3] = A3D_ELEM(inputVol, k, i-1, j);
				aux[4] = A3D_ELEM(inputVol, k, i, j+1);
				aux[5] = A3D_ELEM(inputVol, k, i, j-1);
				//std::cout << " blabla" << std::endl;
				std::sort(aux.begin(), aux.end());

				median = 0.5*(aux[2] + aux[3]);

				A3D_ELEM(FilteredVol, k, i, j) = median;
			}
		}
		//std::cout << "k = " << k << " " << ZSIZE(FilteredVol) << std::endl;
	}

}


void ProgMonogenicSignalRes::FourierAmplitudeMonogenicSignal3D(const MultidimArray<double> &inputVol,
													  const std::vector<MultidimArray<double> > &RieszVector, MultidimArray<double> &amplitude,
													  MultidimArray< std::complex<double> > &fftamplitude)
{
	amplitudeMonogenicSignal3D(inputVol, RieszVector, amplitude);

	FourierTransformer transformer;
	MultidimArray<double> TFamplitude_aux = amplitude;
	TFamplitude_aux.setXmippOrigin();
	transformer.FourierTransform(TFamplitude_aux,fftamplitude, false);
}


void ProgMonogenicSignalRes::passbandfiltervol(const MultidimArray<double> &inputVol, double freq, MultidimArray<double> &filteredVol, double sigma)
{
	filteredVol = inputVol;

	int idx, idx_y_z;
	MultidimArray<double> inputVol_aux = inputVol;
	inputVol_aux.setXmippOrigin();

	FourierTransformer transformer, transformer_inv;
	MultidimArray< std::complex<double> > fftV, fftVfiltered_vol;

	transformer.FourierTransform(inputVol_aux,fftV, false);

	fftVfiltered_vol = fftV;

	double uz, uy, ux, u, invsigma2;
	invsigma2 = 1/(sigma*sigma);

	DIGFREQ2FFT_IDX(freq, ZSIZE(inputVol), idx);

	FOR_ALL_ELEMENTS_IN_ARRAY3D(fftVfiltered_vol)
	{
		FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
		FFT_IDX2DIGFREQ(j,XSIZE(inputVol),ux);
		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol),uz);

		u = sqrt(ux*ux + uy*uy + uz*uz);

		A3D_ELEM(fftVfiltered_vol, k, i, j) = exp(-0.5 *(u-freq)*(u-freq) * invsigma2)*A3D_ELEM(fftV,k,i,j);
	}

	transformer_inv.inverseFourierTransform(fftVfiltered_vol, filteredVol);

}


void ProgMonogenicSignalRes::passbandfiltervolFourier(const MultidimArray<double> &inputVol, MultidimArray< std::complex<double> > fftVol,
																double freq, MultidimArray< std::complex<double> > fftVfiltered_vol, double sigma)
{
	int idx, idx_y_z;

	fftVfiltered_vol = fftVol;

	double uz, uy, ux, u, invsigma2;
	invsigma2 = 1/(sigma*sigma);

	DIGFREQ2FFT_IDX(freq, ZSIZE(inputVol), idx);

	FOR_ALL_ELEMENTS_IN_ARRAY3D(fftVfiltered_vol)
	{
		FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
		FFT_IDX2DIGFREQ(j,XSIZE(inputVol),ux);
		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol),uz);

		u = sqrt(ux*ux + uy*uy + uz*uz);

		A3D_ELEM(fftVfiltered_vol, k, i, j) = exp(-0.5 *(u-freq)*(u-freq) * invsigma2)*A3D_ELEM(fftVol,k,i,j);
	}
}


void ProgMonogenicSignalRes::highpassfiltervol(const MultidimArray<double> &inputVol, double freq, MultidimArray<double> &filteredVol, double raised_w)
{
	filteredVol = inputVol;

	const MultidimArray<double> &v_flt=filteredVol;  //TOCHECK
	int idx, idx_y_z;


	FourierTransformer transformer_inv;
	MultidimArray< std::complex<double> > fftV, fftVfiltered_vol;

	FourierTransformer transformer;
	MultidimArray<double> inputVol_aux = inputVol;
	inputVol_aux.setXmippOrigin();
	transformer.FourierTransform(inputVol_aux,fftV, false);

	fftVfiltered_vol = fftV;

	double uz, uy, ux, u;

	DIGFREQ2FFT_IDX(freq, ZSIZE(inputVol), idx);

	FOR_ALL_ELEMENTS_IN_ARRAY3D(fftVfiltered_vol)
	{
		FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
		FFT_IDX2DIGFREQ(j,XSIZE(inputVol),ux);
		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol),uz);

		u = sqrt(ux*ux + uy*uy + uz*uz);
//		if (u>freq)
//		{
//			A3D_ELEM(fftVfiltered_vol, k, i, j) = A3D_ELEM(fftV,k,i,j);
//		}
//		else{
//			A3D_ELEM(fftVfiltered_vol, k, i, j) = 0;
//		}
	    if (u>freq)
	    	A3D_ELEM(fftVfiltered_vol, k, i, j) =  A3D_ELEM(fftV,k,i,j);
	    else if (u>freq-raised_w)
	    	A3D_ELEM(fftVfiltered_vol, k, i, j) = ( (1+cos((PI/raised_w)*(freq-u)))/2 )*(A3D_ELEM(fftV,k,i,j));
	    else
	    	A3D_ELEM(fftVfiltered_vol, k, i, j) = 0;
	}

	transformer_inv.inverseFourierTransform(fftVfiltered_vol, filteredVol);
}


void ProgMonogenicSignalRes::applymask(const MultidimArray<double> &amplitudeMS, const MultidimArray<double> &mask,
																					  MultidimArray<double> &mask_amplitudeMS)
{
	mask_amplitudeMS.initZeros(amplitudeMS);
	FOR_ALL_ELEMENTS_IN_ARRAY3D(amplitudeMS)
	{
		A3D_ELEM(mask_amplitudeMS, k, i, j) = A3D_ELEM(amplitudeMS, k, i, j)*A3D_ELEM(mask, k, i, j);
	}
}


void ProgMonogenicSignalRes::run()
{
#ifdef METHOD1
	std::cout << "Starting ... " << std::endl;
	std::vector< MultidimArray<double > > RieszVector(3), RieszVector_original(3);
	MultidimArray<double> mask_amplitudeMS, mask_amplitudeMS_original;
	Image<double> imgVol, amplitudeMS, mask, outputResolution, amplitudeMS_original;
	double idx_threshold_energy, idx_threshold_noise, threshold_energy, threshold_noise, step = 0.01;
	size_t len_energy_vector, len_noise_vector, count_freq, N_elems = (size_t) (0.5/step);
	MultidimArray <double> filteredVol;

	//Reading input Volume
	imgVol.read(fnVol);
	std::cout << "Volume read!" << std::endl;

	size_t zdim = ZSIZE(imgVol())*0.5; //radius
	size_t ydim = YSIZE(imgVol())*0.5;
	size_t xdim = XSIZE(imgVol())*0.5;

	// Calculating Riesz transform

	int count = 0;
	double r2, r_kij;
	if (xdim<ydim){
		if (xdim<zdim){
			r2 = xdim*xdim;}
		else{
			r2 = zdim*zdim;}}
	else{
		if (ydim<zdim){
			r2 = ydim*ydim;}
		else{
			r2 = zdim*zdim;}
	}

	std::cout << "r2 = " << r2 << std::endl;

	outputResolution().initZeros(imgVol());
	mask.read(fnMask);

	MultidimArray<double> &p2imgVol = imgVol(), &p2amplitudeMS = amplitudeMS(), &p2mask = mask(), p2amplitudeMS_original = amplitudeMS_original();

	/////////////
	RieszTransform3Dreal(imgVol(), RieszVector_original);

	amplitudeMonogenicSignal3D(imgVol(), RieszVector_original, p2amplitudeMS_original);

	size_t sum_mask = 0;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(p2mask)
	{
		if (A3D_ELEM(p2mask, k, i, j) > 0)
		sum_mask++;
	}

	std::cout << "Analysing frequencies..." << std::endl;
	std::cout << "Sampling rate = " << smpr << std::endl;
	for (size_t idx_freq = 1; idx_freq<N_elems; idx_freq++)
	{
		if (count_freq == (sum_mask))
			break;

		count_freq = 0;
		std::vector<double> noise_vector(1), energy_vector(1), noise_vector_sorted, energy_vector_sorted;

		filteredVol.initZeros(p2amplitudeMS);
		if (fType == 0.0)
			passbandfiltervol(p2imgVol, idx_freq*step, filteredVol, 0.01);
		else
		{
			if (fType == 1.0)
			{
				highpassfiltervol(p2imgVol, idx_freq*step, filteredVol, 0.01);
			}
			else
			{
				std::cout << "Error: Type of filter was not provided" << std::endl;
				exit(0);
			}
		}

		RieszTransform3Dreal(filteredVol, RieszVector);

		amplitudeMonogenicSignal3D(filteredVol, RieszVector, p2amplitudeMS);

		applymask(p2amplitudeMS, p2mask, mask_amplitudeMS);

		FOR_ALL_ELEMENTS_IN_ARRAY3D(p2amplitudeMS)
		{
			if (A3D_ELEM(p2mask, k, i, j) > 0 )
			{
				energy_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
			}
			else
			{
				if (condpremask == true)
				{
					r_kij = (k-zdim)*(k-zdim) + (j-ydim)*(j-ydim) + (i-xdim)*(i-xdim);
					if (r_kij<r2){
						noise_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
						//noise_vector.push_back(A3D_ELEM(p2amplitudeMS_original, k, i, j));
					}
				}
				else{
					noise_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
					//noise_vector.push_back(A3D_ELEM(p2amplitudeMS_original, k, i, j));
				}
			}
		}


		std::sort(energy_vector.begin(), energy_vector.end());
		std::sort(noise_vector.begin(), noise_vector.end());


		//len_energy_vector = energy_vector.size();
		len_noise_vector = noise_vector.size();


		//idx_threshold_energy = floor((double) 0.95*len_energy_vector);
		idx_threshold_noise = floor((double) 0.95*len_noise_vector);
		std::cout << "len_noise = " << len_noise_vector << std::endl;


		//threshold_energy = energy_vector[idx_threshold_energy];
		threshold_noise = noise_vector[idx_threshold_noise];
		std::cout << "threshold_noise = " << threshold_noise << std::endl;

		double inv_idx_freq = (double) 1/(idx_freq*step);
		std::cout << "freq = " << idx_freq << std::endl;
		std::cout << "inversa freq = " << inv_idx_freq << std::endl;
		FOR_ALL_ELEMENTS_IN_ARRAY3D(outputResolution())
		{
			if (A3D_ELEM(mask_amplitudeMS, k, i, j) > threshold_noise)
			{
				A3D_ELEM(outputResolution(), k, i, j) = smpr*inv_idx_freq;
			}
			else
			{
				count_freq++;
			}
		}
	}

	Image<double> outputResolution2;

	medianFilter3x3x3(outputResolution(), outputResolution2());


	outputResolution2.write(fnDir);

	std::cout <<  "Finished!" << std::endl;
#endif

#ifdef METHOD2
	std::cout << "Starting ... " << std::endl;
	std::vector< MultidimArray<double > > RieszVector(3), RieszVector_original(3);
	MultidimArray<double> mask_amplitudeMS, mask_amplitudeMS_original;
	Image<double> imgVol, amplitudeMS, mask, outputResolution, amplitudeMS_original;
	double idx_threshold_energy, idx_threshold_noise, threshold_energy, threshold_noise, step = 0.01;
	size_t len_energy_vector, len_noise_vector, count_freq, N_elems = (size_t) (0.5/step);
	MultidimArray <double> filteredVol;

	//Reading input Volume
	imgVol.read(fnVol);
	std::cout << "Volume read!" << std::endl;

	size_t zdim = ZSIZE(imgVol())*0.5; //radius
	size_t ydim = YSIZE(imgVol())*0.5;
	size_t xdim = XSIZE(imgVol())*0.5;

	// Calculating Riesz transform

	int count = 0;
	double r2, r_kij;
	if (xdim<ydim){
		if (xdim<zdim){
			r2 = xdim*xdim;}
		else{
			r2 = zdim*zdim;}}
	else{
		if (ydim<zdim){
			r2 = ydim*ydim;}
		else{
			r2 = zdim*zdim;}
	}

	std::cout << "r2 = " << r2 << std::endl;

	outputResolution().initZeros(imgVol());
	mask.read(fnMask);

	MultidimArray<double> &p2imgVol = imgVol(), &p2amplitudeMS = amplitudeMS(), &p2mask = mask(), p2amplitudeMS_original = amplitudeMS_original();

	/////////////
	RieszTransform3Dreal(imgVol(), RieszVector_original);

	amplitudeMonogenicSignal3D(imgVol(), RieszVector_original, p2amplitudeMS_original);

	size_t sum_mask = 0;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(p2mask)
	{
		if (A3D_ELEM(p2mask, k, i, j) > 0)
		sum_mask++;
	}

	std::cout << "Analysing frequencies..." << std::endl;
	std::cout << "Sampling rate = " << smpr << std::endl;
	for (size_t idx_freq = 1; idx_freq<N_elems; idx_freq++)
	{
		if (count_freq == (sum_mask))
			break;

		count_freq = 0;
		std::vector<double> noise_vector(1), energy_vector(1), noise_vector_sorted, energy_vector_sorted;

		filteredVol.initZeros(p2amplitudeMS);
		if (fType == 0.0)
			passbandfiltervol(p2imgVol, idx_freq*step, filteredVol, 0.01);
		else
		{
			if (fType == 1.0)
			{
				highpassfiltervol(p2imgVol, idx_freq*step, filteredVol, 0.01);
			}
			else
			{
				std::cout << "Error: Type of filter was not provided" << std::endl;
				exit(0);
			}
		}

		RieszTransform3Dreal(filteredVol, RieszVector);

		amplitudeMonogenicSignal3D(filteredVol, RieszVector, p2amplitudeMS);

		applymask(p2amplitudeMS, p2mask, mask_amplitudeMS);

		FOR_ALL_ELEMENTS_IN_ARRAY3D(p2amplitudeMS)
		{
			if (A3D_ELEM(p2mask, k, i, j) > 0 )
			{
				energy_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
			}
			else
			{
				if (condpremask == true)
				{
					r_kij = (k-zdim)*(k-zdim) + (j-ydim)*(j-ydim) + (i-xdim)*(i-xdim);
					if (r_kij<r2){
						noise_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
						//noise_vector.push_back(A3D_ELEM(p2amplitudeMS_original, k, i, j));
					}
				}
				else{
					noise_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
					//noise_vector.push_back(A3D_ELEM(p2amplitudeMS_original, k, i, j));
				}
			}
		}


		std::sort(energy_vector.begin(), energy_vector.end());
		std::sort(noise_vector.begin(), noise_vector.end());


		//len_energy_vector = energy_vector.size();
		len_noise_vector = noise_vector.size();


		//idx_threshold_energy = floor((double) 0.95*len_energy_vector);
		idx_threshold_noise = floor((double) 0.95*len_noise_vector);
		std::cout << "len_noise = " << len_noise_vector << std::endl;


		//threshold_energy = energy_vector[idx_threshold_energy];
		threshold_noise = noise_vector[idx_threshold_noise];
		std::cout << "threshold_noise = " << threshold_noise << std::endl;

		double inv_idx_freq = (double) 1/(idx_freq*step);
		std::cout << "freq = " << idx_freq << std::endl;
		std::cout << "inversa freq = " << inv_idx_freq << std::endl;
		FOR_ALL_ELEMENTS_IN_ARRAY3D(outputResolution())
		{
			if (A3D_ELEM(mask_amplitudeMS, k, i, j) > threshold_noise)
			{
				A3D_ELEM(outputResolution(), k, i, j) = smpr*inv_idx_freq;
			}
			else
			{
				count_freq++;
			}
		}
	}

	outputResolution.write(fnDir);

	std::cout <<  "Finished!" << std::endl;
#endif


#ifdef METHOD4
	std::cout << "Starting ... " << std::endl;
	std::vector< MultidimArray<double > > RieszVector(3), RieszVector_original(3);
	MultidimArray<double> mask_amplitudeMS, mask_amplitudeMS_original;
	Image<double> imgVol, amplitudeMS, mask, outputResolution, amplitudeMS_original;
	double idx_threshold_energy, idx_threshold_noise, threshold_energy, threshold_noise, step = 0.01;
	size_t len_energy_vector, len_noise_vector, count_freq, N_elems = (size_t) (0.5/step);
	MultidimArray <double> filteredAmpMG;

	//Reading input Volume
	imgVol.read(fnVol);
	std::cout << "Volume read!" << std::endl;

	size_t zdim = ZSIZE(imgVol())*0.5; //radius
	size_t ydim = YSIZE(imgVol())*0.5;
	size_t xdim = XSIZE(imgVol())*0.5;

	// Calculating Riesz transform

	int count = 0;
	double r2, r_kij;
	if (xdim<ydim){
		if (xdim<zdim){
			r2 = xdim*xdim;}
		else{
			r2 = zdim*zdim;}}
	else{
		if (ydim<zdim){
			r2 = ydim*ydim;}
		else{
			r2 = zdim*zdim;}
	}

	std::cout << "r2 = " << r2 << std::endl;

	outputResolution().initZeros(imgVol());
	mask.read(fnMask);

	MultidimArray<double> &p2imgVol = imgVol(), &p2amplitudeMS = amplitudeMS(), &p2mask = mask(), &p2amplitudeMS_original = amplitudeMS_original();

	/////////////
	RieszTransform3Dreal(imgVol(), RieszVector_original);

	amplitudeMonogenicSignal3D(imgVol(), RieszVector_original, p2amplitudeMS_original);

	size_t sum_mask = 0;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(p2mask)
	{
		if (A3D_ELEM(p2mask, k, i, j) > 0)
		sum_mask++;
	}

	std::cout << "Analysing frequencies..." << std::endl;
	std::cout << "Sampling rate = " << smpr << std::endl;
	for (size_t idx_freq = 1; idx_freq<N_elems; idx_freq++)
	{
		if (count_freq == (sum_mask))
			break;

		count_freq = 0;
		std::vector<double> noise_vector(1), energy_vector(1), noise_vector_sorted, energy_vector_sorted;

		filteredAmpMG.initZeros(p2amplitudeMS);
		if (fType == 0.0)
			passbandfiltervol(p2amplitudeMS_original, idx_freq*step, filteredAmpMG, 0.01);
		else
		{
			if (fType == 1.0)
			{
				highpassfiltervol(p2amplitudeMS_original, idx_freq*step, filteredAmpMG, 0.01);
			}
			else
			{
				std::cout << "Error: Type of filter was not provided" << std::endl;
				exit(0);
			}
		}

		applymask(filteredAmpMG, p2mask, mask_amplitudeMS);

		FOR_ALL_ELEMENTS_IN_ARRAY3D(p2amplitudeMS)
		{
			if (A3D_ELEM(p2mask, k, i, j) > 0 )
			{
				energy_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
			}
			else
			{
				if (condpremask == true)
				{
					r_kij = (k-zdim)*(k-zdim) + (j-ydim)*(j-ydim) + (i-xdim)*(i-xdim);
					if (r_kij<r2){
						noise_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
						//noise_vector.push_back(A3D_ELEM(p2amplitudeMS_original, k, i, j));
					}
				}
				else{
					noise_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
					//noise_vector.push_back(A3D_ELEM(p2amplitudeMS_original, k, i, j));
				}
			}
		}


		std::sort(energy_vector.begin(), energy_vector.end());
		std::sort(noise_vector.begin(), noise_vector.end());


		//len_energy_vector = energy_vector.size();
		len_noise_vector = noise_vector.size();


		//idx_threshold_energy = floor((double) 0.95*len_energy_vector);
		idx_threshold_noise = floor((double) 0.95*len_noise_vector);
		std::cout << "len_noise = " << len_noise_vector << std::endl;


		//threshold_energy = energy_vector[idx_threshold_energy];
		threshold_noise = noise_vector[idx_threshold_noise];
		std::cout << "threshold_noise = " << threshold_noise << std::endl;

		double inv_idx_freq = (double) 1/(idx_freq*step);
		std::cout << "freq = " << idx_freq << std::endl;
		std::cout << "inversa freq = " << inv_idx_freq << std::endl;
		FOR_ALL_ELEMENTS_IN_ARRAY3D(outputResolution())
		{
			if (A3D_ELEM(mask_amplitudeMS, k, i, j) > threshold_noise)
			{
				A3D_ELEM(outputResolution(), k, i, j) = smpr*inv_idx_freq;
			}
			else
			{
				count_freq++;
			}
		}
	}

	outputResolution.write(fnDir);

	std::cout <<  "Finished!" << std::endl;
#endif

#ifdef TESTING
	std::cout << "Starting ... " << std::endl;
	std::vector< MultidimArray<double > > RieszVector(3), RieszVector_original(3);
	MultidimArray<double> mask_amplitudeMS, mask_amplitudeMS_original;
	Image<double> imgVol, amplitudeMS, mask, outputResolution, amplitudeMS_original, filteredVol;
	double idx_threshold_energy, idx_threshold_noise, threshold_energy, threshold_noise, max_val, step = 0.01;
	size_t len_energy_vector, len_noise_vector, count_freq, count_freq2, N_elems = (size_t) (0.5/step);
	//MultidimArray <double> filteredVol;
	FileName cmd;


	//Reading input Volume
	imgVol.read(fnVol);
	std::cout << "Volume read!" << std::endl;

	size_t zdim = ZSIZE(imgVol())*0.5; //radius
	size_t ydim = YSIZE(imgVol())*0.5;
	size_t xdim = XSIZE(imgVol())*0.5;


	// Calculating Riesz transform

	int count = 0;
	double r2, r_kij;
	if (xdim<ydim){
		if (xdim<zdim){
			r2 = xdim*xdim;}
		else{
			r2 = zdim*zdim;}}
	else{
		if (ydim<zdim){
			r2 = ydim*ydim;}
		else{
			r2 = zdim*zdim;}
	}

	outputResolution().initZeros(imgVol());
	mask.read(fnMask);

	MultidimArray<double> &p2imgVol = imgVol(), &p2amplitudeMS = amplitudeMS(), &p2mask = mask(), p2amplitudeMS_original = amplitudeMS_original();
	MultidimArray<double> outputResolution_aux, mask_filteredvol;
	MultidimArray<double> &p2filteredVol = filteredVol(), &p2mask_filteredvol = mask_filteredvol;
	MultidimArray<double> &p2outputResolution_aux = outputResolution_aux, &p2outputResolution = outputResolution();
	/////////////
	outputResolution_aux.initZeros(p2imgVol);

	RieszTransform3Dreal(imgVol(), RieszVector_original);

	amplitudeMonogenicSignal3D(imgVol(), RieszVector_original, p2amplitudeMS_original);

	max_val = p2amplitudeMS_original.computeMax();

	size_t sum_mask = 0;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(p2mask)
	{
		if (A3D_ELEM(p2mask, k, i, j) > 0)
		sum_mask++;
	}

	//max_val = imgVol().computeMax();
	std::cout << "max = " << max_val << std::endl;
	count_freq = sum_mask*0.5;
	std::cout << "Analysing frequencies..." << std::endl;
	std::cout << "Sampling rate = " << smpr << std::endl;
	for (size_t idx_freq = 1; idx_freq<N_elems; idx_freq++)
	{
		std::cout << "------------------------------------------------------" << std::endl;
//		std:: cout << "count_freq =  " << count_freq << std::endl;
//		std:: cout << "sum_mask*0.80 =  " << sum_mask*0.80 << std::endl;
//		if (count_freq > (sum_mask*0.80))
//		{
//			std::cout << "BREAK" << std::endl;
//			break;
//		}


		count_freq = 0;
		count_freq2 = 0;
		std::vector<double> noise_vector(1), energy_vector(1), noise_vector_sorted, energy_vector_sorted;

		if (fType == 0)
		{
		cmd = formatString("xmipp_transform_filter -i %s --fourier band_pass %f %f -o %s/filtered_vol.vol", fnVol.c_str(),
							((double) idx_freq*step)-0.002 , ((double) idx_freq*step)+0.002, fnDir.c_str());
		}
		else
		{
		cmd = formatString("xmipp_transform_filter -i %s --fourier high_pass %f -o %s/filtered_vol.vol", fnVol.c_str(),
							((double) idx_freq*step) , fnDir.c_str());
		}
		std::cout << cmd << std::endl;
		system(cmd.c_str());

		filteredVol.read(fnDir+"/filtered_vol.vol");

		//applymask(p2filteredVol, p2mask, p2mask_filteredvol);

		RieszTransform3Dreal(p2filteredVol, RieszVector);

		amplitudeMonogenicSignal3D(p2filteredVol, RieszVector, p2amplitudeMS);

		applymask(p2amplitudeMS, p2mask, mask_amplitudeMS);

		double max_fil = filteredVol().computeMax();
		std::cout << "max = " << max_fil << std::endl;

		FOR_ALL_ELEMENTS_IN_ARRAY3D(p2amplitudeMS)
		{
			if (A3D_ELEM(p2mask, k, i, j) > 0 )
			{
				energy_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
			}
			else
			{
				if (condpremask == true)
				{
					r_kij = (k-zdim)*(k-zdim) + (j-ydim)*(j-ydim) + (i-xdim)*(i-xdim);
					if (r_kij<r2){
						noise_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
						//noise_vector.push_back(A3D_ELEM(p2amplitudeMS_original, k, i, j));
					}
				}
				else{
					noise_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
					//noise_vector.push_back(A3D_ELEM(p2amplitudeMS_original, k, i, j));
				}
			}
		}


		std::sort(energy_vector.begin(), energy_vector.end());
		std::sort(noise_vector.begin(), noise_vector.end());


		len_energy_vector = energy_vector.size();
		len_noise_vector = noise_vector.size();


		idx_threshold_energy = floor((double) 0.95*len_energy_vector);
		idx_threshold_noise = floor((double) 0.99*len_noise_vector);
		//std::cout << "len_noise = " << len_noise_vector << std::endl;


		threshold_energy = energy_vector[idx_threshold_energy];
		threshold_noise = noise_vector[idx_threshold_noise];
		std::cout << "threshold_noise = " << threshold_noise << std::endl;
		std::cout << "threshold_energy = " << threshold_energy << std::endl;

		double inv_idx_freq = (double) 1/(idx_freq*step);
		std::cout << "freq = " << idx_freq << std::endl;
		std::cout << "inversa freq = " << inv_idx_freq << std::endl;

		FOR_ALL_ELEMENTS_IN_ARRAY3D(p2outputResolution_aux)
		{
			if (A3D_ELEM(mask_amplitudeMS, k, i, j) > threshold_noise)
			{
				A3D_ELEM(p2outputResolution_aux, k, i, j) = smpr*inv_idx_freq;
				count_freq2++;
			}
			else
			{
				count_freq++;
			}
		}

		//std::cout << "count_freq = " << count_freq << std::endl;
		//std::cout << "count_freq2 = " << count_freq2 << std::endl;


//		if (count_freq2>0)
//		{
//			std::cout << "opening... " << std::endl;
//			Image<double> aux_aux = p2outputResolution_aux;
//			aux_aux.write(fnDir+"/caca.vol");
//			opening3D(p2outputResolution_aux, p2outputResolution, 6, 0, 1);
//		}

	}

	Image<double> outputResolution2;

	outputResolution() = outputResolution_aux;
	medianFilter3x3x3(outputResolution(), outputResolution2());


	outputResolution2.write(fnDir+"/MGresolution.vol");

	std::cout <<  "Finished!" << std::endl;
#endif
}



//void ProgMonogenicSignalRes::run()
//{
//	std::cout << "Starting ... " << std::endl;
//
//	//Initializing variables
//	std::vector< MultidimArray<double > > RieszVector(3), RieszVector_original(3);
//	MultidimArray<double> mask_amplitudeMS, mask_amplitudeMS_original;
//	Image<double> imgVol, amplitudeMS, mask, outputResolution, amplitudeMS_original;
//	double idx_threshold_energy, idx_threshold_noise, threshold_energy, threshold_noise, r_kij, r2, step = 0.01;
//	size_t len_energy_vector, len_noise_vector, count_freq = 0, N_elems = (size_t) (0.5/step);
//	MultidimArray <double> filteredVol;
//	int count = 0;
//
//
//	//Reading input Volume
//	imgVol.read(fnVol);
//	std::cout << "Volume read!" << std::endl;
//
//	MultidimArray<double> &p2imgVol = imgVol();
//
//	//Volume dimensions
//	size_t zdim = ZSIZE(p2imgVol)*0.5; //radius
//	size_t ydim = YSIZE(p2imgVol)*0.5;
//	size_t xdim = XSIZE(p2imgVol)*0.5;
//
//	//Calculating Riesz transform
//	if (xdim<ydim){
//		if (xdim<zdim){
//			r2 = xdim*xdim;}
//		else{
//			r2 = zdim*zdim;}}
//	else{
//		if (ydim<zdim){
//			r2 = ydim*ydim;}
//		else{
//			r2 = zdim*zdim;}
//	}
//
//	outputResolution().initZeros(p2imgVol);
//	mask.read(fnMask);
//	MultidimArray<double> &p2outputResolution = outputResolution();
//	MultidimArray<double> &p2mask = mask();
//	MultidimArray<double> &p2amplitudeMS = amplitudeMS();
//	MultidimArray<double> &p2mask_amplitudeMS = mask_amplitudeMS;
//	MultidimArray<double> &p2amplitudeMS_original = amplitudeMS_original();
//
//	size_t sum_mask = 0;
//	FOR_ALL_ELEMENTS_IN_ARRAY3D(p2mask)
//	{
//		if (A3D_ELEM(p2mask, k, i, j) > 0)
//		sum_mask++;
//	}
//	std::cout << "len_mask = " << sum_mask << std::endl;
//
//	std::cout << "Analysing frequencies..." << std::endl;
//	std::cout << "Sampling rate = " << smpr << std::endl;
//	filteredVol.initZeros(p2imgVol);
//	MultidimArray<double> &p2filteredVol = filteredVol;
//	MultidimArray<double> &p2amplitudeMS_median = p2amplitudeMS;
//
////	for (size_t idx_freq = 1; idx_freq<N_elems; idx_freq++)
////	{
//	size_t idx_freq = 1;
//
////		if (count_freq > floor(sum_mask*0.30))
////			break;
//
//		count_freq = 0;
//		std::vector<double> noise_vector(1), energy_vector(1), noise_vector_sorted, energy_vector_sorted;
//
//		if (fType == 0.0)
//			passbandfiltervol(p2imgVol, idx_freq*step, p2filteredVol, 0.01);
//		else
//		{
//			if (fType == 1.0)
//			{
//				highpassfiltervol(p2imgVol, idx_freq*step, p2filteredVol, 0.01);
//			}
//			else
//			{
//				std::cout << "Error: Type of filter was not provided" << std::endl;
//				exit(0);
//			}
//		}
//
//		RieszTransform3Dreal(p2filteredVol, RieszVector);
//
//		amplitudeMonogenicSignal3D(p2filteredVol, RieszVector, p2amplitudeMS);
//
//
//		FileName fnMG = "MG_signal_pre.vol";
//		amplitudeMS().write(fnMG);
//
//		medianFilter3x3x3(p2amplitudeMS, p2amplitudeMS_median);
//
//		applymask(p2amplitudeMS_median, p2mask, p2mask_amplitudeMS);
//
//		double mean_energy=0;
//		FOR_ALL_ELEMENTS_IN_ARRAY3D(p2amplitudeMS_median)
//		{
//			if (A3D_ELEM(p2mask, k, i, j) > 0 )
//			{
//				energy_vector.push_back(A3D_ELEM(p2amplitudeMS_median, k, i, j));
//				mean_energy = A3D_ELEM(p2amplitudeMS_median, k, i, j) + mean_energy;
//			}
//			else
//			{
//				if (condpremask == true)
//				{
//					r_kij = (k-zdim)*(k-zdim) + (j-ydim)*(j-ydim) + (i-xdim)*(i-xdim);
//					if (r_kij<r2){
//						noise_vector.push_back(A3D_ELEM(p2amplitudeMS_median, k, i, j));
//						//noise_vector.push_back(A3D_ELEM(p2amplitudeMS_original, k, i, j));
//					}
//				}
//				else{
//					noise_vector.push_back(A3D_ELEM(p2amplitudeMS_median, k, i, j));
//					//noise_vector.push_back(A3D_ELEM(p2amplitudeMS_original, k, i, j));
//				}
//			}
//		}
//
//
//		std::sort(energy_vector.begin(), energy_vector.end());
//		std::sort(noise_vector.begin(), noise_vector.end());
//
//
//		len_energy_vector = energy_vector.size();
//		len_noise_vector = noise_vector.size();
//
//		mean_energy = mean_energy/len_energy_vector;
//		std::cout << "mean_energy = " << mean_energy << std::endl;
//
//		//idx_threshold_energy = floor((double) 0.95*len_energy_vector);
//		idx_threshold_noise = floor((double) 0.95*len_noise_vector);
//		std::cout << "len_noise = " << len_noise_vector << std::endl;
//
//
//		//threshold_energy = energy_vector[idx_threshold_energy];
//		threshold_noise = noise_vector[idx_threshold_noise];
//		std::cout << "threshold_noise = " << threshold_noise << std::endl;
//
//		if (threshold_noise<mean_energy)
//			std::cout << "resolution valida" << std::endl;
//
//		double inv_idx_freq = (double) 1/(idx_freq*step);
//		std::cout << "freq = " << idx_freq << std::endl;
//		std::cout << "resolution = " << inv_idx_freq << "A" << std::endl;
//		FOR_ALL_ELEMENTS_IN_ARRAY3D(p2outputResolution)
//		{
//			if (A3D_ELEM(p2mask_amplitudeMS, k, i, j) > threshold_noise)
//			{
//				A3D_ELEM(p2outputResolution, k, i, j) = smpr*inv_idx_freq;
//			}
//			else
//			{
//				count_freq++;
//			}
//		}
//		std::cout << "------------------------------------" << std::endl;
////	}
//
//	Image<double> outputResolution2;
//
//	medianFilter3x3x3(p2outputResolution, outputResolution2());
//
//	outputResolution.write(fnDir);
//
//	std::cout <<  "Finished!" << std::endl;
//}
