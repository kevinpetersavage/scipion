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

void ProgMonogenicSignalRes::readParams()
{
	fnVol = getParam("--vol");
	fnOut = getParam("-o");
	fnMask = getParam("--mask");
	sampling = getDoubleParam("--sampling_rate");
	minRes = getDoubleParam("--minRes");
	maxRes = getDoubleParam("--maxRes");
	R = getIntParam("--circular_mask");
	stepW = getDoubleParam("--stepW");
}

void ProgMonogenicSignalRes::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("   --vol <vol_file>   : Input volume");
	addParamsLine("  [-o <output=\"MGresolution.vol\">]: Local resolution volume (in Angstroms)");
	addParamsLine("   --mask <vol_file=\"\">   : Mask defining the macromolecule");
	addParamsLine("  [--sampling_rate <s=1>]   : Sampling rate (A/px)");
	addParamsLine("  [--circular_mask <R=-1>]  : The volume has been masked to a sphere of this radius (in pixels)");
	addParamsLine("                            : Use -1 to disable this option");
	addParamsLine("  [--stepW <w=0.005>]       : Step in digital frequency (normalized to 0.5)");
	addParamsLine("  [--minRes <s=30>]         : Minimum resolution (A)");
	addParamsLine("  [--maxRes <s=1>]          : Maximum resolution (A)");
}

void ProgMonogenicSignalRes::produceSideInfo()
{
    Image<double> V; // Input volume
	V.read(fnVol);
	V().setXmippOrigin();

	FourierTransformer transformer;
	MultidimArray<double> &inputVol = V();
	VRiesz.resizeNoCopy(inputVol);

	transformer.FourierTransform(inputVol, fftV);
	iu.initZeros(fftV);

	// Calculate u and first component of Riesz vector
	double uz, uy, ux, uz2, u2, uz2y2;
	long n=0;
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
				if ((k != 0) || (i != 0) || (j != 0))
					DIRECT_MULTIDIM_ELEM(iu,n) = 1.0/sqrt(u2);
				else
					DIRECT_MULTIDIM_ELEM(iu,n) = 1.0;
				++n;
			}
		}
	}
	V.clear();

	// Prepare low pass filter
	lowPassFilter.FilterShape = RAISED_COSINE;
	lowPassFilter.raised_w = 0.01;
	lowPassFilter.do_generate_3dmask = false;
	lowPassFilter.FilterBand = LOWPASS;

	// Prepare mask
	mask.read(fnMask);
	mask().setXmippOrigin();
	MultidimArray<int> &pMask=mask();
	if (R>0)
	{
		double R2=R*R;
		FOR_ALL_ELEMENTS_IN_ARRAY3D(pMask)
		{
			double r2=k*k+i*i+j*j;
			if (r2>=R2)
				A3D_ELEM(pMask,k,i,j)=-1;
		}
	}

	NS=0, NN=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pMask)
	{
		if (DIRECT_MULTIDIM_ELEM(pMask, n)==1)
			NS++;
		else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
			NN++;
	}
}

void ProgMonogenicSignalRes::amplitudeMonogenicSignal3D(double w1, MultidimArray<double> &amplitude)
{
	fftVRiesz.initZeros(fftV);
	amplitude.resizeNoCopy(VRiesz);


	// Filter the input volume and add it to amplitude
	long n=0;
	double iu1=1.0/w1;
	for(size_t k=0; k<ZSIZE(fftV); ++k)
	{
		for(size_t i=0; i<YSIZE(fftV); ++i)
		{
			for(size_t j=0; j<XSIZE(fftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				if (iun<=iu1)
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(fftV, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	DIRECT_MULTIDIM_ELEM(amplitude,n)=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate first component of Riesz vector
	fftVRiesz.initZeros(fftV);
	double uz, uy, ux;
	n=0;
	for(size_t k=0; k<ZSIZE(fftV); ++k)
	{
		for(size_t i=0; i<YSIZE(fftV); ++i)
		{
			for(size_t j=0; j<XSIZE(fftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				if (iun<=iu1)
				{
					FFT_IDX2DIGFREQ(j,XSIZE(amplitude),ux);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (ux*iun)*DIRECT_MULTIDIM_ELEM(fftV, n);
				}
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate second component of Riesz vector
	fftVRiesz.initZeros(fftV);
	n=0;
	for(size_t k=0; k<ZSIZE(fftV); ++k)
	{
		for(size_t i=0; i<YSIZE(fftV); ++i)
		{
			FFT_IDX2DIGFREQ(i,YSIZE(amplitude),uy);
			for(size_t j=0; j<XSIZE(fftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				if (iun<=iu1)
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (uy*iun)*DIRECT_MULTIDIM_ELEM(fftV, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate third component of Riesz vector
	fftVRiesz.initZeros(fftV);
	n=0;
	for(size_t k=0; k<ZSIZE(fftV); ++k)
	{
		FFT_IDX2DIGFREQ(k,ZSIZE(amplitude),uz);
		for(size_t i=0; i<YSIZE(fftV); ++i)
		{
			for(size_t j=0; j<XSIZE(fftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				if (iun<=iu1)
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (uz*iun)*DIRECT_MULTIDIM_ELEM(fftV, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	{
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
		DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));
	}

	// Low pass filter the monogenic amplitude
	lowPassFilter.w1 = w1;
	amplitude.setXmippOrigin();
	lowPassFilter.applyMaskSpace(amplitude);
}

void ProgMonogenicSignalRes::medianFilter3x3x3Thresholding(MultidimArray<double> &m, MultidimArray<double> &out, double threshold)
{
	out = m;
	std::vector<double> values(6);
	double val;



	for (size_t k = 1; k<(ZSIZE(m)-1); k++)
	{
	  for (size_t i= 1; i<(YSIZE(m)-1); i++)
	  {
		for (size_t j= 1; j<(XSIZE(m)-1); j++)
		{
		  values[0] = A3D_ELEM(m, k-1, i, j);
		  values[1] = A3D_ELEM(m, k+1, i, j);
		  values[2] = A3D_ELEM(m, k, i-1, j);
		  values[3] = A3D_ELEM(m, k, i+1, j);
		  values[4] = A3D_ELEM(m, k, i, j-1);
		  values[5] = A3D_ELEM(m, k, i, j+1);

		  std::sort(values.begin(), values.end());

		  val = 0.5* (values[2] + values[3]);
		  if (val>10*threshold)
			  val = 15*threshold;
		  if (val>threshold)
			  A3D_ELEM(out, k, i, j) = val;

		}
	  }
	  std::cout << values[0] << " " << values[1] << " " << values[2] << " " << values[3] << " " << values[4] << " " << values[5] << std::endl;
	}


}


void ProgMonogenicSignalRes::run()
{
	produceSideInfo();

	Image<double> outputResolution;
	outputResolution().initZeros(VRiesz);

	MultidimArray<int> &pMask = mask();
	MultidimArray<double> &pOutputResolution = outputResolution();
	MultidimArray<double> amplitudeMS;

	std::cout << "Looking for maximum frequency ..." << std::endl;
	double criticalZ=icdf_gauss(0.95);
	double criticalW=-1;
	double initFreq = sampling/minRes;
	double endFreq = sampling/maxRes;
	double resolution, last_resolution = 10000;  //A huge value for achieving last_resolution < resolution
	double freq;

	std::cout << "Sampling Rate = " << sampling << std::endl;
	std::cout << "Initial Frequency = " << initFreq << std::endl;
	std::cout << "Final Frequency = " << endFreq << std::endl;
	for (double w=endFreq; w<=initFreq; w+=stepW)
	{
		resolution =  sampling/w;
		std::cout << " Resolution = " << resolution << std::endl;

		if (fabs(resolution-last_resolution)<0.1)
			resolution = floor(10*last_resolution)*0.1 - 0.1;
		if (w > endFreq)
		{
			last_resolution = resolution;
			freq = sampling/last_resolution;
		}
		else
			freq = sampling/resolution;


		std::cout << " last_resolution = " << resolution << std::endl;

		if (last_resolution < resolution)
			break;

		std::cout << "------------------------------------------------------" << std::endl;
		std::cout << "Freq = " << freq << " Resolution = " << resolution << " (A)" << std::endl;

		amplitudeMonogenicSignal3D(freq, amplitudeMS);

		double sumS=0, sumS2=0, sumN=0, sumN2=0;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
		{
			double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
			if (DIRECT_MULTIDIM_ELEM(pMask, n)==1)
			{
				sumS  += amplitudeValue;
				sumS2 += amplitudeValue*amplitudeValue;
			}
			else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
			{
				sumN  += amplitudeValue;
				sumN2 += amplitudeValue*amplitudeValue;
			}
		}
		double meanS=sumS/NS;
		double sigma2S=sumS2/NS-meanS*meanS;
		double meanN=sumN/NN;
		double sigma2N=sumN2/NN-meanN*meanN;

		// Check local resolution
		double thresholdNoise=meanN+criticalZ*sqrt(sigma2N);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
		{
			if (DIRECT_MULTIDIM_ELEM(pMask, n)==1)
				if (DIRECT_MULTIDIM_ELEM(amplitudeMS, n)>thresholdNoise)
					DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = freq;
		}

		// Continue?
		// Is the mean inside the signal significantly different from the noise?
		double z=(meanS-meanN)/sqrt(sigma2S/NS+sigma2N/NN);
		std::cout << "  meanS= " << meanS << " sigma2S= " << sigma2S << " NS= " << NS << std::endl;
		std::cout << "  meanN= " << meanN << " sigma2N= " << sigma2N << " NN= " << NN << std::endl;
		std::cout << "  z=" << z << " (" << criticalZ << ")" << std::endl;
		if (z<criticalZ)
		{
			criticalW = freq;
			break;
		}
	}

//	iu.clear();
//    mask.clear();
//    VRiesz.clear();
//	fftV.clear(); // Fourier transform of the input volume
//	transformer_inv.clear();
//	fftVRiesz.clear();



	// Compute resolution in Angstroms
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
	{
		if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n)>0)
			DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = sampling/DIRECT_MULTIDIM_ELEM(pOutputResolution, n);
	}
	outputResolution.write("MG_premedian_filter_Resolution.vol");


	Image<double> outputResolutionfiltered;


	std::cout << "Applying median filter" << std::endl;
	medianFilter3x3x3Thresholding(outputResolution(), outputResolutionfiltered(), sampling*freq);

	std::cout << "Saving" << std::endl;
	outputResolutionfiltered.write(fnOut);





}
