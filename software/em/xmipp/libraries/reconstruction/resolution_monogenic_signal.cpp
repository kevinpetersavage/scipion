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
//#define DEBUG

void ProgMonogenicSignalRes::readParams()
{
	fnVol = getParam("--vol");
	fnVol1 = getParam("--vol1");
	fnVol2 = getParam("--vol2");
	fnOut = getParam("-o");
	fnMask = getParam("--mask");
	fnchim = getParam("--chimera_volume");
	sampling = getDoubleParam("--sampling_rate");
	minRes = getDoubleParam("--minRes");
	maxRes = getDoubleParam("--maxRes");
	R = getIntParam("--circular_mask");
	N_freq = getDoubleParam("--number_frequencies");
	kValue = getDoubleParam("--trimmed");


}

void ProgMonogenicSignalRes::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("  --vol <vol_file>   : Input volume");
	addParamsLine("  --mask <vol_file>  : Mask defining the macromolecule");
	addParamsLine("  [--vol1 <vol_file=\"\">]: Half volume 1");
	addParamsLine("                          :+ If two half volume are given, the noise is estimated from them");
	addParamsLine("                          :+ Otherwise the noise is estimated outside the mask");
	addParamsLine("  [--vol2 <vol_file=\"\">]: Half volume 2");
	addParamsLine("  [-o <output=\"MGresolution.vol\">]: Local resolution volume (in Angstroms)");
	addParamsLine("  [--chimera_volume <output=\"Chimera_resolution_volume.vol\">]: Local resolution volume for chimera viewer (in Angstroms)");
	addParamsLine("  [--sampling_rate <s=1>]   : Sampling rate (A/px)");
	addParamsLine("  [--circular_mask <R=-1>]  : The volume has been masked to a sphere of this radius (in pixels)");
	addParamsLine("                            : Use -1 to disable this option");
	addParamsLine("  [--number_frequencies <w=50>]       : The resolution is computed at a number of frequencies between mininum and");
	addParamsLine("                            : maximum resolution px/A. This parameter determines that number");
	addParamsLine("  [--minRes <s=30>]         : Minimum resolution (A)");
	addParamsLine("  [--maxRes <s=1>]          : Maximum resolution (A)");
	addParamsLine("  [--trimmed <s=4>]         : Trimming value for smoothing the output resolution");

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
					DIRECT_MULTIDIM_ELEM(iu,n) = 1e38;
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

	if ((fnVol1 !="") && (fnVol2 !=""))
	{
		Image<double> V1, V2;
		V1.read(fnVol1);
		V2.read(fnVol2);

		V1()-=V2();
		V1()/=sqrt(2);

		fftN=new MultidimArray< std::complex<double> >;
		transformer.FourierTransform(inputVol, *fftN);
		halfMapsGiven = true;
	}
	else
	{
		fftN=&fftV;
		halfMapsGiven = false;
	}

	NS=0, NN=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pMask)
	{
		if (DIRECT_MULTIDIM_ELEM(pMask, n)==1)
		{
			NS++;
			if (halfMapsGiven)
				NN++;
		}
		else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
			NN++;
	}
}

void ProgMonogenicSignalRes::amplitudeMonogenicSignal3D(MultidimArray< std::complex<double> > &myfftV,
		double w1, MultidimArray<double> &amplitude, int count)
{
	fftVRiesz.initZeros(myfftV);
	amplitude.resizeNoCopy(VRiesz);


	// Filter the input volume and add it to amplitude
	long n=0;
	double iu1=1.0/w1;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				if (iun<=iu1)
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

	if (verbose>=2)
	{
	Image<double> saveImg2;
	saveImg2() = VRiesz;
	FileName fnSaveImg0;
	fnSaveImg0 = formatString("filteredVolume_%i.vol", count);
	saveImg2.write(fnSaveImg0);
	saveImg2.clear();
	}


	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	DIRECT_MULTIDIM_ELEM(amplitude,n)=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate first component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	double uz, uy, ux;
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				if (iun<=iu1)
				{
					FFT_IDX2DIGFREQ(j,XSIZE(amplitude),ux);
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (ux*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
				}
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate second component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			FFT_IDX2DIGFREQ(i,YSIZE(amplitude),uy);
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				if (iun<=iu1)
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (uy*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate third component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		FFT_IDX2DIGFREQ(k,ZSIZE(amplitude),uz);
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				if (iun<=iu1)
					DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = (uz*iun)*DIRECT_MULTIDIM_ELEM(myfftV, n);
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

	if (verbose>=2)
	{
	Image<double> saveImg;
	saveImg = amplitude;
	FileName fnSaveImg;
	fnSaveImg = formatString("amplitudeMS_%i.vol", count);
	saveImg.write(fnSaveImg);
	saveImg.clear();
	}

	// Low pass filter the monogenic amplitude
	lowPassFilter.w1 = w1;
	amplitude.setXmippOrigin();
	lowPassFilter.applyMaskSpace(amplitude);

	if (verbose>=2)
	{
	Image<double> saveImg2;
	saveImg2 = amplitude;
	FileName fnSaveImg2;
	fnSaveImg2 = formatString("amplitudeMfiltered_%i.vol", count);
	saveImg2.write(fnSaveImg2);
	saveImg2.clear();
//	char c;
//
//	std::cout << "Press any key ";
//	std::cin >> c;
	}

}

void ProgMonogenicSignalRes::medianFilter3x3x3Thresholding(MultidimArray<double> &m, MultidimArray<double> &out, double threshold)
{
	out = m;
	std::vector<double> values(6);

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

		  double val = 0.5* (values[2] + values[3]);
		  if (val>10*threshold)
			  val = 15*threshold;
		  if (val>threshold)
			  A3D_ELEM(out, k, i, j) = val;
		}
	  }
	}
}




void ProgMonogenicSignalRes::run()
{
	produceSideInfo();

	Image<double> outputResolution;
	outputResolution().initZeros(VRiesz);

	MultidimArray<int> &pMask = mask();
	MultidimArray<double> &pOutputResolution = outputResolution();
	MultidimArray<double> amplitudeMS, amplitudeMN;

	std::cout << "Looking for maximum frequency ..." << std::endl;
	double criticalZ=icdf_gauss(0.95);
	double criticalW=-1;
	double resolution, last_resolution = 10000;  //A huge value for achieving last_resolution < resolution
	double freq;
	double max_meanS = -1e38;

	double w0 = sampling/maxRes;
	double wF = sampling/minRes;
	double w=w0;
	double stepW = (wF-w0)/N_freq;
	bool doNextIteration=true;
	int iter=0;
	do
	{
		resolution =  sampling/w;
		std::cout << " Resolution = " << resolution << std::endl;
		freq = sampling/resolution;

		std::cout << "------------------------------------------------------" << std::endl;
		std::cout << "Iteration " << iter << " Freq = " << freq << " Resolution = " << resolution << " (A)" << std::endl;


		amplitudeMonogenicSignal3D(fftV, freq, amplitudeMS, iter);
		if (halfMapsGiven)
			amplitudeMonogenicSignal3D(*fftN, freq, amplitudeMN, iter);

		double sumS=0, sumS2=0, sumN=0, sumN2=0;
		if (halfMapsGiven)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
			{
				double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
				double amplitudeValueN=DIRECT_MULTIDIM_ELEM(amplitudeMN, n);
				if (DIRECT_MULTIDIM_ELEM(pMask, n)==1)
				{
					sumS  += amplitudeValue;
					sumS2 += amplitudeValue*amplitudeValue;
					sumN  += amplitudeValueN;
					sumN2 += amplitudeValueN*amplitudeValueN;
				}
				else if (DIRECT_MULTIDIM_ELEM(pMask, n)==0)
				{
					sumN  += amplitudeValueN;
					sumN2 += amplitudeValueN*amplitudeValueN;
				}
			}
		}
		else
		{
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
		}
		double meanS=sumS/NS;
		double sigma2S=sumS2/NS-meanS*meanS;
		double meanN=sumN/NN;
		double sigma2N=sumN2/NN-meanN*meanN;

		if (meanS>max_meanS)
			max_meanS = meanS;

		if (meanS<0.001*max_meanS)
		{
			std::cout << "Search of resolutions stopped due to too low signal" << std::endl;
			break;
		}

		// Check local resolution
		double thresholdNoise=meanN+criticalZ*sqrt(sigma2N);
		std::cout << "thresholdNoise = " << thresholdNoise << std::endl;
		std::cout << "uncertain = " << criticalZ*sqrt(sigma2N) << std::endl;
		int Entro = 0;

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
		{
			if (DIRECT_MULTIDIM_ELEM(pMask, n)==1)
				if (DIRECT_MULTIDIM_ELEM(amplitudeMS, n)>thresholdNoise)
				{
					DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = freq;
					++Entro;
				}
		}

		std::cout << "mascara = " << Entro << std::endl;
		// Continue?
		// Is the mean inside the signal significantly different from the noise?
		double z=(meanS-meanN)/sqrt(sigma2S/NS+sigma2N/NN);
		std::cout << "  meanS= " << meanS << " sigma2S= " << sigma2S << " NS= " << NS << std::endl;
		std::cout << "  meanN= " << meanN << " sigma2N= " << sigma2N << " NN= " << NN << std::endl;
		std::cout << "  z=" << z << " (" << criticalZ << ")" << std::endl;

		if (z<criticalZ)
		{
			criticalW = freq;
			doNextIteration=false;
		}

		if (doNextIteration)
		{
			last_resolution = resolution;
			w+=stepW;
			resolution = sampling/w;
			if (last_resolution-resolution<0.1)
			{
				resolution=last_resolution-0.1;
				w = sampling/resolution;
			}
			if (w > wF)
			{
				doNextIteration = false;
				std::cout << "Search of resolutions stopped due to out of resolution range" << std::endl;
			}

		}
		iter++;
	} while (doNextIteration);




	double resValue;
	if (kValue == 0)
	{
		// Compute resolution in Angstroms
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
		{
			if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n)>0)
				DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = sampling/DIRECT_MULTIDIM_ELEM(pOutputResolution, n);
		}


		std::cout << "Trimming is not performed" << std::endl;
		outputResolution.write(fnOut);

		//Calculating means and stds
		double SumRes = 0, Nsum = 0, SumRes2 = 0;

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
		{
			resValue = DIRECT_MULTIDIM_ELEM(pOutputResolution, n);
			if (resValue>0)
			{
				SumRes  += resValue;
				SumRes2 += resValue*resValue;
				Nsum += 1;
			}
		}

		double meanRes = SumRes/Nsum;
		double sigmaRes = sqrt(SumRes2/Nsum-meanRes*meanRes);

		std::cout << "mean = " << meanRes << "  sigmaRes = " << sigmaRes << std::endl;

		if (fnchim != "")
		{
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
		{
			resValue = DIRECT_MULTIDIM_ELEM(pOutputResolution, n);
			if (resValue<=0)
				DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = meanRes;
		}

		outputResolution.write(fnchim);
		}
	}
	else{

		// Compute resolution in Angstroms
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
		{
			if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n)>0)
				DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = sampling/DIRECT_MULTIDIM_ELEM(pOutputResolution, n);
		}


	//Trimming volume
	double SumRes = 0, Nsum = 0, SumRes2 = 0;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
	{
		resValue = DIRECT_MULTIDIM_ELEM(pOutputResolution, n);
		if (resValue>0)
		{
			SumRes  += resValue;
			SumRes2 += resValue*resValue;
			Nsum += 1;
		}
	}

	double meanRes = SumRes/Nsum;
	double sigmaRes = sqrt(SumRes2/Nsum-meanRes*meanRes);

	std::cout << "mean = " << meanRes << "  sigmaRes = " << sigmaRes << std::endl;

	double maxValue = meanRes + kValue*sigmaRes;
	double minValue = meanRes - kValue*sigmaRes;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
	{
		resValue = DIRECT_MULTIDIM_ELEM(pOutputResolution, n);
		if (resValue>0)
		{
			if (resValue>maxValue)
				DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = maxValue;
			if (resValue<minValue)
				DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = minValue;
		}
	}

	outputResolution.write(fnOut);




	if (fnchim != "")
	{
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
	{
		resValue = DIRECT_MULTIDIM_ELEM(pOutputResolution, n);
		if (resValue>0)
		{
			if (resValue>maxValue)
				DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = maxValue;
			if (resValue<minValue)
				DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = minValue;
		}
		else
			DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = meanRes;
	}

	outputResolution.write(fnchim);
	}
	}


}
