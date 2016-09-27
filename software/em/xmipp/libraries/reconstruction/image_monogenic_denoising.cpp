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

#include "image_monogenic_denoising.h"
#define DEBUG

void ProgMonogenicDenoising::readParams()
{
	fnStack = getParam("--stack");
	fnOut = getParam("-o");
	sampling = getDoubleParam("--sampling_rate");
	minRes = getDoubleParam("--minRes");
	maxRes = getDoubleParam("--maxRes");
	R = getIntParam("--circular_mask");
	N_freq = getDoubleParam("--number_frequencies");
	linearchk = checkParam("--linear");
	fnSpatial = getParam("--filtered");

}


void ProgMonogenicDenoising::defineParams()
{
	addUsageLine("This function determines the local resolution of a map");
	addParamsLine("  --stack <md_file>   : Metadata Stack with particles");
	addParamsLine("  [-o <output=\"MGresolution.stk\">]: Local resolution volume (in Angstroms)");
	addParamsLine("  [--sampling_rate <s=1>]   : Sampling rate (A/px)");
	addParamsLine("  [--circular_mask <R=-1>]  : The volume has been masked to a sphere of this radius (in pixels)");
	addParamsLine("                            : Use -1 to disable this option");
	addParamsLine("  [--number_frequencies <w=50>]       : The resolution is computed at a number of frequencies between minimum and");
	addParamsLine("                            : maximum resolution px/A. This parameter determines that number");
	addParamsLine("  [--minRes <s=30>]         : Minimum resolution (A)");
	addParamsLine("  [--maxRes <s=1>]          : Maximum resolution (A)");
	addParamsLine("  [--linear]                : The search for resolution is linear (equidistance between resolutions).");
	addParamsLine("  [--filtered <stk_file=\"\">]: Local resolution volume (in Angstroms)");

}

void ProgMonogenicDenoising::produceSideInfo(Image<double> &V)
{

	FourierTransformer transformer;
	MultidimArray<double> &inputVol = V();
	VRiesz.resizeNoCopy(inputVol);
	if (fnSpatial!="")
		VresolutionFiltered().initZeros(V());

	transformer.FourierTransform(inputVol, fftV);
	iu.initZeros(fftV);

	// Calculate u and first component of Riesz vector
	double uz, uy, ux, u2, uz2y2;
	long n=0;
	for(size_t i=0; i<YSIZE(fftV); ++i)
	{
		FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
		uz2y2=uy*uy;

		for(size_t j=0; j<XSIZE(fftV); ++j)
		{
			FFT_IDX2DIGFREQ(j,XSIZE(inputVol),ux);
			u2=uz2y2+ux*ux;
			if ((i != 0) || (j != 0))
				DIRECT_MULTIDIM_ELEM(iu,n) = 1.0/sqrt(u2);
			else
				DIRECT_MULTIDIM_ELEM(iu,n) = 1e38;
			++n;
		}
	}
	V.clear();

	// Prepare low pass filter
	lowPassFilter.FilterShape = RAISED_COSINE;
	lowPassFilter.raised_w = 0.01;
	lowPassFilter.do_generate_3dmask = false;
	lowPassFilter.FilterBand = LOWPASS;

	fftN=&fftV;

}

void ProgMonogenicDenoising::amplitudeMonogenicSignal3D(MultidimArray< std::complex<double> > &myfftV,
		double w1, MultidimArray<double> &amplitude, int count)
{
	fftVRiesz.initZeros(myfftV);
	amplitude.resizeNoCopy(VRiesz);


	// Filter the input volume and add it to amplitude
	long n=0;
	double iu1=1.0/w1;

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

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	if (fnSpatial!="")
		Vfiltered()=VRiesz;

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

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate second component of Riesz vector
	fftVRiesz.initZeros(myfftV);
	n=0;

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

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate third component of Riesz vector
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
	fnSaveImg2 = formatString("amplitudeMfiltered_%i.xmp", count);
	saveImg2.write(fnSaveImg2);
	saveImg2.clear();
//	char c;
//
//	std::cout << "Press any key ";
//	std::cin >> c;
	}
}

void ProgMonogenicDenoising::computingResolutions(const MultidimArray<double> &p2mask, Image<double> &outputResolution, FileName fnVol, Image<double> &out_filtered)
{

	outputResolution().initZeros(VRiesz);

	MultidimArray<double> &pOutputResolution = outputResolution();
	MultidimArray<double> &pVfiltered = Vfiltered();
	MultidimArray<double> amplitudeMS, amplitudeMN;
	MultidimArray<double> &pVresolutionFiltered = VresolutionFiltered();

	std::cout << "Looking for maximum frequency ..." << std::endl;
	double significance = 0.95;
	double criticalZ=icdf_gauss(significance);
	double criticalW=-1;
	double resolution, last_resolution = 10000;  //A huge value for achieving last_resolution < resolution
	double freq;
	double max_meanS = -1e38;

	double range = maxRes-minRes;
	double R_ = range/N_freq;

	double w0 = sampling/maxRes;
	double wF = sampling/minRes;
	double w=w0;
	double stepW = (wF-w0)/N_freq;
	bool doNextIteration=true;
	int iter=0;
	int count_res = 0;

	std::vector<double> noiseValues;
	do
	{
		if (linearchk ==true)
		{
			resolution = maxRes - count_res*R_;
			freq = sampling/resolution;
			++count_res;
		}
		else
		{
			resolution =  sampling/w;
			freq = sampling/resolution;
		}

		std::cout << "Iteration " << iter << " Freq = " << freq << " Resolution = " << resolution << " (A)" << std::endl;


		amplitudeMonogenicSignal3D(fftV, freq, amplitudeMS, iter);

#ifdef DEBUG
//		Image<double> save;
//		save() = amplitudeMS;
//		save.write("PPP_amplitude.xmp");

		Image<double> saveImg2;
		saveImg2 = amplitudeMS;
		FileName fnSaveImg2;
		fnSaveImg2 = formatString("amplitude_%i.xmp", iter);
		saveImg2.write(fnSaveImg2);
#endif

		double sumS=0, sumS2=0, sumN=0, sumN2=0, NN = 0, NS = 0;
		noiseValues.clear();

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
		{
			double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
			if (DIRECT_MULTIDIM_ELEM(p2mask, n)>0)//==1)
			{
				sumS  += amplitudeValue;
				sumS2 += amplitudeValue*amplitudeValue;
				++NS;
			}
			else //if (DIRECT_MULTIDIM_ELEM(p2mask, n)==0)
			{
				sumN  += amplitudeValue;
				sumN2 += amplitudeValue*amplitudeValue;
				++NN;
			}
		}


		if (NS == 0)
		{
			std::cout << "There are no points to compute inside the mask" << std::endl;
			break;
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
		double thresholdNoise;

		thresholdNoise = meanN+criticalZ*sqrt(sigma2N);

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
		{
			if (DIRECT_MULTIDIM_ELEM(p2mask, n)>0)//==1)
			{
				if (DIRECT_MULTIDIM_ELEM(amplitudeMS, n)>thresholdNoise)
				{
					DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = freq;
					if (fnSpatial!="")
						DIRECT_MULTIDIM_ELEM(pVresolutionFiltered,n)=DIRECT_MULTIDIM_ELEM(pVfiltered,n);
				}
			}
			else
			{
				DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = 0;
			}
		}


		// Is the mean inside the signal significantly different from the noise?
		double z=(meanS-meanN)/sqrt(sigma2S/NS+sigma2N/NN);
		if (verbose >=2)
		{
			std::cout << "thresholdNoise = " << thresholdNoise << std::endl;
			std::cout << "  meanS= " << meanS << " sigma2S= " << sigma2S << " NS= " << NS << std::endl;
			std::cout << "  meanN= " << meanN << " sigma2N= " << sigma2N << " NN= " << NN << std::endl;
			std::cout << "  z=" << z << " (" << criticalZ << ")" << std::endl;
		}
		if (z<criticalZ)
		{
			criticalW = freq;
			doNextIteration=false;
		}

		if (doNextIteration)
		{
			if (linearchk == false)
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
			else
			{
				if (resolution == minRes)
					doNextIteration = false;
			}
		}
		iter++;
	} while (doNextIteration);

	if (fnSpatial!="")
	{
//	mask.read(fnMask);
//	mask().setXmippOrigin();
	Vfiltered.read(fnVol);
	pVfiltered=Vfiltered();
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pVfiltered)
	if (DIRECT_MULTIDIM_ELEM(p2mask,n)>0)
		DIRECT_MULTIDIM_ELEM(pVfiltered,n)-=DIRECT_MULTIDIM_ELEM(pVresolutionFiltered,n);
//		else
//			DIRECT_MULTIDIM_ELEM(pVfiltered,n)=0;
	//Vfiltered.write(fnSpatial);

	out_filtered = Vfiltered;

	VresolutionFiltered().clear();
	Vfiltered().clear();
	}


	// Compute resolution in Angstroms
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pOutputResolution)
	{
		if (DIRECT_MULTIDIM_ELEM(pOutputResolution, n)>0)
			DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = sampling/DIRECT_MULTIDIM_ELEM(pOutputResolution, n);
		else
			DIRECT_MULTIDIM_ELEM(pOutputResolution, n) = 1;
	}
	//outputResolution.write(fnOut);
	std::cout << "last_resolutoin = " << resolution << std::endl;
}


void ProgMonogenicDenoising::run()
{
	MetaData mdStack;
	FileName filename_particle;
	Image<double> Img_par;
	size_t Xdim, Ydim, Zdim, Ndim;

	mdStack.read(fnStack);
	mdStack.getValue(MDL_IMAGE, filename_particle, 1);
	Img_par.read(filename_particle);
	Img_par.getDimensions(Xdim, Ydim, Zdim, Ndim);

	//Defining the circular mask

	mask().initZeros((int) Xdim, (int) Ydim);
	mask().setXmippOrigin();

	MultidimArray<double> &p2mask=mask();

	R = (int) Xdim*0.5;
	double R2=R*R;
	FOR_ALL_ELEMENTS_IN_ARRAY2D(p2mask)
	{
		double r2=i*i+j*j;
		if (r2>=R2)
			A2D_ELEM(p2mask,i,j)=0;
		else
			A2D_ELEM(p2mask,i,j)=1;
	}

	Image<double> saveImg2;
	saveImg2() = p2mask;
	saveImg2.write("mask_.xmp");


	Image<double> outputResolution;

	FileName fnout_set = formatString("newset.stk");
	//FileName fnout_filtered = formatString("filtered_set.stk");
	size_t id_part;
	FOR_ALL_OBJECTS_IN_METADATA(mdStack)
	{
		mdStack.getValue(MDL_IMAGE, filename_particle, __iter.objId);
		mdStack.getValue(MDL_PARTICLE_ID, id_part, __iter.objId);
		std::cout << "particula = " << filename_particle << std::endl;
		std::cout << "id_part " << id_part << std::endl;


		Img_par.read(filename_particle);
		Img_par().setXmippOrigin();

		produceSideInfo(Img_par);

		Image<double> out_filtered;
		computingResolutions(p2mask, outputResolution, filename_particle, out_filtered);

		#ifdef DEBUG
		FileName fnout__aux;
		fnout__aux = formatString("amplitude_%i.xmp", id_part);
		outputResolution.write(fnout__aux);
		#endif

//		char c;
//
//		std::cout << "Press any key ";
//		std::cin >> c;
		outputResolution.write(fnout_set, id_part, true, WRITE_APPEND);
		if (fnSpatial!="")
			out_filtered.write(fnSpatial, id_part, true, WRITE_APPEND);
		out_filtered.clear();


	}
}
