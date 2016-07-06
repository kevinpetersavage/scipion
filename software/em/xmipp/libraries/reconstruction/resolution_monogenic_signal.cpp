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
#include "fourier_filter.h"
#include <vector>
#include <math.h>
#include <limits>
#include <complex>
//#define DEBUGGING
#define FILTERING_TEST
//#define TESTING2



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

void ProgMonogenicSignalRes::applymask(const MultidimArray<double> &amplitudeMS, const MultidimArray<double> &mask,
																					  MultidimArray<double> &mask_amplitudeMS)
{
	mask_amplitudeMS.initZeros(amplitudeMS);
	FOR_ALL_ELEMENTS_IN_ARRAY3D(amplitudeMS)
	{
		A3D_ELEM(mask_amplitudeMS, k, i, j) = A3D_ELEM(amplitudeMS, k, i, j)*A3D_ELEM(mask, k, i, j);
	}
}

void ProgMonogenicSignalRes::getHalfDimensions(const  MultidimArray<double> &Vol, size_t &xdim, size_t &ydim, size_t &zdim)
{
	zdim = ZSIZE(Vol)*0.5; //radius
	ydim = YSIZE(Vol)*0.5;
	xdim = XSIZE(Vol)*0.5;
}

double ProgMonogenicSignalRes::calculateMaskRadius(const MultidimArray<double> &Vol, size_t &xdim, size_t &ydim, size_t &zdim)
{
	double r2;
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
	return r2;
}

void ProgMonogenicSignalRes::run()
{
#ifdef FILTERING_TEST
	std::cout << "Starting ... " << std::endl;
	std::vector< MultidimArray<double > > RieszVector(3), RieszVector_original(3);
	MultidimArray<double> mask_amplitudeMS, mask_amplitudeMS_original;
	Image<double> imgVol, amplitudeMS, mask, outputResolution, amplitudeMS_original, filteredVol, filteredMS, MGVol;
	double idx_threshold_energy, idx_freq_threshold = 0, idx_threshold_noise, threshold_energy, threshold_noise, max_val, r2, r_kij;
	size_t xdim, ydim, zdim, len_energy_vector, len_noise_vector, count_freq, count_freq2;
	FileName cmd;

	MultidimArray<double> &p2imgVol = imgVol();
	MultidimArray<double> &p2mask = mask();
	MultidimArray<double> &p2amplitudeMS = amplitudeMS();
	MultidimArray<double> mask_filteredvol;
	MultidimArray<double> &p2mask_filteredvol = mask_filteredvol;

	//Parameters:
	double Nsteps = 50;
	double step = 0.01;
	size_t N_elems = (size_t) (0.5/step);

	std::cout << "Nelems = " << N_elems << std::endl;

	//Reading input Volume
	imgVol.read(fnVol);
	imgVol().setXmippOrigin();

	//Calculating Fourier Transform
	FourierTransformer transformer, transformer2, transformer3, transformer4, transformer5;
	MultidimArray< std::complex<double> > fftV, fftAmpMS;
	transformer.FourierTransform(p2imgVol, fftV, false);

	// Reserve memory for filtered volume
	filteredVol() = imgVol();
	MGVol() = imgVol();
	transformer2.setReal(filteredVol());

	FourierFilter filter, filter2, filter3;
	filter.FilterShape = RAISED_COSINE;
	filter.raised_w = 0;//0.01;
	filter.do_generate_3dmask = false;

	//Obtaining dimensions
	getHalfDimensions(p2imgVol, xdim, ydim, zdim);

	if (condpremask)
		r2 = calculateMaskRadius(p2imgVol, xdim, ydim, zdim);

	//Reading mask
	mask.read(fnMask);
	p2mask.setXmippOrigin();

	std::cout << "Analysing frequencies..." << std::endl;
	double num_energy_mean, num_noise_mean, max_energy_mean=0;
	size_t count, count2;
	std::vector<double> noise_vector, energy_vector;

	for (size_t idx_freq = 1; idx_freq<N_elems; idx_freq++)
	{
		std::cout << "------------------------------------------------------" << std::endl;
		num_energy_mean = 0;
		num_noise_mean = 0;
		count = 0;
		count2 = 0;
		double inv_idx_freq = (double) 1.0/(idx_freq*step);
		std::cout << "freq = " << idx_freq << std::endl;
		std::cout << "Resolution = " << inv_idx_freq << std::endl;

		noise_vector.clear();
		energy_vector.clear();

		if (fType == 0)
		{
			exit(0);
		} else {
			filter.FilterBand = HIGHPASS;
			filter.w1 = idx_freq*step;
			transformer2.fFourier = transformer.fFourier;
		    filter.applyMaskFourierSpace(p2imgVol,transformer2.fFourier);
			transformer2.inverseFourierTransform();
		}

		RieszTransform3Dreal(filteredVol(), RieszVector);
		amplitudeMonogenicSignal3D(filteredVol(), RieszVector, amplitudeMS());

		transformer3.FourierTransform(amplitudeMS(), fftAmpMS, false);
		transformer3.setReal(MGVol());
		filter2.FilterShape = RAISED_COSINE;
		filter2.raised_w = 0;//.01;
		filter2.do_generate_3dmask = false;
		filter2.FilterBand = LOWPASS;
		filter2.w1 = (double) idx_freq*step;
		filter2.applyMaskFourierSpace(amplitudeMS(),transformer3.fFourier);
		transformer3.inverseFourierTransform();

		MultidimArray<double> &p2amplitudeMS = MGVol();

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(p2amplitudeMS)
		{
			if (DIRECT_MULTIDIM_ELEM(p2mask, n) > 0 )
			{
				num_energy_mean += DIRECT_MULTIDIM_ELEM(p2amplitudeMS, n);
				count++;
			}
			else
			{
				if (condpremask)
				{
//					r_kij = (k-zdim)*(k-zdim) + (j-ydim)*(j-ydim) + (i-xdim)*(i-xdim);
//					if (r_kij<r2){
//						num_noise_mean += DIRECT_MULTIDIM_ELEM(p2amplitudeMS,n);}
					exit(0);
				}
				else
				{
				num_noise_mean += DIRECT_MULTIDIM_ELEM(p2amplitudeMS,n);
				count2++;
				}
			}
		}
//		std::cout << "count = " << count << std::endl;
//		std::cout << "count2 = " << count2 << std::endl;

		double threshold_energy = num_energy_mean/count;
		double threshold_noise = num_noise_mean/count2;
		std::cout << "energy_mean = " << threshold_energy << std::endl;

		if (idx_freq == 1){
			max_energy_mean =  threshold_energy;
			std::cout << "max_energy_mean = " << max_energy_mean << std::endl;}

		//std::sort(energy_vector.begin(), energy_vector.end());
		//std::sort(noise_vector.begin(), noise_vector.end());

		//len_energy_vector = energy_vector.size();
		//len_noise_vector = noise_vector.size();

		//idx_threshold_energy = floor((double) 0.95*len_energy_vector);
		//idx_threshold_noise = floor((double) 0.99*len_noise_vector);

		//threshold_energy = energy_vector[idx_threshold_energy];
		//threshold_noise = noise_vector[idx_threshold_noise];

		std::cout << "threshold_noise = " << threshold_noise << std::endl;
		std::cout << "threshold_energy = " << threshold_energy << std::endl;


		if ((threshold_noise>threshold_energy) || (threshold_energy<max_energy_mean*0.1)){
			std::cout << "Maximum resolution = " << inv_idx_freq << std::endl;
			std::cout << "-iteration = " << idx_freq << std::endl;
			idx_freq_threshold = idx_freq*step;
			break;
		}
	}

	std::cout << "idx_freq_threshold = " << idx_freq_threshold << std::endl;

	//Obtaining local resolutions
	if (idx_freq_threshold>0)
	{
	outputResolution().initZeros(p2imgVol);

	MultidimArray<double> outputResolution_aux;
	MultidimArray<double> &p2outputResolution_aux = outputResolution_aux, &p2outputResolution = outputResolution();
	MultidimArray<double> &p2mask_amplitudeMS = mask_amplitudeMS;

	outputResolution_aux.initZeros(p2imgVol);


	std::cout << "------------------------------------------------------" << std::endl;
	std::cout << "Calculating local resolution" << std::endl;
	std::cout << "------------------------------------------------------" << std::endl;

	std::cout << "Analysing frequencies..." << std::endl;
	step = (idx_freq_threshold-0.01)/Nsteps;
	std::cout << "step = " << step << std::endl;
	std::cout << "step*Nsteps = " << step*Nsteps << std::endl;

	for (size_t idx_freq = 0; idx_freq<Nsteps; idx_freq++)
	{
		std::cout << "------------------------------------------------------" << std::endl;
		double inv_idx_freq = (double) 1/((idx_freq*step + 0.01));
		std::cout << "freq = " << idx_freq << std::endl;
		std::cout << "inversa freq = " << inv_idx_freq << std::endl;
		std::vector<double> noise_vector(1), energy_vector(1), noise_vector_sorted, energy_vector_sorted;
		num_energy_mean = 0;
		num_noise_mean = 0;
		count = 0;
		count2 = 0;

		if (fType == 0)
		{
			exit(0);
		}
		else
		{

		filter.FilterBand = HIGHPASS;
		filter.w1 = idx_freq*step;
		transformer2.fFourier = transformer.fFourier;
		filter.applyMaskFourierSpace(p2imgVol,transformer2.fFourier);
		transformer2.inverseFourierTransform();
		}

		RieszTransform3Dreal(filteredVol(), RieszVector);
		amplitudeMonogenicSignal3D(filteredVol(), RieszVector, amplitudeMS());

		transformer3.FourierTransform(amplitudeMS(), fftAmpMS, false);
		transformer3.setReal(MGVol());
		filter2.FilterShape = RAISED_COSINE;
		filter2.raised_w = 0.01;
		filter2.do_generate_3dmask = false;
		filter2.FilterBand = LOWPASS;
		filter2.w1 = (double) idx_freq*step;
		filter2.applyMaskFourierSpace(amplitudeMS(),transformer3.fFourier);
		transformer3.inverseFourierTransform();


		MultidimArray<double> &p2amplitudeMS = MGVol();
		applymask(p2amplitudeMS, p2mask, p2mask_amplitudeMS);

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(p2amplitudeMS)
		{
			if (DIRECT_MULTIDIM_ELEM(p2mask, n) == 0 )
			{
				if (condpremask == true)
				{
//					r_kij = (k-zdim)*(k-zdim) + (j-ydim)*(j-ydim) + (i-xdim)*(i-xdim);
//					if (r_kij<r2){
//						noise_vector.push_back(DIRECT_MULTIDIM_ELEM(p2amplitudeMS, n));
//					}
					exit(0);
				}
				else{
					noise_vector.push_back(DIRECT_MULTIDIM_ELEM(p2amplitudeMS, n));
				}
				num_noise_mean += DIRECT_MULTIDIM_ELEM(p2amplitudeMS, n);
				count2++;
			}
		}
		double noise_mean = num_noise_mean/count2;

		//std::sort(noise_vector.begin(), noise_vector.end());
		//len_noise_vector = noise_vector.size();
		//idx_threshold_noise = floor((double) 0.99*len_noise_vector);
		//threshold_noise = noise_vector[idx_threshold_noise];
		threshold_noise = noise_mean;
		std::cout << "threshold_noise = " << threshold_noise << std::endl;



		FOR_ALL_ELEMENTS_IN_ARRAY3D(p2outputResolution_aux)
		{
			if (A3D_ELEM(p2mask_amplitudeMS, k, i, j) > threshold_noise)
				A3D_ELEM(p2outputResolution_aux, k, i, j) = smpr*inv_idx_freq;
		}

		//Image<double> aux_aux = p2outputResolution_aux;
		//FileName fn_aux_MG;
		//fn_aux_MG = fnDir+formatString("/_MG_partial_mask%i.vol", (int) idx_freq);
		//aux_aux.write(fn_aux_MG);
	}

	Image<double> outputResolution2;

	outputResolution() = outputResolution_aux;

	outputResolution.write(fnDir+"/MGresolution.vol");

	std::cout <<  "Finished!" << std::endl;
	}
#endif

#ifdef TESTING2
	std::cout << "Starting ... " << std::endl;
	std::vector< MultidimArray<double > > RieszVector(3), RieszVector_original(3);
	MultidimArray<double> mask_amplitudeMS, mask_amplitudeMS_original;
	Image<double> imgVol, amplitudeMS, mask, outputResolution, amplitudeMS_original, filteredVol;
	double idx_threshold_energy, idx_freq_threshold = 0, idx_threshold_noise, threshold_energy, threshold_noise, max_val, r2, r_kij;
	size_t xdim, ydim, zdim, len_energy_vector, len_noise_vector, count_freq, count_freq2;
	FileName cmd;

	MultidimArray<double> &p2imgVol = imgVol();

	//Parameters:
	double Nsteps = 50;
	double step = 0.01;
	size_t N_elems = (size_t) (0.5/step);

	//Reading input Volume
	imgVol.read(fnVol);
	std::cout << "Volume read!" << std::endl;

	//Obtaining dimensions
	getHalfDimensions(p2imgVol, xdim, ydim, zdim);

	if (condpremask == true)
		r2 = calculateMaskRadius(p2imgVol, xdim, ydim, zdim);

	//Reading mask
	mask.read(fnMask);
	MultidimArray<double> &p2mask = mask();


	MultidimArray<double> &p2amplitudeMS = amplitudeMS();
	MultidimArray<double> mask_filteredvol;
	MultidimArray<double> &p2filteredVol = filteredVol(), &p2mask_filteredvol = mask_filteredvol;

	std::cout << "Analysing frequencies..." << std::endl;
	double num_energy_mean, num_noise_mean;
	size_t count, count2;
	double max_energy_mean=0;
	for (size_t idx_freq = 1; idx_freq<N_elems; idx_freq++)
	{
		std::cout << "------------------------------------------------------" << std::endl;
		num_energy_mean = 0;
		num_noise_mean = 0;
		count = 0;
		count2 = 0;
		double inv_idx_freq = (double) 1/(idx_freq*step);
		std::cout << "freq = " << idx_freq << std::endl;
		std::cout << "Resolution = " << inv_idx_freq << std::endl;

		std::vector<double> noise_vector(1), energy_vector(1), noise_vector_sorted, energy_vector_sorted;

		if (fType == 0){
			cmd = formatString("xmipp_transform_filter -i %s --fourier band_pass %f %f -o %s/filtered_vol.vol", fnVol.c_str(),
							((double) idx_freq*step)-0.002 , ((double) idx_freq*step)+0.002, fnDir.c_str());}
		else{
			cmd = formatString("xmipp_transform_filter -i %s --fourier high_pass %f 0 -o %s/filtered_vol.vol", fnVol.c_str(),
							((double) idx_freq*step) , fnDir.c_str());}

		std::cout << cmd << std::endl;
		system(cmd.c_str());

		std::cout << "Filtrado" << std::endl;
		filteredVol.read(fnDir+"/filtered_vol.vol");

		RieszTransform3Dreal(p2filteredVol, RieszVector);

		amplitudeMonogenicSignal3D(p2filteredVol, RieszVector, p2amplitudeMS);



		//		std::cout << "min = " << p2amplitudeMS.computeMin() << std::endl;


		#ifdef DEBUGGING
		Image<double> debug_Amp_MG = p2amplitudeMS, debug_volfil = filteredVol;
		FileName fndebug_Amp_Mg, fndebug_fil_vol, fndebug_Amp_Mg_filtered;
		fndebug_Amp_Mg = fnDir+formatString("/Amp_MG%i.vol", (int) idx_freq);
		fndebug_Amp_Mg_filtered =  fnDir+formatString("/filtered_MG%i.vol", (int) idx_freq);
		fndebug_fil_vol = fnDir+formatString("/filvol%i.vol", (int) idx_freq);
		debug_Amp_MG.write(fndebug_Amp_Mg);
		debug_volfil.write(fndebug_fil_vol);

		cmd = formatString("xmipp_transform_filter -i %s --fourier low_pass %f 0 -o %s", fndebug_Amp_Mg.c_str(),
											((double) idx_freq*step), fndebug_Amp_Mg_filtered.c_str(), (int) idx_freq);

		std::cout << cmd << std::endl;
		system(cmd.c_str());


		Image<double> MG_filtered;

		MG_filtered.read(fndebug_Amp_Mg_filtered);

		MultidimArray<double> &p2amplitudeMS = MG_filtered();
#endif
		p2amplitudeMS.printShape();
		p2mask.printShape();

		FOR_ALL_ELEMENTS_IN_ARRAY3D(p2amplitudeMS)
		{
			if (A3D_ELEM(p2mask, k, i, j) > 0 )
			{
				energy_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
				num_energy_mean = A3D_ELEM(p2amplitudeMS, k, i, j) + num_energy_mean;
				count++;
			}
			else
			{
				if (condpremask == true)
				{
					r_kij = (k-zdim)*(k-zdim) + (j-ydim)*(j-ydim) + (i-xdim)*(i-xdim);
					if (r_kij<r2)
						noise_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
				}
				else
					noise_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));

				num_noise_mean = A3D_ELEM(p2amplitudeMS, k, i, j) + num_noise_mean;
				count2++;
			}
		}


		std::cout << "count = " << count << std::endl;
		std::cout << "count2 = " << count2 << std::endl;

		double energy_mean = num_energy_mean/count;
		double noise_mean = num_noise_mean/count2;
		std::cout << "energy_mean = " << energy_mean << std::endl;

		if (idx_freq == 1)
			max_energy_mean =  energy_mean;

		//std::sort(energy_vector.begin(), energy_vector.end());
		//std::sort(noise_vector.begin(), noise_vector.end());

		//len_energy_vector = energy_vector.size();
		//len_noise_vector = noise_vector.size();

		//idx_threshold_energy = floor((double) 0.95*len_energy_vector);
		//idx_threshold_noise = floor((double) 0.99*len_noise_vector);

		//threshold_energy = energy_vector[idx_threshold_energy];
		//threshold_noise = noise_vector[idx_threshold_noise];
		threshold_energy = energy_mean;
		threshold_noise = noise_mean;

		std::cout << "threshold_noise = " << threshold_noise << std::endl;
		std::cout << "threshold_energy = " << threshold_energy << std::endl;


		if ((threshold_noise>threshold_energy) || (threshold_energy<max_energy_mean*0.1)){
			std::cout << "Maximum resolution = " << inv_idx_freq << std::endl;
			std::cout << "-iteration = " << idx_freq << std::endl;
			idx_freq_threshold = idx_freq*step;
			break;
		}
	}

	std::cout << "idx_freq_threshold = " << idx_freq_threshold << std::endl;

	//Obtaining local resolutions
	if (idx_freq_threshold>0)
	{
	outputResolution().initZeros(p2imgVol);

	MultidimArray<double> outputResolution_aux;
	MultidimArray<double> &p2outputResolution_aux = outputResolution_aux, &p2outputResolution = outputResolution();
	MultidimArray<double> &p2mask_amplitudeMS = mask_amplitudeMS;

	outputResolution_aux.initZeros(p2imgVol);

	std::cout << "------------------------------------------------------" << std::endl;
	std::cout << "Calculating local resolution" << std::endl;
	std::cout << "------------------------------------------------------" << std::endl;

	std::cout << "Analysing frequencies..." << std::endl;
	step = (idx_freq_threshold-0.01)/Nsteps;
	std::cout << "step = " << step << std::endl;
	std::cout << "step*Nsteps = " << step*Nsteps << std::endl;

	for (size_t idx_freq = 0; idx_freq<Nsteps; idx_freq++)
	{
		std::cout << "------------------------------------------------------" << std::endl;
		double inv_idx_freq = (double) 1/((idx_freq*step + 0.01));
		std::cout << "freq = " << idx_freq << std::endl;
		std::cout << "inversa freq = " << inv_idx_freq << std::endl;
		std::vector<double> noise_vector(1), energy_vector(1), noise_vector_sorted, energy_vector_sorted;
		num_energy_mean = 0;
		num_noise_mean = 0;
		count = 0;
		count2 = 0;

		if (fType == 0)
		{
		cmd = formatString("xmipp_transform_filter -i %s --fourier band_pass %f %f -o %s/filtered_vol.vol", fnVol.c_str(),
							((double) (idx_freq*step +0.01))-0.002 , ((double) (idx_freq*step +0.01))+0.002, fnDir.c_str());
		}
		else
		{
		cmd = formatString("xmipp_transform_filter -i %s --fourier high_pass %f 0 -o %s/filtered_vol.vol", fnVol.c_str(),
							((double) (idx_freq*step +0.01)) , fnDir.c_str());
		}
		std::cout << cmd << std::endl;
		system(cmd.c_str());

		filteredVol.read(fnDir+"/filtered_vol.vol");

		RieszTransform3Dreal(p2filteredVol, RieszVector);

		amplitudeMonogenicSignal3D(p2filteredVol, RieszVector, p2amplitudeMS);

		#ifdef DEBUGGING
		Image<double> debug_Amp_MG = p2amplitudeMS, debug_volfil = p2filteredVol;
		FileName fndebug_Amp_Mg, fndebug_fil_vol, fndebug_Amp_Mg_filtered;;
		fndebug_Amp_Mg = fnDir+formatString("/_Amp_MG%i.vol", (int) idx_freq);
		fndebug_Amp_Mg_filtered =  fnDir+formatString("/_filtered_MG%i.vol", (int) idx_freq);
		fndebug_fil_vol = fnDir+formatString("/_filvol%i.vol", (int) idx_freq);
		debug_Amp_MG.write(fndebug_Amp_Mg);
		debug_volfil.write(fndebug_fil_vol);



		cmd = formatString("xmipp_transform_filter -i %s --fourier low_pass %f 0 -o %s", fndebug_Amp_Mg.c_str(),
											((double) idx_freq*step) , fndebug_Amp_Mg_filtered.c_str());

		std::cout << cmd << std::endl;
		system(cmd.c_str());

		Image<double> MG_filtered;

		MG_filtered.read(fndebug_Amp_Mg_filtered);
		MultidimArray<double> &p2amplitudeMS = MG_filtered();

		#endif
		applymask(p2amplitudeMS, p2mask, p2mask_amplitudeMS);

		FOR_ALL_ELEMENTS_IN_ARRAY3D(p2amplitudeMS)
		{
			if (A3D_ELEM(p2mask, k, i, j) == 0 )
			{
				if (condpremask == true)
				{
					r_kij = (k-zdim)*(k-zdim) + (j-ydim)*(j-ydim) + (i-xdim)*(i-xdim);
					if (r_kij<r2){
						noise_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
					}
				}
				else{
					noise_vector.push_back(A3D_ELEM(p2amplitudeMS, k, i, j));
				}
				num_noise_mean = A3D_ELEM(p2amplitudeMS, k, i, j) + num_noise_mean;
				count2++;
			}
		}
		double noise_mean = num_noise_mean/count2;

		//std::sort(noise_vector.begin(), noise_vector.end());

		//len_noise_vector = noise_vector.size();

		//idx_threshold_noise = floor((double) 0.99*len_noise_vector);

		//threshold_noise = noise_vector[idx_threshold_noise];
		threshold_noise = noise_mean;
		std::cout << "threshold_noise = " << threshold_noise << std::endl;



		FOR_ALL_ELEMENTS_IN_ARRAY3D(p2outputResolution_aux)
		{
			if (A3D_ELEM(p2mask_amplitudeMS, k, i, j) > threshold_noise)
				A3D_ELEM(p2outputResolution_aux, k, i, j) = smpr*inv_idx_freq;
		}

//		if (count_freq2>0)
//		{
//			std::cout << "opening... " << std::endl;
			Image<double> aux_aux = p2outputResolution_aux;
			FileName fn_aux_MG;
			fn_aux_MG = fnDir+formatString("/_MG_partial_mask%i.vol", (int) idx_freq);
			aux_aux.write(fn_aux_MG);
//			opening3D(p2outputResolution_aux, p2outputResolution, 6, 0, 1);
//		}

	}

	Image<double> outputResolution2;

	outputResolution() = outputResolution_aux;
	//medianFilter3x3x3(outputResolution(), outputResolution2());


	outputResolution.write(fnDir+"/MGresolution.vol");

	std::cout <<  "Finished!" << std::endl;
	}
#endif

}
