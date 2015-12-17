/***************************************************************************
 * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
 *
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

#include "angular_tilt_pair_assignment.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <data/xmipp_image.h>
#include <data/filters.h>
#include <data/projection.h>
#include <data/transformations.h>
#include <geometry.h>
#include <data/matrix1d.h>
#include <data/matrix2d.h>
#include <data/xmipp_image_base.h>
#include <data/projection.h>

//#define DEBUG
//#define DEBUG_CHECK_REF_VOLUME
//#define DEBUG_SHIFT


void ProgTiltPairassignment::readParams()
{
	fnUntilt = getParam("--untiltparticles");
	fnTilt = getParam("--tiltparticles");
	alphaT = getDoubleParam("--alphaT");
	tilt_mic = getDoubleParam("--tiltmic");
	alphaU = getDoubleParam("--alphaU");
	fnVol = getParam("--vol");
	fnDir = getParam("--odir");
	smprt = getDoubleParam("--angular_sampling");
	maxshift = getDoubleParam("--maxshift");
	fnRefVol = getParam("--ref");
}


void ProgTiltPairassignment::defineParams()
{
	//usage
	addUsageLine("From an untilted particles stack and an initial volume, an angular assignment of the untilted particles is performed."
			"Then, knowing the tilt axis and the tilt angle of the micrography, tiltmic, and the rotational angles of the untilted and"
			"tilted micrographs, alphaU and alphaT respectively, an angular assignment of the tilted particles is determined");
	//params
	addParamsLine("  [--untiltparticles <md_file=\"\">]    : Untilt particles stack");
	addParamsLine("  [--tiltparticles <md_file=\"\">]    : Tilt particles stack");
	addParamsLine("  [--vol <img_file=\"\">]    : Input reference volume");
	addParamsLine("  [--ref <img_file=\"\">]    : Reference volume");
	addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
	addParamsLine("  [--angular_sampling <s=5>]   : Angular sampling rate in degrees. This sampling represents the accuracy for assigning directions");
	addParamsLine("  [--alphaU <s=0>]   : Rotational angle of the untilted micrograph (degrees)");
	addParamsLine("  [--alphaT <s=0>]   : Rotational angle of the tilted micrograph (degrees)");
	addParamsLine("  [--tiltmic <s=0>]   : Tilt angle of the micrograph (degrees)");
	addParamsLine("  [--maxshift <s=10>]   : Maximum shift for aligning images (in pixels)");
}


void ProgTiltPairassignment::generateProjections(FileName fnVol, double smprt)
{
	FileName fnGallery, fnGalleryMetaData;

	// Generate projections
	fnGallery=formatString("%s/gallery.stk",fnDir.c_str());
	//fnAngles=formatString("%s/angles.xmd",fnDir.c_str());
	String args=formatString("-i %s -o %s --sampling_rate %f",
			fnVol.c_str(),fnGallery.c_str(),smprt);
			//We are considering the psi sampling = angular sampling rate

	std::cout << args << std::endl;
	String cmd=(String)"xmipp_angular_project_library " + args;
	system(cmd.c_str());
}


void ProgTiltPairassignment::run()
{
	std::cout << "Starting..." << std::endl;

	// Generating projection from the volume, fnVol, with an angular sampling rate of smprt degrees
	std::cout << "Generating projections" << std::endl;
	generateProjections(fnVol, smprt);
		//State: Finished

	//Reading particles from the untilted stack and projections
	MetaData mduntilt_exp, mdtilt_exp, mdproj;
	FileName fnprojection = fnDir+"/gallery.doc";
	mduntilt_exp.read(fnUntilt);
	mdtilt_exp.read(fnTilt);
	mdproj.read(fnprojection);
	// allGalleryProjection;

	FileName fnprojection_stk = fnDir+"/gallery.stk";
	Image<double> imgstack;
	imgstack.read(fnprojection_stk);
	const MultidimArray<double> allGalleryProjection = imgstack();

	#ifdef DEBUG
	size_t size_mduntilt_exp = mduntilt_exp.size();
	size_t size_mdtilt_exp = mdtilt_exp.size();
	size_t size_mdproj = mdproj.size();
	size_t Mdim_Xdim, Mdim_Ydim, Mdim_Zdim, Mdim_Ndim;
	allGalleryProjection.getDimensions(Mdim_Xdim, Mdim_Ydim, Mdim_Zdim, Mdim_Ndim);
	std::cout << "  " << std::endl;
	std::cout << "---------------------------------" << std::endl;
	std::cout << "READING PARTICLES FROM STACKS AND PROJECTIONS" << std::endl;
	std::cout << "MetaData Size (Untilted     mduntilt_exp)  : " << size_mduntilt_exp <<  std::endl;
	std::cout << "MetaData Size (Tilted       mdtilt_exp)    : " << size_mdtilt_exp <<  std::endl;
	std::cout << "MetaData Size (Projections  mdproj)        : " << size_mdproj <<  std::endl;

	std::cout << "MultidimArray Size             (allGalleryProjection)  : " << Mdim_Xdim << " x " <<  Mdim_Ydim << " x " << Mdim_Zdim << std::endl;
	std::cout << "MultidimArray Number of images (allGalleryProjection)  : " << Mdim_Ndim << std::endl;

	if (size_mduntilt_exp != size_mdtilt_exp)
	{
		std::cerr << "WARNING: Untilted (mduntilt_exp) and Tilted (mdtilt_exp) MetaData dont have the same number of particles."
				"They may not be tilt pairs of particles" << std::endl;
	}
	std::cout << "---------------------------------" << std::endl;
	std::cout << "  " << std::endl;
	#endif
		//State: Finished


	//Reading untilted, tilted particles and the projections of the volume,then they are stored them into string vectors.
	//Thus we are avoiding accessing once and again to metadata
	std::vector<std::string> Untilted_filenames, Tilted_filenames, proj_filenames;
	FileName fnuntilt_exp, fntilt_exp, fnproj;
	double rot, tilt;

	Matrix2D<double> angles_rot_tilt;
	angles_rot_tilt.initZeros(mdproj.size(),2); //first column is the rot angle and second column is the tilt angle

	FOR_ALL_OBJECTS_IN_METADATA(mduntilt_exp)
	{
		mduntilt_exp.getValue(MDL_IMAGE, fnuntilt_exp, __iter.objId);
		Untilted_filenames.push_back(fnuntilt_exp);
	}
	std::cout << "Untilt MetaData read" << std::endl;

	FOR_ALL_OBJECTS_IN_METADATA(mdtilt_exp)
	{
		mdtilt_exp.getValue(MDL_IMAGE, fntilt_exp, __iter.objId);
		Tilted_filenames.push_back(fntilt_exp);
	}
	std::cout << "Tilt MetaData read" << std::endl;
	FOR_ALL_OBJECTS_IN_METADATA(mdproj)
	{
		mdproj.getValue(MDL_IMAGE, fnproj, __iter.objId);
		mdproj.getValue(MDL_ANGLE_ROT, rot, __iter.objId);
		mdproj.getValue(MDL_ANGLE_TILT, tilt, __iter.objId);
		proj_filenames.push_back(fnproj);

		MAT_ELEM(angles_rot_tilt, (__iter.objId)-1, 0) = rot;
		MAT_ELEM(angles_rot_tilt, (__iter.objId)-1, 1) = tilt;
	}
	std::cout << "Projection MetaData gallery read" << std::endl;
	size_t len_u = Untilted_filenames.size();
	size_t len_t = Tilted_filenames.size();
	size_t len_p = proj_filenames.size();
//	std::cout << "Projecction filenames = " << proj_filenames[0] << std::endl;
//	std::cout << "Projecction filenames = " << proj_filenames[len_p-1] << std::endl;
//
//	for (size_t j =0; j<(len_p); j++)
//	{
//		MultidimArray<double> imgGallery;
//		imgGallery.aliasImageInStack(allGalleryProjection,j);
//		Image<double> imag_projected;
//		imag_projected = imgGallery;
//		imag_projected.write("imag_projected.xmp");
//		std::cout << "Imagen = " << proj_filenames[j] << std::endl;
//		std::cout << "Press any key" << std::endl;
//		int ccc;
//		ccc= getchar();
//	}


	#ifdef DEBUG
	size_mduntilt_exp = mduntilt_exp.size();
	size_mdtilt_exp = mdtilt_exp.size();
	size_mdproj = mdproj.size();

	std::cout << "  " << std::endl;
	std::cout << "---------------------------------" << std::endl;
	std::cout << "CONVERSION FROM METADATA TO MATRIX2D" << std::endl;
	std::cout << "MetaData Size (Untilted     mduntilt_exp)  : " << size_mduntilt_exp <<  std::endl;
	std::cout << "MetaData Size (Tilted       mdtilt_exp)    : " << size_mdtilt_exp <<  std::endl;
	std::cout << "MetaData Size (Projections  mdproj)        : " << size_mdproj <<  std::endl;

	std::cout << "Matrix2D Size (Untilted filenames     Untilted_filenames)  : " << len_u <<  std::endl;
	std::cout << "Matrix2D Size (Tilted filenames       Tilted_filenames)    : " << len_t <<  std::endl;
	std::cout << "Matrix2D Size (Projections filenames  proj_filenames)      : " << len_p <<  std::endl;

	std::cout << "Matrix2D Size (Angles - Rot Tilt     angles_rot_tilt)      : " << MAT_YSIZE(angles_rot_tilt) <<
			" x " << MAT_XSIZE(angles_rot_tilt) <<  std::endl;


	if (size_mduntilt_exp != size_mdtilt_exp)
	{
		std::cerr << "WARNING: Untilted (mduntilt_exp) and Tilted (mdtilt_exp) MetaData dont have the same number of particles."
				"They may not be tilt pairs of particles" << std::endl;
	}
	if ( (size_mduntilt_exp != len_u) || (size_mdtilt_exp != len_t) || (size_mdproj != len_p) )
		{
			std::cerr << "WARNING: Size of MetaData Untilted/Tilted/projection is not the same that their correspondent Untilted/Tilted/projection Matrix2D" << std::endl;
		}
	if ( MAT_YSIZE(angles_rot_tilt) != len_p )
			{
				std::cerr << "WARNING: Size of MetaData or Matrix2D does not match to the number of rows of Matrix2D angles_rot_tilt." << std::endl;
			}
	std::cout << "---------------------------------" << std::endl;
	std::cout << "  " << std::endl;
	#endif


	//Cleaning memory
	mduntilt_exp.clear();
	mdtilt_exp.clear();
	mdproj.clear();
	imgstack.clear();
		//State: Finished


	//For each experimental untilted image, an angular assignment is performed
	std::cout << "Searching correlations" << std::endl;
	Image<double> ImgUn_exp, ImgT_exp, Img_proj, tilt_proj;
	MultidimArray<double> Multidim_exp, Multidim_proj, Multidim_vol;
	Matrix2D<double> transformation_matrix;
	Matrix2D<double> position_u_gallery_and_psi, ZYZ_u, ZYZ_angles, ZYZ_t;

	position_u_gallery_and_psi.initZeros(14,len_u);  //Row1 k position of the untilted vector.
													//Row2 j position of the gallery vector proj_filenames.
													//Row3 rot_u angle from the alignment.
													//Row4 tilt_u angle from the alignment.
													//Row5 psi_u angle from the alignment.
													//Row6 rot_t angle from the alignment.
													//Row7 tilt_t angle from the alignment.
													//Row8 psi_t angle from the alignment.

	double rot_t, tilt_t, psi_t, shiftX, shiftY;

	Image<double> Img_vol;
	Img_vol.read(fnVol);
	Img_vol().setXmippOrigin();

	ImgUn_exp.read(Untilted_filenames[0]);
	size_t Ydim, Xdim, Zdim, Ndim;
	ImgUn_exp.getDimensions(Xdim, Ydim, Zdim, Ndim);

	CorrelationAux Auxcorr_bestshift;

#ifdef DEBUG
	size_t XTdim, YTdim, ZTdim, NTdim, XVoldim, YVoldim, ZVoldim, NVoldim;
	ImgT_exp.read(Tilted_filenames[0]);
	ImgT_exp.getDimensions(XTdim, YTdim, ZTdim, NTdim);
	Img_vol.getDimensions(XVoldim, YVoldim, ZVoldim, NVoldim);
	//allGalleryProjection.getDimensions(Mdim_Xdim, Mdim_Ydim, Mdim_Zdim, Mdim_Ndim);
	std::cout << "  " << std::endl;
	std::cout << "---------------------------------" << std::endl;
	std::cout << "PREPARING DATA FOR ITERATIONS AND CORRELATIONS" << std::endl;
	std::cout << "Untilted Particle Size : " << Xdim << " x " << Ydim << std::endl;
	std::cout << "Tilted Particle Size   : " << XTdim << " x " << YTdim << std::endl;
	std::cout << "Volume Size            : " << XVoldim << " x " << YVoldim  << " x " << ZVoldim << std::endl;

	std::cout << "Matrix2D Size (Untilted filenames     Untilted_filenames)  : " << len_u <<  std::endl;
	std::cout << "Matrix2D Size (Tilted filenames       Tilted_filenames)    : " << len_t <<  std::endl;

	std::cout << "Matrix2D Size (Outputdata     position_u_gallery_and_psi)  : " << MAT_YSIZE(position_u_gallery_and_psi) <<
			" x " << MAT_XSIZE(position_u_gallery_and_psi) <<  std::endl;

	if ( (Xdim != XTdim) || (Ydim != YTdim) )
	{
		std::cerr << "WARNING: Untilted and Tilted particles not have the same size." << std::endl;
	}
	if ( (Xdim != XVoldim) || (Ydim != YVoldim) || (Xdim != ZVoldim) || (Ydim != ZVoldim))
		{
			std::cerr << "WARNING: Particles and volume may not have the same size." << std::endl;
		}
	if ( (len_u != MAT_XSIZE(position_u_gallery_and_psi)))
			{
				std::cerr << "WARNING: position_u_gallery_and_psi matrix does not have the same number columns "
						"than the number of particles" << std::endl;
			}

	std::cout << "---------------------------------" << std::endl;
	std::cout << "  " << std::endl;
	#endif

#ifdef DEBUG
	//Preparing counters for checking for loop
	size_t iteration_number_k = 0;
	size_t iteration_number_j;

#endif

	Matrix1D<double> Shiftvec(2);
	MultidimArray<double> imgGallery;

	for (size_t k=0; k<(len_u); k++)
	{
		#ifdef DEBUG
			iteration_number_k = iteration_number_k + 1;
			iteration_number_j = 0;
			bool flag = false;
		#endif
		ImgUn_exp.read(Untilted_filenames[k]);	//Reading image
		ImgT_exp.read(Tilted_filenames[k]);
		ImgUn_exp().setXmippOrigin();
		ImgT_exp().setXmippOrigin();
		MAT_ELEM(position_u_gallery_and_psi,0,k) = k;
		std::cout << "Iteration:  " << k << "/" << len_u-1 << std::endl;
		Projection tilt_proj;
		double corr1 = 0, corr2=0, corr=0, bestcorr1 = 0, bestcorr = 0;
		MultidimArray<double> Improj_traslated;

		for (size_t j =0; j<(len_p); j++)
		{
			//Untilt assignment
			imgGallery.aliasImageInStack(allGalleryProjection,j);
			MultidimArray<double> imgGallery_orig = imgGallery;
			imgGallery_orig.setXmippOrigin();
			#ifdef DEBUG
			Image<double> prueba;
			prueba = imgGallery;
			prueba.write("antes_alignImage.xmp");
			std::cout << "Projection   Image = " << proj_filenames[j] << std::endl;
			std::cout << "j = " << j << std::endl;


				iteration_number_j = iteration_number_j + 1;
			#endif

			//////////////////////////////////////////
			//CORRELATION UNTILT AND PROJECTIONS
			//corr1 = alignImages(ImgUn_exp(), imgGallery, transformation_matrix, true);
			corr1 = alignImages(ImgUn_exp(), imgGallery_orig, transformation_matrix, true);


			if ((fabs(MAT_ELEM(transformation_matrix, 0, 2)) > maxshift) || (fabs(MAT_ELEM(transformation_matrix, 1, 2)) > maxshift))
				continue;

			if ((corr1 <0.7*bestcorr1) || (corr1<0))
				continue;


			//////////////////////////////////////////
			//UNTILT ASSIGNMENT
			double psi = acos( 0.5*( MAT_ELEM(transformation_matrix,0,0) + MAT_ELEM(transformation_matrix,1,1) ) )*180/PI;
			//rot  = MAT_ELEM(angles_rot_tilt, j, 0);
			//tilt = MAT_ELEM(angles_rot_tilt, j, 1);
			//Tilted angular assignment
//#ifdef ONLYUNTILT
			//////////////////////////////////////////
			//TILT ASSIGNMENT
			Euler_angles2matrix(MAT_ELEM(angles_rot_tilt, j, 0), MAT_ELEM(angles_rot_tilt, j, 1), psi, ZYZ_u);
			Euler_angles2matrix(alphaT, tilt_mic, -alphaU, ZYZ_angles);
			ZYZ_t = ZYZ_angles*ZYZ_u;
			Euler_matrix2angles(ZYZ_t, rot_t, tilt_t, psi_t);


			projectVolume(Img_vol(), tilt_proj, (int) Ydim, (int) Xdim, rot_t, tilt_t, psi_t);
			tilt_proj().statisticsAdjust(0,1);
			tilt_proj().setXmippOrigin();

			//////////////////////////////////////////
			//CORRELATION TILT AND PROJECTIONS
			corr2 = bestShift(ImgT_exp(), tilt_proj(),  shiftX, shiftY, Auxcorr_bestshift, NULL, maxshift);
			VEC_ELEM(Shiftvec,0) = shiftX;
			VEC_ELEM(Shiftvec,1) = shiftY;

			Improj_traslated = tilt_proj();
			selfTranslate(1, Improj_traslated, Shiftvec);

			corr2 = correlationIndex(ImgT_exp(), Improj_traslated, NULL, NULL);
			if (corr2<0)
				continue;

//#endif
			//std::cout << "CORR2 = " << corr2 << std::endl;
			corr = 0.5*(corr1 + corr2);

			//////////////////////////////////////////
			//OUTPUT ANGLES AND ASSIGNMENTS
			if (corr>bestcorr)
			{

				bestcorr = corr;
				MAT_ELEM(position_u_gallery_and_psi, 1, k) = j;
				MAT_ELEM(position_u_gallery_and_psi, 2, k) = MAT_ELEM(angles_rot_tilt, j, 0);
				MAT_ELEM(position_u_gallery_and_psi, 3, k) = MAT_ELEM(angles_rot_tilt, j, 1);
				MAT_ELEM(position_u_gallery_and_psi, 4, k) = psi;
				MAT_ELEM(position_u_gallery_and_psi, 5, k) = rot_t;
				MAT_ELEM(position_u_gallery_and_psi, 6, k) = tilt_t;
				MAT_ELEM(position_u_gallery_and_psi, 7, k) = psi_t;
				MAT_ELEM(position_u_gallery_and_psi, 8, k) = corr1;
				MAT_ELEM(position_u_gallery_and_psi, 9, k) = corr2;
				MAT_ELEM(position_u_gallery_and_psi, 10, k) = -MAT_ELEM(transformation_matrix,0,2);  // X Shift  in untilted image
				MAT_ELEM(position_u_gallery_and_psi, 11, k) = -MAT_ELEM(transformation_matrix,1,2);  // Y Shift
				MAT_ELEM(position_u_gallery_and_psi, 12, k) = -VEC_ELEM(Shiftvec,0);  // X Shift in tilted image
				MAT_ELEM(position_u_gallery_and_psi, 13, k) = -VEC_ELEM(Shiftvec,1);  // Y Shift
				#ifdef DEBUG
					Image<double> save;
					save()=ImgUn_exp();
					save.write("PPPuntilted.xmp");

					save()=ImgT_exp();
					save.write("PPPtilted.xmp");

					save()=imgGallery;
					save.write("PPPprojection_untilted.xmp");

					save()=imgGallery_orig;
					save.write("PPPprojection_untilted_orig.xmp");

					save()=Improj_traslated;
					save.write("PPPprojection_tilted.xmp");

					std::cout << "Corr1 =" << corr1 << std::endl;
					std::cout << "Corr2 =" << corr2 << std::endl;
					std::cout << "Corr  =" << corr << std::endl;
					std::cout << "                          " << std::endl;

					std::cout << "Untilted ROT = " << MAT_ELEM(angles_rot_tilt, j, 0) << std::endl;
					std::cout << "Untilted TILT = " << MAT_ELEM(angles_rot_tilt, j, 1) << std::endl;
					std::cout << "Untilted PSI = " << psi << std::endl;

					std::cout << "                          " << std::endl;
					std::cout << "Projection   Image = " << proj_filenames[j] << std::endl;
					std::cout << "Experimental Image = " << Untilted_filenames[k] << std::endl;
//					std::cout << "Tilted ROT = " << rot_t << std::endl;
//					std::cout << "Tilted TILT = " << tilt_t << std::endl;
//					std::cout << "Tilted PSI = " << psi_t << std::endl;
					std::cout << "------------------------------------" << std::endl;
					std::cout << "------------------------------------" << std::endl;


//					String cmd=(String)"scipion viewer PPPuntilted.xmp  PPPtilted.xmp PPPprojection_untilted.xmp PPPprojection_tilted.xmp &";
//					system(cmd.c_str());
					std::cout << "Press any key" << std::endl;
					int c;
					c= getchar();
				#endif
			}
		}
	}

#ifdef DEBUG
std::cout << "  " << std::endl;
std::cout << "---------------------------------" << std::endl;
std::cout << "Matrix2D Size (Untilted filenames     Untilted_filenames)  : " << len_u <<  std::endl;
std::cout << "Matrix2D Size (Tilted filenames       Tilted_filenames)    : " << len_t <<  std::endl;
std::cout << "Number of iterations for projections          : " << iteration_number_j <<  std::endl;
std::cout << "Number of iterations for experimental data    : " << iteration_number_k <<  std::endl;


if (iteration_number_j != len_p)
{
	std::cerr << "WARNING: Number of projections and number of iterations do not match" << std::endl;
}
if (iteration_number_k != len_u)
{
	std::cerr << "WARNING: Number of projections and number of iterations do not match" << std::endl;
}

std::cout << "---------------------------------" << std::endl;
std::cout << "  " << std::endl;
#endif

	std::cout << "Correlation ended " << std::endl;


	//Storing assignments into output metadata
	MetaData mduntilt_output, mdtilt_output;
	size_t objId_un, objId_t;
	for (size_t k=0; k<(len_u); k++)
	{
		objId_un = mduntilt_output.addObject();
		mduntilt_output.setValue(MDL_IMAGE, Untilted_filenames[k], objId_un);
		mduntilt_output.setValue(MDL_PARTICLE_ID, k, objId_un);
		mduntilt_output.setValue(MDL_ANGLE_ROT, MAT_ELEM(position_u_gallery_and_psi, 2, k), objId_un);
		mduntilt_output.setValue(MDL_ANGLE_TILT, MAT_ELEM(position_u_gallery_and_psi, 3, k), objId_un);
		mduntilt_output.setValue(MDL_ANGLE_PSI, MAT_ELEM(position_u_gallery_and_psi, 4, k), objId_un);
		mduntilt_output.setValue(MDL_MAXCC, MAT_ELEM(position_u_gallery_and_psi, 8, k), objId_un);
		mduntilt_output.setValue(MDL_SHIFT_X, MAT_ELEM(position_u_gallery_and_psi, 10, k), objId_un);
		mduntilt_output.setValue(MDL_SHIFT_Y, MAT_ELEM(position_u_gallery_and_psi, 11, k), objId_un);

		objId_t = mdtilt_output.addObject();
		mdtilt_output.setValue(MDL_IMAGE, Tilted_filenames[k], objId_t);
		mdtilt_output.setValue(MDL_PARTICLE_ID, k, objId_t);
		mdtilt_output.setValue(MDL_ANGLE_ROT, MAT_ELEM(position_u_gallery_and_psi, 5, k), objId_t);
		mdtilt_output.setValue(MDL_ANGLE_TILT, MAT_ELEM(position_u_gallery_and_psi, 6, k), objId_t);
		mdtilt_output.setValue(MDL_ANGLE_PSI, MAT_ELEM(position_u_gallery_and_psi, 7, k), objId_t);
		mdtilt_output.setValue(MDL_MAXCC, MAT_ELEM(position_u_gallery_and_psi, 9, k), objId_t);
		mdtilt_output.setValue(MDL_SHIFT_X, MAT_ELEM(position_u_gallery_and_psi, 12, k), objId_t);
		mdtilt_output.setValue(MDL_SHIFT_Y, MAT_ELEM(position_u_gallery_and_psi, 13, k), objId_t);
	}

	mduntilt_output.write((String)"particles@"+fnDir+'/'+fnUntilt.getBaseName() +"_angular_assignment"+ ".xmd" );
	mdtilt_output.write((String)"particles@"+fnDir+'/'+fnTilt.getBaseName() +"_angular_assignment"+ ".xmd" );
}
