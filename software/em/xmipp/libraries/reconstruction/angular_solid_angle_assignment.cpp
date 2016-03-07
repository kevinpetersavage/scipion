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

#include "angular_solid_angle_assignment.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <data/xmipp_image.h>
#include <data/multidim_array.h>
#include <data/projection.h>
#include "fourier_projection.h"
#include <data/transformations.h>
#include <geometry.h>
#include <data/matrix1d.h>
#include <data/matrix2d.h>
#include <data/xmipp_image_base.h>
#include <data/projection.h>
#include <time.h>
#include "fourier_filter.h"
//#define DEBUG



void ProgSolidAngleAssignment::readParams()
{
	fnparticles = getParam("--particles");
	alphaT = getDoubleParam("--alphaT");
	tilt_mic = getDoubleParam("--tiltmic");
	alphaU = getDoubleParam("--alphaU");
	fnDir = getParam("--odir");
	fnSym = getParam("--sym");
	smprt = getDoubleParam("--angular_sampling");
	maxshift = getDoubleParam("--maxshift");
	fnVol = getParam("--initvol");
	pad = getIntParam("--pad");
	fnmic = getParam("--micangles");
}


void ProgSolidAngleAssignment::defineParams()
{
	//usage
	addUsageLine("From an untilted particles stack and an initial volume, an angular assignment of the untilted particles is performed."
			"Then, knowing the tilt axis and the tilt angle of the micrography, tiltmic, and the rotational angles of the untilted and"
			"tilted micrographs, alphaU and alphaT respectively, an angular assignment of the tilted particles is determined");
	//params
	addParamsLine("  [--particles <md_file=\"\">]    : Particles stack");
	addParamsLine("  [--alphaT <s=0>]   : Rotational angle of the tilted micrograph (degrees)");
	addParamsLine("  [--tiltmic <s=0>]   : Tilt angle of the micrograph (degrees)");
	addParamsLine("  [--alphaU <s=0>]   : Rotational angle of the untilted micrograph (degrees)");
	addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
	addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
	addParamsLine("  [--angular_sampling <s=5>]   : Angular sampling rate in degrees. This sampling represents the accuracy for assigning directions");
	addParamsLine("  [--maxshift <s=10>]   : Maximum shift for aligning images (in pixels)");
	addParamsLine("  [--initvol <md_file=\"\">]   : Initial Volume");
	addParamsLine("  [--pad <s=10>]   : Padding factor");
	addParamsLine("  [--micangles <md_file=\"\">]   : Metadata with micrograph angles, rot, tilt, psi and the micrograph_Id");
}


void ProgSolidAngleAssignment::assignAngles(const MetaData mduntilt_exp, FileName fnun_out)
{

	size_t len_p, idxp, Xdim, Ydim, Zdim, Ndim, len_u;
	Matrix2D<double> angles_rot_tilt, transformation_matrix, position_u_gallery_and_psi;
	MetaData mdproj;
	Image<double> imgstack, Img_vol, ImgUn_exp;
	MultidimArray<double> ImgUn_exp_copy, imgGallery, imgGallery_orig, ImgUn_exp_copy2;
	double rot, tilt, psi;
	FileName Untilted_filename_aux, fnuntilt_exp;
	std::vector<std::string> Untilted_filenames;
	std::vector< AlignmentTransforms> galleryTransform_;
	AlignmentAux aux_u;
	RotationalCorrelationAux aux3;
	CorrelationAux aux2;

	std::cout << "Generating projections" << std::endl;
	String args=formatString("-i %s -o gallery.stk --sampling_rate %f", fnVol.c_str(), smprt);
			//We are considering the psi sampling = angular sampling rate

	String cmd=(String)"xmipp_angular_project_library " + args;
	system(cmd.c_str());

	FileName fnprojection = formatString("%s/gallery.doc",fnDir.c_str());
	mdproj.read(fnprojection);
	len_p = mdproj.size();


	imgstack.read(fnDir+"/gallery.stk");

	const MultidimArray<double> allGalleryProjection = imgstack();
		//State: Finished

	angles_rot_tilt.initZeros(len_p,2); //first column is the rot angle and second column is the tilt angle

	idxp = 0;
	FOR_ALL_OBJECTS_IN_METADATA(mdproj)
	{
		mdproj.getValue(MDL_ANGLE_ROT, rot, __iter.objId);
		mdproj.getValue(MDL_ANGLE_TILT, tilt, __iter.objId);

		MAT_ELEM(angles_rot_tilt, idxp, 0) = rot;
		MAT_ELEM(angles_rot_tilt, idxp, 1) = tilt;
		++idxp;
	}

	//For each experimental untilted image, an angular assignment is performed
	std::cout << "Searching correlations" << std::endl;

	Img_vol.read(fnVol);
	Img_vol().setXmippOrigin();

	mduntilt_exp.getValue(MDL_IMAGE, Untilted_filename_aux, 1);
	ImgUn_exp.read(Untilted_filename_aux);
	ImgUn_exp.getDimensions(Xdim, Ydim, Zdim, Ndim);

	//Xdim = NSIZE(Img_vol);
	int Xdim_int = (int) Xdim;
	int Ydim_int = (int) Ydim;

	generateFourierStack(allGalleryProjection, galleryTransform_); //Untilted

	size_t idx =0;

	len_u = mduntilt_exp.size();

	position_u_gallery_and_psi.initZeros(7,len_u);

	FOR_ALL_OBJECTS_IN_METADATA(mduntilt_exp)
	{
		mduntilt_exp.getValue(MDL_IMAGE, fnuntilt_exp, __iter.objId);
		Untilted_filenames.push_back(fnuntilt_exp);

		ImgUn_exp.read(fnuntilt_exp);	//Reading image
		ImgUn_exp().setXmippOrigin();

		ImgUn_exp_copy = ImgUn_exp();

		std::cout << "-----------------" <<std::endl;
		std::cout << fnuntilt_exp << "  " << std::endl;
		std::cout << "-----------------" <<std::endl;

		double corr=0, bestcorr = 0;

		for (size_t j = 1; j< len_p; j++)
		{
			//std::cout << "j = " << j << std::endl;
			imgGallery.aliasImageInStack(allGalleryProjection,j);
			imgGallery_orig = imgGallery;
			imgGallery_orig.setXmippOrigin();
			ImgUn_exp_copy2 = ImgUn_exp();

//			CORRELATION UNTILT AND PROJECTIONS
			corr = alignImages(imgGallery_orig, galleryTransform_[j], ImgUn_exp_copy2, transformation_matrix, true, aux_u, aux2, aux3);

			if ((fabs(MAT_ELEM(transformation_matrix, 0, 2)) > maxshift) || (fabs(MAT_ELEM(transformation_matrix, 1, 2)) > maxshift))
				continue;

			if ((corr <0.7*bestcorr) || (corr<0))
				continue;

//			//UNTILT ASSIGNMENT
			double psi = atan2( MAT_ELEM(transformation_matrix,1,0), MAT_ELEM(transformation_matrix,0,0) )*180/PI;

			if (corr>bestcorr)
			{
//					std::cout << "j = " << j << std::endl;
				#ifdef DEBUG
				std::cout << "j = " << j << std::endl;
				std::cout << fnuntilt_exp << "  " << std::endl;
				//UNTILTED PARTICLE
				Image<double> save;
				save()=ImgUn_exp();
				save.write("PPPExp_img.xmp");  //Experimental image
				save()=ImgUn_exp_copy2;
				save.write("PPPExp_img_rot.xmp");  //Experimental rotated image
				save()=imgGallery_orig;
				save.write("PPPGal_img.xmp");  //Gallery image

				std::cout << "rot = " << MAT_ELEM(angles_rot_tilt, j, 0) << std::endl;
				std::cout << "tilt = " << MAT_ELEM(angles_rot_tilt, j, 1) << std::endl;
				std::cout << "psi = " << psi << std::endl;
				std::cout << "corr = " << corr << std::endl;

				char c;
				getchar();
				#endif
				bestcorr = corr;
				MAT_ELEM(position_u_gallery_and_psi, 0, idx) = j;
				MAT_ELEM(position_u_gallery_and_psi, 1, idx) = MAT_ELEM(angles_rot_tilt, j, 0);
				MAT_ELEM(position_u_gallery_and_psi, 2, idx) = MAT_ELEM(angles_rot_tilt, j, 1);
				MAT_ELEM(position_u_gallery_and_psi, 3, idx) = psi;
				MAT_ELEM(position_u_gallery_and_psi, 4, idx) = corr;
				MAT_ELEM(position_u_gallery_and_psi, 5, idx) = -MAT_ELEM(transformation_matrix,0,2);  // X Shift  in untilted image
				MAT_ELEM(position_u_gallery_and_psi, 6, idx) = -MAT_ELEM(transformation_matrix,1,2);  // Y Shift

			}
		}
		++idx;
	}

	MetaData mduntilt_output;
	size_t objId_un;

	for (size_t k=0; k<(len_u); k++)
	{
		//std::cout << "Iteration = " << k << std::endl;

		objId_un = mduntilt_output.addObject();
		//std::cout << "ObjeId_un = " << objId_un << std::endl;
		mduntilt_output.setValue(MDL_IMAGE, Untilted_filenames[k], objId_un);
		mduntilt_output.setValue(MDL_PARTICLE_ID, k+1, objId_un);
		mduntilt_output.setValue(MDL_ANGLE_ROT, MAT_ELEM(position_u_gallery_and_psi, 1, k), objId_un);
		mduntilt_output.setValue(MDL_ANGLE_TILT, MAT_ELEM(position_u_gallery_and_psi, 2, k), objId_un);
		mduntilt_output.setValue(MDL_ANGLE_PSI, MAT_ELEM(position_u_gallery_and_psi, 3, k), objId_un);
		mduntilt_output.setValue(MDL_MAXCC, MAT_ELEM(position_u_gallery_and_psi, 4, k), objId_un);
		mduntilt_output.setValue(MDL_SHIFT_X, MAT_ELEM(position_u_gallery_and_psi, 5, k), objId_un);
		mduntilt_output.setValue(MDL_SHIFT_Y, MAT_ELEM(position_u_gallery_and_psi, 6, k), objId_un);
	}

	mduntilt_output.write(fnDir+"/"+fnun_out);
}


void ProgSolidAngleAssignment::generateProjections(FileName &fnVol, double &smprt)
{
	FileName fnGallery, fnGalleryMetaData;

	// Generate projections
	char gal[300];

	//sprintf(gal, "%s/%s",fnDir.c_str(), fngallerystk.c_str());
	fnGallery=formatString("%s/gallery.stk",fnDir.c_str());

//	String args=formatString("-i %s -o %s --sampling_rate %f",
//			fnVol.c_str(),gal,smprt);
	String args=formatString("-i %s -o %s --sampling_rate %f", fnVol.c_str(),fnGallery.c_str(),smprt);
			//We are considering the psi sampling = angular sampling rate

	std::cout << args << std::endl;
	String cmd=(String)"xmipp_angular_project_library " + args;
	system(cmd.c_str());
}


void ProgSolidAngleAssignment::generateFourierStack(const MultidimArray<double> &input_stack,
															std::vector< AlignmentTransforms> &galleryTransforms_Test)
{
	FourierTransformer transformer;
	AlignmentAux aux2;
	MultidimArray<double> mGalleryProjection;

	// Calculate transforms of this stack
	size_t kmax = NSIZE(input_stack);

	galleryTransforms_Test.resize(kmax);

	for (size_t k=0; k<kmax; ++k)
	{
		mGalleryProjection.aliasImageInStack(input_stack,k);
		mGalleryProjection.setXmippOrigin();

		AlignmentTransforms imageTransforms_Test;
		transformer.FourierTransform(mGalleryProjection, imageTransforms_Test.FFTI, true);
		size_t Ydim, Xdim, Zdim, Ndim;
		imageTransforms_Test.FFTI.getDimensions(Xdim, Ydim, Zdim, Ndim);
		normalizedPolarFourierTransform(mGalleryProjection, imageTransforms_Test.polarFourierI, false,
										XSIZE(mGalleryProjection) / 5, XSIZE(mGalleryProjection) / 2, aux2.plans, 1);

		galleryTransforms_Test[k] = imageTransforms_Test;
	}
}


void ProgSolidAngleAssignment::run()
{
	std::cout << "Starting..." << std::endl;

	MetaData md_exp, md_ang_assign;
	FileName fnAng_assignment = "Angular_Assignment.xmd";

	md_exp.read(fnparticles);

	std::cout << "Performing Angular Assignment" << std::endl;
	assignAngles(md_exp, fnAng_assignment);

	md_ang_assign.read(fnAng_assignment);


	//Particles_with same direction, splitting in groups.




}
