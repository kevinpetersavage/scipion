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
//#include <data/filters.h>
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
#include "morphology.h"
//#define DEBUG


void ProgAngularTiltPairAssignment::readParams()
{
	fnUntilt = getParam("--untiltparticles");
	fnTilt = getParam("--tiltparticles");
	alphaT = getDoubleParam("--alphaT");
	tilt_mic = getDoubleParam("--tiltmic");
	alphaU = getDoubleParam("--alphaU");
	fnDir = getParam("--odir");
	fnOutUntiltedAssign = getParam("--untiltassignment");
	fnOutTiltedAssign =	getParam("--tiltassignment");
	fnSym = getParam("--sym");
	smprt = getDoubleParam("--angular_sampling");
	maxshift = getDoubleParam("--maxshift");
	fnVol = getParam("--initvol");
	pad = getIntParam("--pad");
	fnmic = getParam("--micangles");
}
ProgAngularTiltPairAssignment::ProgAngularTiltPairAssignment()
{
	rank=0;
	Nprocessors=1;
}

void ProgAngularTiltPairAssignment::defineParams()
{
	//usage
	addUsageLine("From an untilted particles stack and an initial volume, an angular assignment of the untilted particles is performed."
			"Then, knowing the tilt axis and the tilt angle of the micrography, tiltmic, and the rotational angles of the untilted and"
			"tilted micrographs, alphaU and alphaT respectively, an angular assignment of the tilted particles is determined");
	//params
	addParamsLine("  [--untiltparticles <md_file=\"\">]    : Untilt particles stack");
	addParamsLine("  [--tiltparticles <md_file=\"\">]    : Tilt particles stack");
	addParamsLine("  [--alphaT <s=0>]   : Rotational angle of the tilted micrograph (degrees)");
	addParamsLine("  [--tiltmic <s=0>]   : Tilt angle of the micrograph (degrees)");
	addParamsLine("  [--alphaU <s=0>]   : Rotational angle of the untilted micrograph (degrees)");
	addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
	addParamsLine("  [--untiltassignment <outputDir=\".\">]   : Name of the metadata file for the untilted assignment");
	addParamsLine("  [--tiltassignment <outputDir=\".\">]   : Name of the metadata file for the tilted assignment");
	addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
	addParamsLine("  [--angular_sampling <s=5>]   : Angular sampling rate in degrees. This sampling represents the accuracy for assigning directions");
	addParamsLine("  [--maxshift <s=10>]   : Maximum shift for aligning images (in pixels)");
	addParamsLine("  [--initvol <md_file=\"\">]   : Initial Volume");
	addParamsLine("  [--pad <s=10>]   : Padding factor");
	addParamsLine("  [--micangles <md_file=\"\">]   : Metadata with micrograph angles, rot, tilt, psi and the micrograph_Id");
}


void ProgAngularTiltPairAssignment::generateProjections(FileName &fnVol, double &smprt)
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


void ProgAngularTiltPairAssignment::generateFourierStack(const MultidimArray<double> &input_stack,
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


void ProgAngularTiltPairAssignment::generateInitialBall(const MetaData &md_u,const MetaData &md_t,
																 MetaData &md_u_assign_iter0, MetaData &md_t_assign_iter0, FileName &fnVol)
{
	randomize_random_generator();
	FileName vol_iter, img_fn_u, img_fn_t;
	if (rank==0)
		std::cout << "Generating initial volume" << std::endl;
	size_t len_u =md_u.size();


	double rand_u, rand_v, rand_theta, rand_phi, rand_psi, rot_t, tilt_t, psi_t;
	Matrix2D<double> ZYZ_u, ZYZ_angles, ZYZ_t;

	for (size_t k=0; k<len_u; k++)
	{
		rand_u = double(std::rand())/RAND_MAX;
		rand_v = double(std::rand())/RAND_MAX;
		rand_theta = 360*rand_u;
		rand_phi =  acos(2*rand_v-1)*180/PI;
		rand_psi = (double(std::rand())*180/PI)/RAND_MAX;

		#ifdef DEBUG
		std::cout << "-----------------------" << std::endl;
		std::cout << rand_u << "  " << rand_v << std::endl;
		std::cout << rand_theta << std::endl;
		std::cout << rand_phi << std::endl;
		std::cout << rand_psi << std::endl;
		std::cout << "-----------------------" << std::endl;
		#endif

		size_t objId_un = md_u_assign_iter0.addObject();
		md_u_assign_iter0.setValue(MDL_ITEM_ID, k,objId_un);
		md_u.getValue(MDL_IMAGE,img_fn_u,objId_un);
		md_u_assign_iter0.setValue(MDL_IMAGE, img_fn_u,objId_un);
		md_u_assign_iter0.setValue(MDL_ANGLE_ROT,rand_phi,objId_un);
		md_u_assign_iter0.setValue(MDL_ANGLE_TILT,rand_theta,objId_un);
		md_u_assign_iter0.setValue(MDL_ANGLE_PSI,rand_psi,objId_un);
		md_u_assign_iter0.setValue(MDL_SHIFT_X,(double) 0,objId_un);
		md_u_assign_iter0.setValue(MDL_SHIFT_Y,(double) 0,objId_un);

		Euler_angles2matrix(rand_phi, rand_theta, rand_psi, ZYZ_u);
		Euler_angles2matrix(-alphaU, tilt_mic, alphaT, ZYZ_angles);
		ZYZ_t = ZYZ_angles*ZYZ_u;
		Euler_matrix2angles(ZYZ_t, rot_t, tilt_t, psi_t);

		size_t objId_t = md_t_assign_iter0.addObject();
		md_t_assign_iter0.setValue(MDL_ITEM_ID, k,objId_t);
		md_t.getValue(MDL_IMAGE,img_fn_t,objId_t);
		md_t_assign_iter0.setValue(MDL_IMAGE, img_fn_t, objId_t);
		md_t_assign_iter0.setValue(MDL_ANGLE_ROT,rot_t,objId_t);
		md_t_assign_iter0.setValue(MDL_ANGLE_TILT,tilt_t,objId_t);
		md_t_assign_iter0.setValue(MDL_ANGLE_PSI,psi_t,objId_t);
		md_t_assign_iter0.setValue(MDL_SHIFT_X,(double) 0,objId_t);
		md_t_assign_iter0.setValue(MDL_SHIFT_Y,(double) 0,objId_t);
	}

	FileName filnm_md_u = fnDir+"/angular_assignment_u_0.xmd";
	FileName filnm_md_t = fnDir+"/angular_assignment_t_0.xmd";

	std::cout << "----------------- " << std::endl;
	std::cout << filnm_md_u << std::endl;

	md_u_assign_iter0.write((String) "particles@"+filnm_md_u);
	md_t_assign_iter0.write((String) "particles@"+filnm_md_t);

	String args=formatString("-i %s -o %s/volume_iter_0.vol --sym %s",filnm_md_u.c_str(), fnDir.c_str(),
																	fnSym.c_str());

	String cmd =(String) "xmipp_reconstruct_fourier "+ args+
	" --max_resolution 0.5 --padding 2";
	std::cout << cmd << std::endl;
	system(cmd.c_str());
	fnVol = fnDir+"/volume_iter_0.vol";
}


void ProgAngularTiltPairAssignment::run()
{
	if (rank==0)
		std::cout << "Starting..." << std::endl;

	//Initial Parameters
	MetaData mduntilt_exp, mdtilt_exp, md_u_assign_iter0, md_t_assign_iter0, mdproj, md_mic;
	FileName fngallerystk, fnprojection, Untilted_filename_aux, fnuntilt_exp, fntilt_exp, fnall_md;
	Image<double> ImgUn_exp, ImgT_exp, Img_proj, Img_vol, imgstack;
	MultidimArray< std::complex<double> > FImgT_exp;
	MultidimArray<double> Multidim_exp, Multidim_proj, Multidim_vol, Improj_traslated, imgGallery_orig, imgGallery,
							ImgUn_exp_copy, ImgUn_exp_copy2, ImgT_exp_copy;
	double rot, tilt, rot_t, tilt_t, psi_t, shiftX, shiftY;
	Matrix1D<double> Shiftvec(2);
	Matrix2D<double> angles_rot_tilt, transformation_matrix, ZYZ_u, ZYZ_angles, ZYZ_t;
	size_t len_p, len_u, len_t, idxp, Ydim, Xdim, Zdim, Ndim, mic_id_u, mic_id_t, idx_mpi;
	double Ts = 5;
	double maxResol = 20;

	Projection projection, tilt_proj;
	CorrelationAux aux_corr, aux2, Auxcorr_bestshift;
	AlignmentTransforms imageTransforms_T, test_transform;
	FourierTransformer transformer_T;
	AlignmentAux aux_u;
	RotationalCorrelationAux aux3;
	std::vector<std::string> Untilted_filenames, Tilted_filenames;
	std::vector< AlignmentTransforms> galleryTransform_;
	//TODO: Comprobar el uso de len_t

	mduntilt_exp.read(fnUntilt);
	mdtilt_exp.read(fnTilt);

	if (fnVol =="")
	{
		if (rank == 0)
		{
		generateInitialBall(mduntilt_exp, mdtilt_exp, md_u_assign_iter0, md_t_assign_iter0, fnVol);
		}
	}
	synchronize();

	//Common parameters
	len_u = mduntilt_exp.size();
	len_t = mdtilt_exp.size();
	Untilted_filenames.resize(len_u);
	Tilted_filenames.resize(len_t);

	Matrix2D<double> position_u_gallery_and_psi;
	position_u_gallery_and_psi.initZeros(14,len_u);


	// Generating projection from the volume, fnVol, with an angular sampling rate of smprt degrees
	if (rank==0)
	{
		std::cout << "Generating projections" << std::endl;
	    generateProjections(fnVol, smprt);
		std::cout << "Projections finished" << std::endl;
	}
	synchronize();

	fnprojection = formatString("%s/gallery.doc",fnDir.c_str());
	mdproj.read(fnprojection);
	len_p = mdproj.size();

	imgstack.read(fnDir+"/gallery.stk");

	const MultidimArray<double> &allGalleryProjection = imgstack();
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
	if (rank==0)
		std::cout << "Searching correlations" << std::endl;

	Img_vol.read(fnVol);
	Img_vol().setXmippOrigin();

	mduntilt_exp.getValue(MDL_IMAGE, Untilted_filename_aux, 1);
	ImgUn_exp.read(Untilted_filename_aux);
	ImgUn_exp.getDimensions(Xdim, Ydim, Zdim, Ndim);

	//Xdim = NSIZE(Img_vol);
	int Xdim_int = (int) Xdim;
	int Ydim_int = (int) Ydim;

	//STACKS FOURIER TRANSFORM
	generateFourierStack(allGalleryProjection, galleryTransform_); //Untilted
	FourierProjector *projectr = new FourierProjector(Img_vol(),pad,Ts/maxResol,BSPLINE3);

	size_t idx =0, lastidx=-1;

	if (fnmic !="")
	{
		md_mic.read(fnmic);
	}

	time_t initial, ending;
	time (&initial);

	//synchronize();
	idx_mpi = 0;
	size_t objId_un, objId_t, ensure_assignment=0;
	bool flag_disable =false;

	FOR_ALL_OBJECTS_IN_METADATA2(mduntilt_exp, mdtilt_exp)
	{
		if ((idx_mpi+1)%Nprocessors==rank)
		{

			mduntilt_exp.getValue(MDL_IMAGE, fnuntilt_exp, __iter.objId);
			mduntilt_exp.getValue(MDL_MICROGRAPH_ID, mic_id_u, __iter.objId);

			if (fnmic !="")
			{
				md_mic.getValue(MDL_ANGLE_ROT, alphaU, mic_id_u);
				md_mic.getValue(MDL_ANGLE_TILT, tilt_mic, mic_id_u);
				md_mic.getValue(MDL_ANGLE_PSI, alphaT, mic_id_u);
			}
			Untilted_filenames.push_back(fnuntilt_exp);

			ImgUn_exp.read(fnuntilt_exp);	//Reading image
			ImgUn_exp().setXmippOrigin();
			mdtilt_exp.getValue(MDL_IMAGE, fntilt_exp, __iter2.objId);
			Tilted_filenames.push_back(fntilt_exp);
			ImgT_exp.read(fntilt_exp);	//Reading image
			ImgT_exp().setXmippOrigin();
			transformer_T.FourierTransform(ImgT_exp(),FImgT_exp,true);

			double corr1 = 0, corr2 = 0, corr = 0, bestcorr1 =0, bestcorr = 0;

			if (rank==0)
			{
				std::cout << fnuntilt_exp << std::endl;
			}

				for (size_t j = 0; j< len_p; j++)
				{
					imgGallery.aliasImageInStack(allGalleryProjection,j);
					imgGallery_orig = imgGallery;
					imgGallery_orig.setXmippOrigin();
					ImgUn_exp_copy2 = ImgUn_exp();

					//CORRELATION UNTILT AND PROJECTIONS
					corr1 = alignImages(imgGallery_orig, galleryTransform_[j], ImgUn_exp_copy2, transformation_matrix, true, aux_u, aux2, aux3);

					if ((fabs(MAT_ELEM(transformation_matrix, 0, 2)) > maxshift) || (fabs(MAT_ELEM(transformation_matrix, 1, 2)) > maxshift))
					{
						ensure_assignment++;
						if (ensure_assignment == len_p)
						{
							ensure_assignment = 0;
							flag_disable = true;
						}
						else
						{
							continue;
						}
					}


					if ((corr1 <0.7*bestcorr1) || (corr1<0))
						continue;

					//UNTILT ASSIGNMENT
					double psi = atan2( MAT_ELEM(transformation_matrix,1,0), MAT_ELEM(transformation_matrix,0,0) )*180/PI;

					Euler_angles2matrix(MAT_ELEM(angles_rot_tilt, j, 0), MAT_ELEM(angles_rot_tilt, j, 1), psi, ZYZ_u);
					Euler_angles2matrix(-alphaU, tilt_mic, alphaT, ZYZ_angles);
					ZYZ_t = ZYZ_angles*ZYZ_u;
					Euler_matrix2angles(ZYZ_t, rot_t, tilt_t, psi_t);

					projectVolume( *projectr, projection, Ydim_int, Xdim_int, rot_t, tilt_t, psi_t, NULL);

					corr2 = bestShift(ImgT_exp(), FImgT_exp, projection(),  shiftX, shiftY, Auxcorr_bestshift, NULL, maxshift);

					VEC_ELEM(Shiftvec,0) = shiftX;
					VEC_ELEM(Shiftvec,1) = shiftY;
					Improj_traslated = projection();
					selfTranslate(1, Improj_traslated, Shiftvec);

					corr2 = correlationIndex(ImgT_exp(), Improj_traslated);

					if (corr2<0)
						continue;

					corr = 0.5*(corr1 + corr2);

						if (corr>bestcorr)
						{
							bestcorr = corr;

							if (idx!=lastidx)
							{
								objId_un = mdPartial_u.addObject();
								objId_t = mdPartial_t.addObject();
								lastidx = idx;
							}

							mdPartial_u.setValue(MDL_IMAGE, fnuntilt_exp, objId_un);
							mdPartial_u.setValue(MDL_PARTICLE_ID, idx+1, objId_un);
							mdPartial_u.setValue(MDL_ANGLE_ROT, MAT_ELEM(angles_rot_tilt, j, 0), objId_un);
							mdPartial_u.setValue(MDL_ANGLE_TILT, MAT_ELEM(angles_rot_tilt, j, 1), objId_un);
							mdPartial_u.setValue(MDL_ANGLE_PSI, psi, objId_un);
							mdPartial_u.setValue(MDL_MAXCC, corr1, objId_un);
							mdPartial_u.setValue(MDL_SHIFT_X, -MAT_ELEM(transformation_matrix,0,2), objId_un);
							mdPartial_u.setValue(MDL_SHIFT_Y, -MAT_ELEM(transformation_matrix,1,2), objId_un);
							if (flag_disable == true)
								mdPartial_u.setValue(MDL_ENABLED, 0, objId_un);
							else
								mdPartial_u.setValue(MDL_ENABLED, 1, objId_un);

//							objId_t = mdPartial_t.addObject();
							mdPartial_t.setValue(MDL_IMAGE, fntilt_exp, objId_t);
							mdPartial_t.setValue(MDL_PARTICLE_ID, idx+1, objId_t);
							mdPartial_t.setValue(MDL_ANGLE_ROT, rot_t, objId_t);
							mdPartial_t.setValue(MDL_ANGLE_TILT, tilt_t, objId_t);
							mdPartial_t.setValue(MDL_ANGLE_PSI, psi_t, objId_t);
							mdPartial_t.setValue(MDL_MAXCC, corr2, objId_t);
							mdPartial_t.setValue(MDL_SHIFT_X, -VEC_ELEM(Shiftvec,0), objId_t);
							mdPartial_t.setValue(MDL_SHIFT_Y, -VEC_ELEM(Shiftvec,1), objId_t);
							if (flag_disable == true)
								mdPartial_t.setValue(MDL_ENABLED, 0, objId_un);
							else
								mdPartial_t.setValue(MDL_ENABLED, 1, objId_un);

							flag_disable = false;


							#ifdef DEBUG
												//UNTILTED PARTICLE
													std::cerr << "rank = " << rank << std::endl;
													Image<double> save;
													save()=ImgUn_exp();
													FileName fn_PPPuntilted = formatString("PPPuntilted%i.xmp", (int) rank);
													save.write(fn_PPPuntilted);  //Experimental untilted image

													save()=imgGallery_orig;
													FileName fn_PPPprojection_untilted = formatString("PPPprojection_untilted%i.xmp", (int) rank);
													save.write(fn_PPPprojection_untilted);   //Gallery untilted particle

													save()=ImgUn_exp_copy2;
													FileName fn_PPPuntilted_rotated = formatString("PPPuntilted_rotated%i.xmp", (int) rank);
													save.write(fn_PPPuntilted_rotated);   //Experimental particle rotated

												//TILTED PARTICLE
													save()=ImgT_exp();
													FileName fn_PPPtilted = formatString("PPPtilted%i.xmp", (int) rank);
													save.write(fn_PPPtilted);    //Experimental tilted image

													save()=Improj_traslated;
													FileName fn_PPPtilted_copy = formatString("PPPtilted_copy%i.xmp", (int) rank);
													save.write(fn_PPPtilted_copy);    //Experimental tilted image copy

													save()=projection();
													FileName fn_PPPprojection_tilted = formatString("PPPprojection_tilted%i.xmp", (int) rank);
													save.write(fn_PPPprojection_tilted);  //Galeria después del align

//													std::cout << "Corr1 =" << VEC_ELEM(corr_vec, kk) << std::endl;
//													std::cout << "Corr2 =" << corr2 << std::endl;
													std::cout << "Corr  =" << corr << std::endl;
													std::cout << "                          " << std::endl;

													std::cout << "Untilted ROT_u = " << MAT_ELEM(angles_rot_tilt, kk, 0) << std::endl;
													std::cout << "Untilted TILT_u = " << MAT_ELEM(angles_rot_tilt, kk, 1) << std::endl;
													std::cout << "Untilted PSI_u = " << VEC_ELEM(psi_aux, kk) << std::endl;
													std::cout << "Tilted ROT_t = " << rot_t << std::endl;
													std::cout << "Tilted TILT_t = " << tilt_t << std::endl;
													std::cout << "Tilted PSI_t = " << psi_t << std::endl;

													std::cout << "                          " << std::endl;
													std::cout << "--------------------------" << std::endl;
													std::cout << "Press any key" << std::endl;
													char c;
													c = getchar();
							#endif
					}
				}
		}
		idx++;
		idx_mpi++;
	}

		FileName fn_upartialmd = formatString("particles@%s/Partial_u_angular_assignment%i.xmd", fnDir.c_str(), (int) rank);
		FileName fn_tpartialmd = formatString("particles@%s/Partial_t_angular_assignment%i.xmd", fnDir.c_str(), (int) rank);
		mdPartial_u.write(fn_upartialmd);
		mdPartial_t.write(fn_tpartialmd);
		synchronize();

		delete projectr;
		mdPartial_u.clear();
		mdPartial_t.clear();
		time (&ending);
		double seconds;
		seconds = difftime(initial,ending);

		if (rank==0)
			std::cerr << "Time employed = " << seconds << std::endl;
		MetaData mdu_all, mdt_all, mdu_aux, mdt_aux;

		if (rank ==0 )
		{
			fn_upartialmd = formatString("particles@%s/Partial_u_angular_assignment%i.xmd", fnDir.c_str(), 0);
			mdu_all.read(fn_upartialmd);
			fn_tpartialmd = formatString("particles@%s/Partial_t_angular_assignment%i.xmd", fnDir.c_str(), 0);
			mdt_all.read(fn_tpartialmd);
			for (int rns = 1; rns<Nprocessors; rns++)
			{
				fn_upartialmd = formatString("particles@%s/Partial_u_angular_assignment%i.xmd", fnDir.c_str(), (int) rns);
				mdu_aux.read(fn_upartialmd);
				mdu_all.unionAll(mdu_aux);
				fn_tpartialmd = formatString("particles@%s/Partial_t_angular_assignment%i.xmd", fnDir.c_str(), (int) rns);
				mdt_aux.read(fn_tpartialmd);
				mdt_all.unionAll(mdt_aux);
			}

//		//Storing assignments into output metadata
		MetaData mduntilt_output, mdtilt_output, md_union, md_union_aux;

		FileName filnm_md_u = formatString("particles@%s", fnOutUntiltedAssign.c_str());
		FileName filnm_md_t = formatString("particles@%s", fnOutTiltedAssign.c_str());

		std::cout << "Angular assignment finished" << std::endl;
		md_union_aux = mdt_all;
		md_union = mdu_all;

		md_union_aux.operate(formatString("particleId = particleId + %d", len_u));
		md_union.unionAll(md_union_aux);

		fnall_md = formatString("particles@%s/All_angular_assignment_lastiter.xmd", fnDir.c_str());
		md_union.write(fnall_md);

		mdu_all.write(filnm_md_u);
		mdt_all.write(filnm_md_t);
		}
		synchronize();
}
