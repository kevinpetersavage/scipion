/***************************************************************************
 *
 * Authors:    Jose Luis Vilas          (jlvilas@cnb.csic.es)
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
#include "validation_tilt_pairs.h"
#include <data/metadata.h>
#include <data/matrix2d.h>
#include <data/metadata_extension.h>
#include <complex>
//#include <external/alglib/src/linalg.h>
//#define DEBUG
//#define DEBUG1

/*
 * xmipp_validation_tilt_pairs --tilt /home/vilas/ScipionUserData/projects/rct/Runs/001623_XmippProtValidateTilt/extra/tilted/angles_iter001_00.xmd --untilt /home/vilas/ScipionUserData/projects/rct/Runs/001623_XmippProtValidateTilt/extra/untilted/angles_iter001_00.xmd -o caca
 * */

//Define Program parameters
void ProgValidationTiltPairs::defineParams()
{
    //Usage
    addUsageLine("Takes two coordinates sets and determines the transformation between them");
	addUsageLine("First set defines the untilted coordinates, second set defines the tilted coordinates");
	addParamsLine(" --tilt <md_file=\"\"> : Metadata with angular assignment for the tilted images");
	addParamsLine(" --untilt <md_file=\"\"> : Metadata with angular assignment for the untilted images");

	addParamsLine(" --untilt2assign <md_file=\"\"> : Metadata with untilted particles to be angular assigned");
	addParamsLine(" --tilt2assign <md_file=\"\"> : Metadata with tilted particles to be angular assigned");

	addParamsLine(" -o <md_file=\"\"> : Metadata with matrix transformation");

	addParamsLine("  [--vol <img_file=\"\">]    : Input reference volume");
	addParamsLine("  [--angular_sampling <s=5>]   : Angular sampling rate in degrees. This sampling "
			"represents the accuracy for assigning directions");
	addParamsLine("  [--maxshift <s=10>]   : Maximum shift for aligning images (in pixels)");
}


//Read params
void ProgValidationTiltPairs::readParams()
{
    fntilt = getParam("--tilt");  //Set of tilted coordinates
    fnuntilt = getParam("--untilt");
    fnuntilt2assign = getParam("--untilt2assign");
    fntilt2assign = getParam("--tilt2assign");
	fnOut = getParam("-o");  //Output file
	fnVol = getParam("--vol");
	smprt = getDoubleParam("--angular_sampling");
	maxshift = getIntParam("--maxshift");
}

void ProgValidationTiltPairs::generateFourierStackTP(const MultidimArray<double> &input_stack,
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


void ProgValidationTiltPairs::assignAngles(const MetaData mduntilt_exp, FileName fnun_out)
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

	FileName fnprojection = formatString("%s/gallery.doc",fnOut.c_str());
	mdproj.read(fnprojection);
	len_p = mdproj.size();


	imgstack.read(fnOut+"/gallery.stk");

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

	generateFourierStackTP(allGalleryProjection, galleryTransform_); //Untilted

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
			#ifdef DEBUG1
			std::cout << "rot = " << MAT_ELEM(angles_rot_tilt, j, 0) << std::endl;
			std::cout << "tilt = " << MAT_ELEM(angles_rot_tilt, j, 1) << std::endl;
			std::cout << "psi = " << psi << std::endl;
			std::cout << "corr = " << corr << std::endl;
			char ch;
			ch = getchar();
			#endif
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

	mduntilt_output.write(fnOut+"/"+fnun_out);
}

void ProgValidationTiltPairs::validate(MetaData &md_u, MetaData &md_t, MetaData &md_out, MetaData &md_validation)
{
	size_t len_u, len_t;

	len_u = md_u.size();  //Metadata length
	len_t = md_t.size();

	if (len_u!=len_t)
	{
		std::cerr << "ERROR: Mismatch dimensions: Number of untilted and tilted particles is not the same" << std::endl;
		exit(0);
	}

	MetaData md_u_sorted, md_t_sorted;
	md_u_sorted.sort(md_u, MDL_ITEM_ID, true);
	md_t_sorted.sort(md_t, MDL_ITEM_ID, true);

	size_t elementId_u, elementId_t, objId_out;
	double rot_u, rot_t, tilt_u, tilt_t, psi_u, psi_t;
	Matrix2D<double> ZYZ_u, ZYZ_t, ZYZ_angles, B_aux, P_mat;
	Matrix1D<double> Eu_dir, D_mat;


	ZYZ_u.initZeros(4,4);
	ZYZ_t.initZeros(4,4);
	ZYZ_angles.initZeros(4,4);

	double alpha, beta, gamma, tr, rotated_angle;

	FOR_ALL_OBJECTS_IN_METADATA2(md_u_sorted, md_t_sorted)
	{
		md_u_sorted.getValue(MDL_ITEM_ID, elementId_u ,__iter.objId);
		md_u_sorted.getValue(MDL_ANGLE_ROT, rot_u ,__iter.objId);
		md_u_sorted.getValue(MDL_ANGLE_TILT, tilt_u,__iter.objId);
		md_u_sorted.getValue(MDL_ANGLE_PSI, psi_u,__iter.objId);

		md_t_sorted.getValue(MDL_ITEM_ID, elementId_t,__iter2.objId);
		md_t_sorted.getValue(MDL_ANGLE_ROT, rot_t,__iter2.objId);
		md_t_sorted.getValue(MDL_ANGLE_TILT, tilt_t,__iter2.objId);
		md_t_sorted.getValue(MDL_ANGLE_PSI, psi_t,__iter2.objId);

		Euler_angles2matrix(rot_u, tilt_u, psi_u, ZYZ_u);
		Euler_angles2matrix(rot_t, tilt_t, psi_t, ZYZ_t);

		ZYZ_angles = ZYZ_t* (ZYZ_u.inv());

		#ifdef DEBUG
		std::cout << elementId_u << "  " << elementId_t << std::endl;
		std::cout << "rot_u  " << rot_u << std::endl;
		std::cout << "tilt_u " << tilt_u << std::endl;
		std::cout << "psi_u  " << psi_u << std::endl;
		std::cout << "--------------------" << std::endl;
		std::cout << "rot_t  " << rot_t << std::endl;
		std::cout << "tilt_t " << tilt_t << std::endl;
		std::cout << "psi_t  " << psi_t << std::endl;
		std::cout << "--------------------" << std::endl;
		std::cout << "ZYZ_u = " << ZYZ_u << std::endl;
		std::cout << "ZYZ_t = " << ZYZ_t << std::endl;
		std::cout << "ZYZ_matrix = " << ZYZ_angles << std::endl;
		#endif

		Euler_matrix2angles(ZYZ_angles, alpha, beta, gamma);
		alpha *= -1;

		if (beta>90)
		{
			beta = beta - 90;
			alpha = alpha - 180;
			gamma = gamma + 180;
		}
		tr = MAT_ELEM(ZYZ_angles, 0,0) + MAT_ELEM(ZYZ_angles, 1,1) + MAT_ELEM(ZYZ_angles, 2,2) + MAT_ELEM(ZYZ_angles, 3,3);
		rotated_angle = acos((tr-1)*0.5)*180/PI;

		B_aux.initIdentity(4);

		generalizedEigs(ZYZ_angles, B_aux, D_mat, P_mat);

		std::cout << "Autovalores = " << D_mat << std::endl;
		std::cout << "Autovectores = " << D_mat << std::endl;

//	    if isreal(eig_angles(1,1))
//	        control(k)=1;
//	        ejevector(:,k) = eig_vectors(:,1);
//	        eje_ang_x(k) = acosd(ex*eig_vectors(:,1));
//	        eje_ang_y(k) = acosd(ey*eig_vectors(:,1));
//	        eje_ang_z(k) = acosd(ez*eig_vectors(:,1));
//	    else
//	        if isreal(eig_angles(2,2))
//	            control(k)=2;
//	            ejevector(:,k) = eig_vectors(:,2);
//	            eje_ang_x(k) = acosd(ex*eig_vectors(:,2));
//	            eje_ang_y(k) = acosd(ey*eig_vectors(:,2));
//	            eje_ang_z(k) = acosd(ez*eig_vectors(:,2));
//	        else
//	            if isreal(eig_angles(3,3))
//	                control(k)=3;
//	                ejevector(:,k) = eig_vectors(:,3);
//	                eje_ang_x(k) = acosd(ex*eig_vectors(:,3));
//	                eje_ang_y(k) = acosd(ey*eig_vectors(:,3));
//	                eje_ang_z(k) = acosd(ez*eig_vectors(:,3));
//	            else
//	                count_error = count_error + 1;
//	            end
//	        end
//	    end



		#ifdef DEBUG
		std::cout << "alpha = " << alpha << std::endl;
		std::cout << "betaa = " << beta << std::endl;
		std::cout << "gamma = " << gamma << std::endl;
		std::cout << "        " << std::endl;
		std::cout << "        " << std::endl;
		#endif

		objId_out = md_out.addObject();
		md_out.setValue(MDL_ITEM_ID, elementId_u,objId_out);
		md_out.setValue(MDL_ANGLE_ROT, alpha, objId_out);
		md_out.setValue(MDL_ANGLE_TILT, beta, objId_out);
		md_out.setValue(MDL_ANGLE_PSI, gamma, objId_out);
		/*
		objId_out = md_out.addObject();
		md_out.setValue(MDL_ITEM_ID, elementId_u,objId_out);
		md_out.setValue(MDL_ROTAION_ANGLE, alpha, objId_out);
		md_out.setValue(MDL_TILT_AXIS_X, alpha, objId_out);
		md_out.setValue(MDL_TILT_AXIS_X, beta, objId_out);
		md_out.setValue(MDL_TILT_AXIS_X, gamma, objId_out);
		*/
	}
	md_out.write((String)"particles@"+fnOut+"_angular_assignment"+".xmd");
	md_validation.write((String)"particles@"+fnOut+"_angular_validation"+".xmd");
}


void ProgValidationTiltPairs::powerIterationMethod(const Matrix2D<double> M)
{
	Matrix1D <double> seed, c;
    double d=0,temp;
    int i,j;

    size_t n =M.mdimx;
    //x is the input vector

    c.initZeros(n);
    seed.initZeros(n);

    VEC_ELEM(seed,0) = 0;
    VEC_ELEM(seed,1) = 1;
    VEC_ELEM(seed,2) = 0;

    do
    {
        for(i=0;i<n;i++)
        {
        	VEC_ELEM(c,i)=0;
            for(j=0;j<n;j++)
            	VEC_ELEM(c,i)+=MAT_ELEM(M,i,j)*VEC_ELEM(seed,j);
        }
        for(i=0;i<n;i++)
        	VEC_ELEM(seed,i)=VEC_ELEM(c,i);

        temp=d;
        d=0;

        for(i=0;i<n;i++)
        {
            if(fabs(VEC_ELEM(seed,i))>fabs(d))
                d=VEC_ELEM(seed,i);
        }
        for(i=0;i<n;i++)
        	VEC_ELEM(seed,i)/=d;

    }while(fabs(d-temp)>0.00001);

    std::cout<<"Eigen value is : "<< d <<std::endl;

    std::cout<<"Eigenvector is: \n";
    for(i=0;i<n;i++)
    	std::cout<<VEC_ELEM(seed,i)<<std::endl;

}



void ProgValidationTiltPairs::run()
{
	std::cout << "Starting..." << std::endl;
	MetaData md_u, md_t, md_mic, md_out, md_val;;

	Matrix2D<double> M;
	M.initZeros(3,3);

	MAT_ELEM(M,0,0) = 0.6056;
	MAT_ELEM(M,0,1) = -0.3741;
	MAT_ELEM(M,0,2) = -0.7023;
	MAT_ELEM(M,1,0) = -0.3110;
	MAT_ELEM(M,1,1) = 0.7012;
	MAT_ELEM(M,1,2) = -0.6416;
	MAT_ELEM(M,2,0) = 0.7325;
	MAT_ELEM(M,2,1) = 0.6070;
	MAT_ELEM(M,2,2) = 0.3083;

	powerIterationMethod(M);


//	if (fnuntilt2assign != "" &&  fntilt2assign != "")
//	{
//		MetaData mduntilt_exp, mdtilt_exp;
//		mduntilt_exp.read(fnuntilt2assign);
//		mdtilt_exp.read(fntilt2assign);
//
//		std::cout << "Performing Untilted Assigning" << std::endl;
//		assignAngles(mduntilt_exp, "untilted_assigned.xmd");
//
//		std::cout << "Performing Tilted Assigning" << std::endl;
//		assignAngles(mdtilt_exp, "tilted_assigned.xmd");
//
//		md_u.read(fnuntilt);
//		md_t.read(fntilt);
//		validate(md_u, md_t, md_out, md_val);
//
////		rmatrixevd(const real_2d_array &a, const ae_int_t n, const ae_int_t vneeded, real_1d_array &wr, real_1d_array &wi, ...
////				real_2d_array &vl, real_2d_array &vr);
//
////		FOR_ALL_OBJECTS_IN_METADATA(md_out)
////		{
////			md_u_sorted.getValue(MDL_ITEM_ID, elementId_u ,__iter.objId);
////			md_u_sorted.getValue(MDL_ANGLE_ROT, rot_u ,__iter.objId);
////			md_u_sorted.getValue(MDL_ANGLE_TILT, tilt_u,__iter.objId);
////			md_u_sorted.getValue(MDL_ANGLE_PSI, psi_u,__iter.objId);
////
////			md_t_sorted.getValue(MDL_ITEM_ID, elementId_t,__iter2.objId);
////			md_t_sorted.getValue(MDL_ANGLE_ROT, rot_t,__iter2.objId);
////			md_t_sorted.getValue(MDL_ANGLE_TILT, tilt_t,__iter2.objId);
////			md_t_sorted.getValue(MDL_ANGLE_PSI, psi_t,__iter2.objId);
////		}
//	}
//	else
//	{
//		md_u.read(fnuntilt);
//		md_t.read(fntilt);
//
//
//		validate(md_u, md_t, md_out,md_val);
//	}






}


