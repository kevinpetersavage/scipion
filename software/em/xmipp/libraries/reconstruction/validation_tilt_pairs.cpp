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
#include <data/matrix1d.h>
#include <data/metadata_extension.h>
#include <complex>
#include <cmath>
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
	addParamsLine(" --untilt <md_file=\"\">     : Metadata with angular assignment for the untilted images");
	addParamsLine(" --tilt <md_file=\"\">       : Metadata with angular assignment for the tilted images");
	addParamsLine("  [--sym <symfile=c1>]       : Enforce symmetry in projections");
	addParamsLine(" --odir <md_file=\"\">        : Metadata with matrix transformation");
	addParamsLine("  [--vol <img_file=\"\">]    : Input reference volume");
}


//Read params
void ProgValidationTiltPairs::readParams()
{
	fnuntilt = getParam("--untilt");
	fntilt = getParam("--tilt");  //Set of tilted coordinates
	fnSym = getParam("--sym");
    fnOut = getParam("--odir");  //Output file
}


void ProgValidationTiltPairs::validate(MetaData &md_u, MetaData &md_t, MetaData &md_out, MetaData &md_validation)
{
	MetaData md_u_sorted, md_t_sorted;
	size_t len_u, len_t, elementId_u, elementId_t, objId_out, objId_out2;
	double rot_u, rot_t, tilt_u, tilt_t, psi_u, psi_t, aux_axis_y, alpha, beta, gamma, tr, rotated_angle, eigenvalueapprox;
	Matrix2D<double> ZYZ_u, ZYZ_t, ZYZ_angles, B_aux, P_mat;
	Matrix1D<double> Eu_dir, eigenvector, angs_axis, axis;

	axis.initZeros(3);

	len_u = md_u.size();  //Metadata length
	len_t = md_t.size();

	if (len_u!=len_t)
	{
		std::cerr << "ERROR: Mismatch dimensions: Number of untilted and tilted particles is not the same" << std::endl;
		exit(0);
	}

	md_u_sorted.sort(md_u, MDL_PARTICLE_ID, true);
	md_t_sorted.sort(md_t, MDL_PARTICLE_ID, true);

	angs_axis.initZeros(3);

	ZYZ_u.initZeros(4,4);
	ZYZ_t.initZeros(4,4);
	ZYZ_angles.initZeros(4,4);

	FOR_ALL_OBJECTS_IN_METADATA2(md_u_sorted, md_t_sorted)
	{
		md_u_sorted.getValue(MDL_PARTICLE_ID, elementId_u ,__iter.objId);
		md_u_sorted.getValue(MDL_ANGLE_ROT, rot_u ,__iter.objId);
		md_u_sorted.getValue(MDL_ANGLE_TILT, tilt_u,__iter.objId);
		md_u_sorted.getValue(MDL_ANGLE_PSI, psi_u,__iter.objId);

		md_t_sorted.getValue(MDL_PARTICLE_ID, elementId_t,__iter2.objId);
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

		eigenvalueapprox = 0.99;

		powerIterationMethod(ZYZ_angles, eigenvector, eigenvalueapprox);

		VEC_ELEM(axis,0) = VEC_ELEM(eigenvector,0);
		VEC_ELEM(axis,2) = VEC_ELEM(eigenvector,2);

		aux_axis_y = VEC_ELEM(eigenvector,1);
		if (aux_axis_y<0){
			VEC_ELEM(axis,1) = -aux_axis_y;}
		else{
			VEC_ELEM(axis,1) = aux_axis_y;
		}

		#ifdef DEBUG
		std::cout << "alpha = " << alpha << std::endl;
		std::cout << "betaa = " << beta << std::endl;
		std::cout << "gamma = " << gamma << std::endl;
		std::cout << "        " << std::endl;
		std::cout << "        " << std::endl;
		#endif

		objId_out = md_out.addObject();
		md_out.setValue(MDL_PARTICLE_ID, elementId_u,objId_out);
		md_out.setValue(MDL_ANGLE_ROT, alpha, objId_out);
		md_out.setValue(MDL_ANGLE_TILT, beta, objId_out);
		md_out.setValue(MDL_ANGLE_PSI, gamma, objId_out);

		objId_out2 = md_validation.addObject();
		md_validation.setValue(MDL_PARTICLE_ID, elementId_u,objId_out2);
		md_validation.setValue(MDL_ROTATION_ANGLE, rotated_angle, objId_out2);
		md_validation.setValue(MDL_TILT_AXIS_X, VEC_ELEM(axis,0), objId_out2);
		md_validation.setValue(MDL_TILT_AXIS_Y, VEC_ELEM(axis,1), objId_out2);
		md_validation.setValue(MDL_TILT_AXIS_Z, VEC_ELEM(axis,2), objId_out2);
	}
	md_out.write((String)"particles@"+fnOut+"/_angular_assignment"+".xmd");
	md_validation.write((String)"particles@"+fnOut+"/_angular_validation"+".xmd");
}


void ProgValidationTiltPairs::powerIterationMethod(const Matrix2D<double> A, Matrix1D<double> &eigenvector, double &eigenvalueapprox)
{
	Matrix2D<double> eye, aux;
	Matrix1D<double> vec_seed, c, b, b_check;
	double Cte_norm, threshold = 0.0001;
	int count = 0;

	eye.initIdentity(A.mdimx);
	b.initZeros(A.mdimx);
	c.initZeros(A.mdimx);

	randomize_random_generator();
	VEC_ELEM(b,0) = double(std::rand())/RAND_MAX;
	VEC_ELEM(b,1) = double(std::rand())/RAND_MAX;
	VEC_ELEM(b,2) = sqrt(1 - VEC_ELEM(b,0)*VEC_ELEM(b,0) - VEC_ELEM(b,1)*VEC_ELEM(b,1));

	//std::cout << "Seed vector = " << b << std::endl;

	b_check = (A*b-b);

//	std::cout << "b_check" << b_check << std::endl;
//	std::cout << "A*b" << A*b << std::endl;

	while ((fabs(VEC_ELEM(b_check, 0))>threshold) || (fabs(VEC_ELEM(b_check, 1))>threshold) || (fabs(VEC_ELEM(b_check, 2))>threshold))
	{
		aux = A-eye*eigenvalueapprox;

		c = aux.inv()* b;
		Cte_norm = c.module();
		c = c/Cte_norm;
		b = c;
		b_check = (A*b-b);
	}
	eigenvector = b;
	//std::cout << "eigenvector = " << eigenvector << std::endl;

	if ((isnan(VEC_ELEM(eigenvector, 0))) || (isnan(VEC_ELEM(eigenvector, 1))) || (isnan(VEC_ELEM(eigenvector, 2))))
	{
		//std::cout << "Recalculating" << std::endl;
		powerIterationMethod(A, eigenvector, eigenvalueapprox);
		count++;
	}

	if (count++>100)
	{
		std::cout << "Eigenvector could not be resolved" << std::endl;
		exit(0);
	}
}


double ProgValidationTiltPairs::R0_ThresholdFisher(int N)
{
	double R0=sqrt(N*3.782);
	return R0;
}

void ProgValidationTiltPairs::R_Fisher(MetaData md, MetaData &md_scattering)
{
	double SumRx, SumRy, SumRz, R_aux, R, theta;
	size_t objId_scttr;
	Matrix1D<double> orth_vector, axis_mean, axis_aux, rot_vector;
	Matrix2D<double> Rodrigues_mat;
	orth_vector.initZeros(3);
	rot_vector.initZeros(3);

	R_aux = R_Value(md, SumRx, SumRy, SumRz);

	axis_mean.initZeros(3);
	axis_aux.initZeros(3);

	VEC_ELEM(axis_mean,0) = SumRx*(1/R_aux);
	VEC_ELEM(axis_mean,1) = SumRy*(1/R_aux);
	VEC_ELEM(axis_mean,2) = SumRz*(1/R_aux);

	FOR_ALL_OBJECTS_IN_METADATA(md)
	{
		md.getValue(MDL_ROTATION_ANGLE, theta, __iter.objId);
		md.getValue(MDL_TILT_AXIS_X, VEC_ELEM(axis_aux,0), __iter.objId);
		md.getValue(MDL_TILT_AXIS_Y, VEC_ELEM(axis_aux,1), __iter.objId);
		md.getValue(MDL_TILT_AXIS_Z, VEC_ELEM(axis_aux,2), __iter.objId);

		orth_vector = vectorProduct(axis_mean, axis_aux);

		Rodrigues_mat = RodriguesMatrix(theta, axis_aux);
		rot_vector = Rodrigues_mat*orth_vector;

		objId_scttr = md_scattering.addObject();
		md_scattering.setValue(MDL_TILT_AXIS_X, VEC_ELEM(rot_vector,0), objId_scttr);
		md_scattering.setValue(MDL_TILT_AXIS_Y, VEC_ELEM(rot_vector,1), objId_scttr);
		md_scattering.setValue(MDL_TILT_AXIS_Z, VEC_ELEM(rot_vector,2), objId_scttr);
	}
}

Matrix2D<double> ProgValidationTiltPairs::RodriguesMatrix( const double theta, const Matrix1D<double> axis)
{
	double k1 = VEC_ELEM(axis,0);
	double k2 = VEC_ELEM(axis,1);
	double k3 = VEC_ELEM(axis,2);

	Matrix2D<double> K, Id, R;
	Id.initIdentity(3);

	K.initZeros(3,3);

	MAT_ELEM(K, 0, 1) = -k3;
	MAT_ELEM(K, 1, 0) = k3;
	MAT_ELEM(K, 0, 2) = k2;
	MAT_ELEM(K, 2, 0) = -k2;
	MAT_ELEM(K, 1, 2) = -k1;
	MAT_ELEM(K, 2, 1) = k1;

	double theta_rad = theta*PI/180;
	R = Id + sin(theta_rad)*K + (1-cos(theta_rad))*K*K;

	return R;
}

double ProgValidationTiltPairs::R_Value(MetaData md, double &SumRx, double &SumRy, double &SumRz)
{
	//double R_fisher = ;
	SumRx = 0;
	SumRy = 0;
	SumRz = 0;

	double Rx, Ry, Rz;

	FOR_ALL_OBJECTS_IN_METADATA(md)
	{
		md.getValue(MDL_TILT_AXIS_X, Rx, __iter.objId);
		md.getValue(MDL_TILT_AXIS_Y, Ry, __iter.objId);
		md.getValue(MDL_TILT_AXIS_Z, Rz, __iter.objId);

		SumRx = SumRx + Rx;
		SumRy = SumRy + Ry;
		SumRz = SumRz + Rz;
	}
	double R = sqrt(SumRx*SumRx + SumRy*SumRy + SumRz*SumRz);
	return R;
}


double ProgValidationTiltPairs::confidenceSolidAngle(double R, int Nparticles, double confidence)
{
	double aux1 = pow((1/confidence),((1/Nparticles)-1)) - 1;
	double aux2 = (Nparticles-R)/R;
	double alpha_th = acos(1-aux2*aux1);

	return alpha_th;
}



void ProgValidationTiltPairs::run()
{
	std::cout << "Starting..." << std::endl;
	MetaData md_u, md_t, md_mic, md_out, md_val, md_scattering;

	Matrix2D<double> M;
	Matrix1D<double> eigenvector;

	md_u.read(fnuntilt);
	md_t.read(fntilt);

	validate(md_u, md_t, md_out, md_val);

	R_Fisher(md_val, md_scattering);

	FileName fn_scattering = formatString("particles@%s/scatter_md.xmd", fnOut.c_str());
	FileName fn_val = formatString("particles@%s/validate_md.xmd", fnOut.c_str());
	md_scattering.write(fn_scattering);
	md_val.write(fn_scattering);


	std::cout << "Finished!" << std::endl;
}
