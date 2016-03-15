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

#include "angular_tilt_pair_heterogeneity.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <data/xmipp_image.h>
#include <geometry.h>
#include <data/matrix1d.h>
#include <data/matrix2d.h>
//#include <data/transformations.h>
#define DEBUG


void ProgAngularTiltPairHeterogeneity::readParams()
{
	fnUntilt = getParam("--untiltparticles");
	fnTilt = getParam("--tiltparticles");
	fnDir = getParam("--odir");
	ang_dist = getDoubleParam("--angular_distance");
	ang_acc = getDoubleParam("--angular_accuracy");
	N_cls = getIntParam("--classes");
}


void ProgAngularTiltPairHeterogeneity::defineParams()
{
	//usage
	addUsageLine("From an untilted particles stack and an initial volume, an angular assignment of the untilted particles is performed."
			"Then, knowing the tilt axis and the tilt angle of the micrography, tiltmic, and the rotational angles of the untilted and"
			"tilted micrographs, alphaU and alphaT respectively, an angular assignment of the tilted particles is determined");
	//params
	addParamsLine("  [--untiltparticles <md_file=\"\">]    : Untilt particles stack");
	addParamsLine("  [--tiltparticles <md_file=\"\">]    : Tilt particles stack");
	addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
	addParamsLine("  [--angular_distance <s=5>]   : Angular distance");
	addParamsLine("  [--angular_accuracy <s=5>]   : Angular distance");
	addParamsLine("  [--classes <s=2>]   : Number of classes");
}


void ProgAngularTiltPairHeterogeneity::md2matrix2d(const MetaData md, Matrix2D<double> &output_matrix)
{
	size_t len =md.size();

	output_matrix.initZeros(3,len);
	size_t idx =0;
	FOR_ALL_OBJECTS_IN_METADATA(md)
	{
		md.getValue(MDL_ANGLE_ROT, MAT_ELEM(output_matrix, 0, idx), __iter.objId);
		md.getValue(MDL_ANGLE_TILT, MAT_ELEM(output_matrix, 1, idx), __iter.objId);
		md.getValue(MDL_ANGLE_PSI, MAT_ELEM(output_matrix, 2, idx), __iter.objId);
		idx++;
	}
}

void ProgAngularTiltPairHeterogeneity::md2matrix2dclasses(const MetaData md, const MetaData mdangles, Matrix2D<double> &output_matrix)
{
	size_t len =md.size();

	output_matrix.initZeros(3,len);
	size_t idx =0;
	FOR_ALL_OBJECTS_IN_METADATA(md)
	{
		md.getValue(MDL_ANGLE_ROT, MAT_ELEM(output_matrix, 0, idx), __iter.objId);
		md.getValue(MDL_ANGLE_TILT, MAT_ELEM(output_matrix, 1, idx), __iter.objId);
		md.getValue(MDL_ANGLE_PSI, MAT_ELEM(output_matrix, 2, idx), __iter.objId);
		idx++;
	}
}


void ProgAngularTiltPairHeterogeneity::spheresplit(const double ang_acc, MetaData &mdangles, int &number_of_regions)
{
	double number_ang_latt, number_ang_long, n_long, rest;
	size_t obId;

	//Number of regions on the projection sphere
	number_ang_latt = floor((180-ang_acc)/(2*ang_acc))+1;
	n_long = 2;  //if k=0 or k=180, then number_ang_long=0, it needs to set n_long=2 for counting properly
	size_t count = 0;
//((180-ang_acc)/(2*ang_acc))
	//TODO: check int and double in for-params
	for (int k = 0; k<181 ; k=k+2*ang_acc)
	{
		number_ang_long = floor(360*fabs(sin(k*PI/180))/(2*ang_acc));  //+-angular_acc, therefore, 2k-steps
		if (number_ang_long ==0)
		{
			rest = 0;
			count = count + 1;
			obId = mdangles.addObject();
			mdangles.setValue(MDL_ANGLE_TILT, (double) k, obId);
			mdangles.setValue(MDL_ANGLE_ROT, (double) 0, obId);
			mdangles.setValue(MDL_IDX, count, obId);
		}
		else
		{
			rest = (360-number_ang_long*(2*ang_acc))/(number_ang_long);
		}
		n_long = n_long + number_ang_long;

		for (int j =0; j<number_ang_long; j++)
		{
			count = count + 1;
			obId = mdangles.addObject();
			mdangles.setValue(MDL_ANGLE_TILT, (double) k, obId);
			mdangles.setValue(MDL_ANGLE_ROT, (2*ang_acc+rest)*j, obId);
			mdangles.setValue(MDL_IDX, count, obId);
		}
	}
	number_of_regions = n_long;
	std::cout << " number of circles in meridian = " << number_ang_latt << std::endl;
	std::cout << " number of circles in parallel = " << number_ang_long << std::endl;
	std::cout << " number of regions = " << n_long << std::endl;
	std::cout << " number of rows = " << mdangles.size() << std::endl;

	double Ffactor=n_long*sin(ang_acc*PI/360)*sin(ang_acc*PI/360);
	std::cout << " Fill_factor_regs = " << Ffactor << std::endl;
	std::cout << " Fill_factor_rows = " << mdangles.size()*sin(ang_acc*PI/360)*sin(ang_acc*PI/360) << std::endl;
}


void ProgAngularTiltPairHeterogeneity::run()
{
	std::cout << "Starting..." << std::endl;

	MetaData md_untilted_exp, md_tilted_exp, mdangles, md_untilted_exp_sorted, md_tilted_exp_sorted, mdu_aux, md_all, md_all_sorted;
	Image<double> img1, img2, img_out;
	double rot_u, tilt_u, psi_u, rot_t, tilt_t, psi_t, rot_gal, tilt_gal, corr;
	Matrix2D<double> angles_u, angles_t, transformation_matrix;
	size_t len_u, len_t, idx_u, idx_t;
	int number_of_regions;
	FileName fnallassignSph, fnimg, fnimg_corr1, fnimg_corr2;
	std::vector<FileName> adressesfn;
	AlignmentAux aux;
	CorrelationAux aux2;
	RotationalCorrelationAux aux3;
	AlignmentTransforms galleryTransform_;

	//Reading MD and storing them into a matrix
	md_untilted_exp.read(fnUntilt);
	md_tilted_exp.read(fnTilt);

	len_u = md_untilted_exp.size();
	len_t = md_tilted_exp.size();

	md2matrix2d(md_untilted_exp, angles_u);
	md2matrix2d(md_tilted_exp, angles_t);

	spheresplit(ang_acc, mdangles, number_of_regions);

	//Assigning directions to all particles (untilted and tilted)
	//Input  -  untiltparticles
	//		    tiltparticles
	//			Init Volume
	//Output -  metata with angles

	//Creating groups with same direction
	size_t idx, idx_reading = 0, idx_reading_last;
	FOR_ALL_OBJECTS_IN_METADATA2(md_untilted_exp, md_tilted_exp)
	{
		md_untilted_exp.getValue(MDL_ANGLE_TILT, tilt_u, __iter.objId);
		md_tilted_exp.getValue(MDL_ANGLE_TILT, tilt_t, __iter2.objId);

		for (size_t j = 0; j<mdangles.size(); j++)
		{
			mdangles.getValue(MDL_IDX, idx, j);
			mdangles.getValue(MDL_ANGLE_TILT, tilt_gal, j);
			mdangles.getValue(MDL_ANGLE_ROT, rot_gal, j);
			mdangles.getValue(MDL_IDX, idx, j);
			if ((tilt_u<tilt_gal+ang_acc) && (tilt_gal-ang_acc<tilt_u))
			{
				if((rot_u<rot_gal+ang_acc) && (rot_gal-ang_acc<rot_u))
					md_untilted_exp.getValue(MDL_IDX, idx, __iter.objId);
			}
			if ((tilt_t<tilt_gal+ang_acc) && (tilt_gal-ang_acc<tilt_t))
			{
				if((rot_t<rot_gal+ang_acc) && (rot_gal-ang_acc<rot_t))
					md_tilted_exp.getValue(MDL_IDX, idx, __iter2.objId);
			}
		}
	}
	// Merging both metadata
	mdu_aux = md_untilted_exp;
	md_all = md_tilted_exp;
	md_all.unionAll(mdu_aux);
	fnallassignSph = formatString("particles@%s/All_angular_assignments_in_sphere.xmd", fnDir.c_str());
	md_all.sort(md_all_sorted, MDL_IDX);
	md_all_sorted.write(fnallassignSph);

	//Defining classes by direction
	for (size_t kk = 1; kk<(number_of_regions+1) ; kk++)
	{
		//Storing image with same direction into a vector
		for (size_t kk_md = 1; kk_md<md_all_sorted.size() ; kk_md++)
		{
			md_all_sorted.getValue(MDL_IDX, idx_reading, kk_md);
			if (kk == idx_reading)
			{
				md_all_sorted.getValue(MDL_IMAGE, fnimg, kk_md);
				adressesfn.push_back(fnimg);
			}
		}

		MultidimArray<double> img1_marr, img2_marr;
		//Once all image in that region have been selected, a correlation is performed
		for (size_t kk_ang1 = 0; adressesfn.size(); kk_ang1++)
		{
			for (size_t kk_ang2 = kk_ang1+1; adressesfn.size(); kk_ang2++)
			{
				fnimg_corr1 = adressesfn[kk_ang1];
				fnimg_corr1 = adressesfn[kk_ang2];
				img1.read(fnimg_corr1);
				img1_marr = img1();
				img1_marr.setXmippOrigin();
				img2.read(fnimg_corr2);
				img2_marr = img2();
				img2_marr.setXmippOrigin();
				corr = alignImages(img1_marr, img2_marr, transformation_matrix, true, aux, aux2, aux3);


			}
		}
	}










}
