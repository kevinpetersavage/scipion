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

#include "image_tilt_pair_heterogeneity.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <data/xmipp_image.h>
#include <geometry.h>
#include <data/matrix1d.h>
#include <data/matrix2d.h>
#define DEBUG


void ProgImageTiltPairHeterogeneity::readParams()
{
	fnUntilt = getParam("--untiltparticles");
	fnTilt = getParam("--tiltparticles");
	fnDir = getParam("--odir");
	ang_dist = getDoubleParam("--angular_distance");
	ang_acc = getDoubleParam("--angular_accuracy");
	N_cls = getIntParam("--classes");
}


void ProgImageTiltPairHeterogeneity::defineParams()
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


void ProgImageTiltPairHeterogeneity::md2matrix2d(const MetaData md, Matrix2D<double> &output_matrix)
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


void ProgImageTiltPairHeterogeneity::md2matrix2dclasses(const MetaData md, const MetaData mdangles, Matrix2D<double> &output_matrix)
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


void ProgImageTiltPairHeterogeneity::run()
{
	std::cout << "Starting..." << std::endl;

	MetaData md_untilted_exp, md_tilted_exp, mdangles, md_all_particles_aux, md_all_particles, md_all_particles_set1, md_all_particles_set2;
	double rot_u, tilt_u, psi_u, rot_t, tilt_t, psi_t, rot_gal, tilt_gal, threshold_split;
	Matrix2D<double> angles_u, angles_t;
	size_t len_u, len_t, len_all, idx, idx_0;
	std::vector<size_t> randomvector;
	MDRow mdall_row;
	FileName fn_md_all_particles_set1, fn_md_all_particles_set2;

	//Reading MD and storing them into a matrix
	md_untilted_exp.read(fnUntilt);
	md_tilted_exp.read(fnTilt);

	len_u = md_untilted_exp.size();
	len_t = md_tilted_exp.size();

	md_all_particles_aux = md_untilted_exp;
	md_all_particles = md_tilted_exp;

	md_all_particles_aux.operate(formatString("particleId = particleId + %d", len_u));
	md_all_particles.unionAll(md_all_particles_aux);

	len_all = md_all_particles.size();

	//Splitting the initial set of particle in two sets
	randomize_random_generator();

	for (size_t i=0; i<len_all; ++i) randomvector.push_back(i+1); // 1 2 3 4 5 6 7 8 9

	std::random_shuffle ( randomvector.begin(), randomvector.end());

	//Instead of performing a random assignment, split in three orthogonal directions and create 2^3 volumes

	std::cout << "len_all = " << len_all << std::endl;
	threshold_split = floor(len_all/2);
	for (size_t i=0; i<len_all; ++i)
	{
		if (i < threshold_split)
		{
			md_all_particles.getRow(mdall_row, randomvector[i]);
			idx = md_all_particles_set1.addObject();
			md_all_particles_set1.setRow(mdall_row, idx);
		}
		else
		{
			md_all_particles.getRow(mdall_row, randomvector[i]);
			idx = md_all_particles_set2.addObject();
			md_all_particles_set2.setRow(mdall_row, idx);
		}
	}

	fn_md_all_particles_set1 = formatString("particles@%s/inputparticles_set1.xmd", fnDir.c_str());
	fn_md_all_particles_set2 = formatString("particles@%s/inputparticles_set2.xmd", fnDir.c_str());

	md_all_particles_set1.write(fn_md_all_particles_set1);
	md_all_particles_set2.write(fn_md_all_particles_set2);
}



