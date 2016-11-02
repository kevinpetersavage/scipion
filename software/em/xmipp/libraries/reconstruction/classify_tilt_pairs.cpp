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

#include "classify_tilt_pairs.h"
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


void ProgClassifyTiltPairs::readParams()
{
	fnMdUntilted1 = getParam("--md_Untilted1");
	fnMdTilted1 = getParam("--md_Tilted1");
	fnMdUntilted2 = getParam("--md_Untilted2");
	fnMdTilted2 = getParam("--md_Tilted2");
	fnUntiltedOdir1 = getParam("--odirMdUntilted_Vol1");
	fnTiltedOdir1 = getParam("--odirMdTilted_Vol1");
	fnUntiltedOdir2 = getParam("--odirMdUntilted_Vol2");
	fnTiltedOdir2 = getParam("--odirMdTilted_Vol2");
}
ProgClassifyTiltPairs::ProgClassifyTiltPairs()
{
	rank=0;
	Nprocessors=1;
}

void ProgClassifyTiltPairs::defineParams()
{
	//usage
	addUsageLine("From an untilted particles stack and an initial volume, an angular assignment of the untilted particles is performed."
			"Then, knowing the tilt axis and the tilt angle of the micrography, tiltmic, and the rotational angles of the untilted and"
			"tilted micrographs, alphaU and alphaT respectively, an angular assignment of the tilted particles is determined");
	//params
	addParamsLine("  [--md_Untilted1 <md_file=\"\">]    : Untilt particles stack");
	addParamsLine("  [--md_Tilted1 <md_file=\"\">]    : Tilt particles stack");
	addParamsLine("  [--md_Untilted2 <md_file=\"\">]    : Tilt particles stack");
	addParamsLine("  [--md_Tilted2 <md_file=\"\">]    : Tilt particles stack");
	addParamsLine("  [--odirMdUntilted_Vol1 <outputDir=\".\">]      : Output directory with a metadata file which contains the untilted particles assigned to volume 1");
	addParamsLine("  [--odirMdTilted_Vol1 <outputDir=\".\">]      : Output directory with a metadata file which contains the tilted particles assigned to volume 2");
	addParamsLine("  [--odirMdUntilted_Vol2 <outputDir=\".\">]      : Output directory with a metadata file which contains the untilted particles assigned to volume 1");
	addParamsLine("  [--odirMdTilted_Vol2 <outputDir=\".\">]      : Output directory with a metadata file which contains the tilted particles assigned to volume 2");
}


void ProgClassifyTiltPairs::run()
{
	MetaData mdUntilted1, mdTilted1, mdUntilted2, mdTilted2, mdUntiltedout1, mdTiltedout1, mdUntiltedout2, mdTiltedout2;
	MDRow mdUntilted1_row, mdTilted1_row, mdUntilted2_row, mdTilted2_row;
	size_t len_mdUntilted1, len_mdTilted1, len_mdUntilted2, len_mdTilted2, len;

	mdUntilted1.read(fnMdUntilted1);
	mdTilted1.read(fnMdTilted1);
	mdUntilted2.read(fnMdUntilted2);
	mdTilted2.read(fnMdTilted2);

	len_mdUntilted1 = mdUntilted1.size();
	len_mdTilted1 = mdTilted1.size();
	len_mdUntilted2 = mdUntilted2.size();
	len_mdTilted2 = mdTilted2.size();

	if (len_mdUntilted1>len_mdTilted1)
	{
		len = len_mdUntilted1;
	}
	else
	{
		len = len_mdTilted1;
	}

	size_t particle_1, particle_2, particle_3, particle_4, objId=1, idx=1;

	double corr1, corr2, corr_vol1, corr_vol2;
	MetaData mdUntilted1_sorted, mdUntilted2_sorted, mdTilted1_sorted, mdTilted2_sorted;

	mdUntilted1_sorted.sort(mdUntilted1, MDL_PARTICLE_ID);
	mdUntilted2_sorted.sort(mdUntilted2, MDL_PARTICLE_ID);
	mdTilted1_sorted.sort(mdTilted1, MDL_PARTICLE_ID);
	mdTilted2_sorted.sort(mdTilted2, MDL_PARTICLE_ID);

	for (size_t k = 0; k<len+1; k++)
	{
		mdUntilted1_sorted.getRow(mdUntilted1_row, objId);
		mdTilted1_sorted.getRow(mdTilted1_row, objId);
		mdUntilted2_sorted.getRow(mdUntilted2_row, objId);
		mdTilted2_sorted.getRow(mdUntilted2_row, objId);

		mdUntilted1_sorted.getValue(MDL_PARTICLE_ID, particle_1, objId);
		mdTilted1_sorted.getValue(MDL_PARTICLE_ID, particle_2, objId);
		mdUntilted2_sorted.getValue(MDL_PARTICLE_ID, particle_3, objId);
		mdTilted2_sorted.getValue(MDL_PARTICLE_ID, particle_4, objId);

		if ((particle_1 == particle_2) && (particle_3 == particle_4) && (particle_1 == particle_3))
		{
			mdUntilted1_sorted.getValue(MDL_MAXCC, corr1, objId);
			mdTilted1_sorted.getValue(MDL_MAXCC, corr2, objId);
			corr_vol1 = 0.5*(corr1 + corr2);

			mdUntilted2_sorted.getValue(MDL_MAXCC, corr1, objId);
			mdTilted2_sorted.getValue(MDL_MAXCC, corr2, objId);
			corr_vol2 = 0.5*(corr1 + corr2);

			if (corr_vol1>corr_vol2)
			{
				mdUntiltedout1.addRow(mdUntilted1_row);
				mdTiltedout1.addRow(mdTilted1_row);
			}
			else
			{
				//mdout2.setRow(md2_row, idx);
				mdUntiltedout2.addRow(mdUntilted2_row);
				mdTiltedout2.addRow(mdTilted2_row);
			}

		}
		//std::cout << "iter %f" << k << std::endl;
		objId++;
		idx++;
	}

	mdUntiltedout1.write(fnUntiltedOdir1);
	mdTiltedout1.write(fnTiltedOdir1);
	mdUntiltedout2.write(fnUntiltedOdir2);
	mdTiltedout2.write(fnTiltedOdir2);
}
