/***************************************************************************
 * Authors:     AUTHOR_NAME (jlvilas@cnb.csic.es)
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

#ifndef CLASSIFY_TILT_PAIRS_H_
#define CLASSIFY_TILT_PAIRS_H_
#define PI 3.14159265



#include <data/xmipp_program.h>
#include <data/filters.h>
#include <math.h>
#include <alglib/src/ap.h>
//#include <fourier_filter.h>

class ProgClassifyTiltPairs: public XmippProgram
{
public:
	// Internal members
	size_t rank, Nprocessors;

	// Auxialiar metadata for storing data
	MetaData mdPartial_u, mdPartial_t;

public:
    /** Filenames */
    FileName fnMdUntilted1, fnMdTilted1, fnMdUntilted2, fnMdTilted2, fnUntiltedOdir1, fnTiltedOdir1,
    fnUntiltedOdir2, fnTiltedOdir2, fnVol, fnVol2, fnSym, fnmic;

    /** Particle size, sampling rate*/
    double smprt, alphaU, alphaT, tilt_mic;

    /** Maximum shif for aligning particles*/
    int maxshift, fnOutVol, pad;

public:
    /// Empty constructor
    ProgClassifyTiltPairs();

public:

    void defineParams();

    void readParams();

    void run();

    void generateProjections(const FileName &fnVol, double &smprt, FileName galleryfn);

    void generateFourierStack(const MultidimArray<double> &input_stack,	std::vector< AlignmentTransforms> &galleryTransforms_Test);

//    /// Gather alignment
//	virtual void gatherAlignment() {}

	/// Synchronize with other processors
	virtual void synchronize() {}

};
#endif
