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

#ifndef VALIDATE_TILT_PAIRS_H
#define VALIDATE_TILT_PAIRS_H

#include <data/xmipp_fftw.h>
#include <data/args.h>
#include <data/xmipp_funcs.h>

#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <data/xmipp_image.h>
#include <data/geometry.h>
#include <data/filters.h>
#include <data/xmipp_program.h>
#include <complex>

/**@defgroup Centilt align_tilt_pairs (Align tilted and untilted images in a random conical tilt experiment)
   @ingroup ReconsLibrary */
//@{
class ProgValidationTiltPairs: public XmippProgram
{
public:
    /**  Filename untilt and tilt particles, reference volume */
    FileName fntilt, fnuntilt, fnOut, fnVol, fnuntilt2assign, fntilt2assign;

    /** Double for sampling projections, each projection will be performed each smprt degrees*/
    double smprt;

    /*** Maximum shift for aligning */
    int maxshift;

public:
    /// Define parameters in the command line
    void defineParams();

    /// Read parameters from the command line
    void readParams();

    /// Generating a galllery of projections
    void generateProjections(FileName fnVol, double smprt);

    void generateFourierStackTP(const MultidimArray<double> &input_stack, std::vector< AlignmentTransforms> &galleryTransforms_Test);

    void assignAngles(const MetaData mduntilt_exp, FileName fnun_out);

    /// Execute de program
    void run();


};
//@}
#endif
