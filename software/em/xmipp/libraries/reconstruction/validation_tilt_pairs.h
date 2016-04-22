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
    FileName fntilt, fnuntilt, fnOut, fnSym;

public:
    /// Define parameters in the command line
    void defineParams();

    /// Read parameters from the command line
    void readParams();

    void validate(MetaData &md_u, MetaData &md_t, MetaData &md_out, MetaData &md_validation);

    void powerIterationMethod(const Matrix2D<double> A, Matrix1D<double> &eigenvector, double &eigenvalueapprox);

    double R0_ThresholdFisher(int N);

    void R_Fisher(MetaData md, MetaData &md_scattering);

    double R_Value(MetaData md, double &SumRx, double &SumRy, double &SumRz);

    Matrix2D<double> RodriguesMatrix( const double theta, const Matrix1D<double> axis);

    double confidenceSolidAngle(double R, int Nparticles, double confidence = 0.01);

    /// Execute the program
    void run();


};
//@}
#endif
