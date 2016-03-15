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

#ifndef ANGULAR_TILT_PAIR_ASSIGNMENT_H_
#define ANGULAR_TILT_PAIR_ASSIGNMENT_H_
#define PI 3.14159265



#include <data/xmipp_program.h>
#include <data/filters.h>
#include <math.h>
#include <alglib/src/ap.h>
//#include <fourier_filter.h>

class ProgAngularTiltPairHeterogeneity: public XmippProgram
{


public:
    /** Filenames */
    FileName fnUntilt, fnTilt, fnDir, fnVol;

    /** sampling rate*/
    double ang_dist, ang_acc;

    /** number of classes*/
    int N_cls;



public:

    void defineParams();

    void readParams();

    void run();

    void md2matrix2d(const MetaData md, Matrix2D<double> &output_matrix);

    void md2matrix2dclasses(const MetaData md, const MetaData mdangles, Matrix2D<double> &output_matrix);

    void spheresplit(double and_acc, MetaData &mdangles, int &number_of_regions);


};
#endif
