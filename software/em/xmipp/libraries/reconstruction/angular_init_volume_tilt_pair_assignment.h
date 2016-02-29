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

#ifndef INIT_VOLUME_TILT_PAIR_ASSIGNMENT_H_
#define INIT_VOLUME_TILT_PAIR_ASSIGNMENT_H_
#define PI 3.14159265



#include <data/xmipp_program.h>
#include <data/filters.h>
#include <math.h>
#include <alglib/src/ap.h>
//#include <fourier_filter.h>

class ProgInitVolumeTiltPairassignment: public XmippProgram
{


public:
    /** Filenames */
    FileName fnUntilt, fnTilt, fnDir, fnVol, fnSym, fnmic;

    /** Particle size, sampling rate*/
    double smprt, alphaU, alphaT, tilt_mic;

    /** Maximum shif for aligning particles*/
    int maxshift, fnOutVol, pad;


public:

    void defineParams();

    void readParams();

    void run();

    void generateProjections(FileName &fnVol, double &smprt);

    void generateFourierStack(const MultidimArray<double> &input_stack,	std::vector< AlignmentTransforms> &galleryTransforms_Test);

    void generateInitialBall(const MetaData &md_u,const MetaData &md_t, MetaData &md_u_assign_iter0, MetaData &md_t_assign_iter0, FileName &fnVol);

    void adaptativeMask(const FileName fnVol, const MetaData &md_angles, MetaData &md_out_filtered, double elongatefactor = 0.15);


};
#endif
