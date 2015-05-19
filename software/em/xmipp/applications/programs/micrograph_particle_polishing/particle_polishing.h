/***************************************************************************
*
* Authors:   Vahid Abrishami         vabrishami@cnb.csic (2015)
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

#ifndef PARTICLE_POLISHING_H
#define PARTICLE_POLISHING_H

#include <data/xmipp_program.h>

class ProgParticlePolishing: public XmippProgram
{
public:
    FileName fnPolishedParticles;
    /// Movie filename
    FileName fnMovie, fnParticleStack;
    /// Shifts from Optical Flow
    FileName fnAvg;
    // File name for the initial estimated optical flows
    FileName fnInitFlow;
    /// Particle coordinates from corrected average (Optical Flow)
    FileName fnparticleCoords;

    /// For reading the input movie and particle coordinates
    ImageGeneric movieStack;
    MetaData particleCoords;

    /// Structure for computing the affine transform
    MultidimArray<double> currentParticle, avgParticle, extParticles;
    /// To keep the average of instances of a particle in different frames
    MultidimArray<double> particleAvgStack;
    Matrix2D<double> A;

    size_t particleRadius, particleSize;
    size_t frameNum, particleNum, iterationNum;
    double shiftLimit, scaleLimit, shearLimit;
    std::vector< Matrix1D<double> > transformationVec;

public:

    ProgParticlePolishing();

    ~ProgParticlePolishing();

    /// Read parameters
    void readParams();

    /// Show parameters
    void show();

    void produceSideInfo();

    /// Remove the shifts from the particles and put them in a stack
    void extractParticels();

    // Extract a noise image from a specific frame (frameImage).
    void random_noise_areas(const MultidimArray<double> frameImage, MultidimArray<double> &noiseStack);

    //void computing_SNR();

    void extractAverageParticle();

    /// Extract a particle from an input frame
    void micExtractParticle(const int x, const int y,
                            const MultidimArray<double> frameImage,
                            MultidimArray<double> &particleImage);

    /// To compute the average of each stack
    void computeAvgStack();

    /// To compute the radial average of each stack
    void computeRadialAvg();

    /// To do the 3D polishing of particles
    void particlePolishing3D();

    /// Define Parameters
    void defineParams();

    /// Compute de signal noise ratio
    void computing_SNR();

    void write_SNR(FileName fnfile, MultidimArray<double> SNR);

    /** Run */
    void run();
};

//@}
#endif
