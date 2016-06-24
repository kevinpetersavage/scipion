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

#ifndef PARTICLE_UNROCK3D_H
#define PARTICLE_UNROCK3D_H

#include <data/xmipp_program.h>

//#include "opencv2/imgcodecs.hpp"
//#include "opencv2/core/core.hpp"
//#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/imgproc/imgproc.hpp"
//#include "opencv2/video/video.hpp"
//#include "opencv2/core/utility.hpp"

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/video/video.hpp"
//#include "opencv2/gpu/gpu.hpp"

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
	double shiftLimit, scaleLimit, shearLimit, w1;
	double rotLimit;
	//std::vector< cv::Mat > transformationVec;
	std::vector<  Matrix1D<double> > transformationVec;
	std::vector<  cv::Mat > transformationVecCv;
	bool useECC;
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

    void extractAverageParticle();

    /// Extract a particle from an input frame
    void micExtractParticle(const int x, const int y,
    		const MultidimArray<double> frameImage,
			MultidimArray<double> &particleImage);

    /// To compute the average of each stack
    void computeAvgStack();

    /// To compute the radial average of each stack
    void computeRadialAvg();

    /// Compute the average of each group of the particles
    void computeGroupAvg(size_t start, size_t particleIndex,size_t end, MultidimArray<double> &groupAverage);

    /// To do the 3D polishing of particles
    void particlePolishing3D();

    /// Define Parameters
    void defineParams();

    /** Run */
    void run();
};

//@}
#endif
