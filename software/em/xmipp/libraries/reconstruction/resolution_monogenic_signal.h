/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 					  jlvilas@cnb.csic.es
 * 			   Carlos Oscar S. Sorzano            coss@cnb.csic.es (2016)
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

#ifndef _PROG_MONOGENIC_SIGNAL_RES
#define _PROG_MONOGENIC_SIGNAL_RES

#include <iostream>
#include <data/xmipp_program.h>
#include <data/xmipp_image.h>
#include <data/metadata.h>
#include <data/xmipp_fft.h>
#include <data/xmipp_fftw.h>
#include <math.h>
#include <limits>
#include <complex>


/**@defgroup SSNR resolution_ssnr (Spectral Signal to Noise Ratio)
   @ingroup ReconsLibrary */
//@{
/** SSNR parameters. */

class ProgMonogenicSignalRes : public XmippProgram
{
public:
	 /** Filenames */
	FileName fnDir, fnVol, fnMask;

	/** Filter type */
	double fType, smpr;

	/** Is the volume previously masked?*/
	bool condpremask;

public:

    void defineParams();
    void readParams();
    void filterVolume();
    void RieszTransform3Dreal(const MultidimArray<double> &inputVol, std::vector<MultidimArray<double> > &RieszVector);
    void amplitudeMonogenicSignal3D(const MultidimArray<double> &inputVol,
			  	  	  	  	  	  const std::vector<MultidimArray<double> > &RieszVector,
			  	  	  	  	  	  	  	  	  	  	  MultidimArray<double> &amplitude);
    void passbandfiltervol(const MultidimArray<double> &inputVol, double freq, MultidimArray<double> &filteredVol, double sigma = 1);

    void passbandfiltervolFourier(const MultidimArray<double> &inputVol, MultidimArray< std::complex<double> > fftVol,
			double freq, MultidimArray< std::complex<double> > fftVfiltered_vol, double sigma);

    void highpassfiltervol(const MultidimArray<double> &inputVol, double freq, MultidimArray<double> &filteredVol, double raised_w);

    void FourierAmplitudeMonogenicSignal3D(const MultidimArray<double> &inputVol, const std::vector<MultidimArray<double> > &RieszVector,
    											MultidimArray<double> &amplitude, MultidimArray< std::complex<double> > &fftamplitude);

    void medianFilter3x3x3(const MultidimArray<double> &inputVol, MultidimArray<double> &FilteredVol);

    void applymask(const MultidimArray<double> &amplitudeMS, const MultidimArray<double> &mask, MultidimArray<double> &mask_amplitudeMS);

    void getHalfDimensions(const  MultidimArray<double> &Vol, size_t &xdim, size_t &ydim, size_t &zdim);

    double calculateMaskRadius(const MultidimArray<double> &Vol, size_t &xdim, size_t &ydim, size_t &zdim);
    void run();

};
//@}
#endif
