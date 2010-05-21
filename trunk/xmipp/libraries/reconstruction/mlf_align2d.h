/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2004)
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

#ifndef MLFALIGN2D_H
#define MLFALIGN2D_H

#include <data/fftw.h>
#include <data/args.h>
#include <data/funcs.h>
#include <data/metadata.h>
#include <data/image.h>
#include <data/geometry.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/ctf.h>
#include <vector>
#include <numeric>

/**@defgroup MLFalign2D mlf_align2d (Maximum likelihood in 2D in Fourier space)
   @ingroup ReconsLibraryPrograms */
//@{

#define FOR_ALL_MODELS() for (int refno=0;refno<n_ref; refno++)
#define FOR_ALL_ROTATIONS() for (int ipsi=0; ipsi<nr_psi; ipsi++ )
#define FOR_ALL_FLIPS() for (int iflip=0; iflip<nr_flip; iflip++)
#define FOR_ALL_LIMITED_TRANSLATIONS() for (int itrans=0; itrans<nr_trans; itrans++)
#define FOR_ALL_DEFOCUS_GROUPS() for (int ifocus=0; ifocus<nr_focus; ifocus++)
#define SIGNIFICANT_WEIGHT_LOW 1e-8
#define SMALLVALUE 1e-4
#define SMALLANGLE 1.75
#define HISTMIN -6.
#define HISTMAX 6.
#define HISTSTEPS 120
#define DATALINELENGTH 11

#define FOR_ALL_GLOBAL_IMAGES() \
    for (int imgno = 0; imgno < nr_images_global; imgno++)
#define FOR_ALL_LOCAL_IMAGES() \
    for (int imgno = myFirstImg; imgno <= myLastImg; imgno++)
#define IMG_LOCAL_INDEX (imgno - myFirstImg)


/** MLFalign2D parameters. */
class Prog_MLFalign2D_prm
{
public:
    /** Filenames reference selfile/image, fraction docfile & output rootname */
    FileName fn_sel, fn_ref, fn_root, fn_frac, fn_sig, fn_doc, fn_ctfdat, fn_oext;
    /** Command line */
    std::string cline;
    /** sigma-value for origin offsets */
    double sigma_offset;
    /** Vector containing estimated fraction for each model */
    std::vector<double> alpha_k;
    /** Vector containing estimated fraction for mirror of each model */
    std::vector<double> mirror_fraction;
    /** Flag for checking mirror images of all references */
    bool do_mirror;
    /** Flag whether to fix estimates for model fractions */
    bool fix_fractions;
    /** Flag whether to fix estimate for sigma of origin offset */
    bool fix_sigma_offset;
    /** Flag whether to fix estimate for sigma of noise */
    bool fix_sigma_noise;
    /** Starting iteration */
    int istart;
    /** Number of iterations to be performed */
    int Niter;
    /** dimension of the images */
    int dim, dim2, hdim;
    /** Number of steps to sample in-plane rotation in 90 degrees */
    int nr_psi, max_nr_psi;
    /** Number of operations in "flip-array" (depending on do_mirror) */
    int nr_flip;
    /** Sampling rate for in-plane rotation */
    float psi_step;
    /** Total degrees in FOR_ALL_ROTATIONS */
    double psi_max;
    /** Vary psi and translational sampling with resolution */
    bool do_variable_psi, do_variable_trans;
    /** Total number of no-mirror rotations in FOR_ALL_FLIPS */
    int nr_nomirror_flips;
    /** Number of reference images */
    int n_ref;
    /** Total number of experimental images */
    int nr_images_global;
    /** Total number of local mpi images */
    int nr_images_local;
    /** First and last images, useful for mpi*/
    int myFirstImg, myLastImg;
    /** Verbose level:
        1: gives progress bar (=default)
        0: gives no output to screen at all */
    int verb;
    /** Stopping criterium */
    double eps;
    /** SelFile images (working, test and reference set) */
    MetaData MDimg, MDref;
    //Vector of image IDs in the MetaData object (change order for randomize)
    std::vector<long int> img_id;
    /** vector for flipping (i.e. 90/180-degree rotations) matrices */
    std::vector<Matrix2D<double> > F;
    /** Vector for images to hold references (new & old) */
    std::vector< Image<double> > Iref, Iold, Ictf;
    /** Matrices for calculating PDF of (in-plane) translations */
    MultidimArray<double> P_phi, Mr2;
    /** Fast mode */
    bool fast_mode;
    /** Fast mode */
    double C_fast;
    /** Limit translational searches */
    bool limit_trans;
    /** Number of limited translations */
    int nr_trans;
    /** Number for which limited translation is zero */
    int zero_trans;
    /** Offsets for limited translations */
    std::vector<Matrix1D<double> > Vtrans;
    /** Limited search range for origin offsets */
    int search_shift;
    /** Limit orientational searches */
    bool limit_rot;
    /** Limited search range for projection directions */
    double search_rot;
    /** Vectors to store old phi and theta for all images */
    std::vector<float> imgs_oldphi, imgs_oldtheta;
    /** Number of subdirectories to keep for unique offsets filenames */
    int offsets_keepdir;
    /** Flag for using ML3D */
    bool do_ML3D;
    /** Flag for generation of initial models from random subsets */
    bool do_generate_refs;
    /** Vector to store optimal origin offsets (if not written to disc) */
    std::vector<std::vector<double> > imgs_offsets;

    /** CTFDat file for all images */
    CTFDat ctfdat;
    /** Flag whether to include CTFs in the image formation model */
    bool do_ctf_correction;
    /** Pixel size in Angstroms */
    double sampling;
    /** Vector with number of images per defocuss group */
    std::vector<int> count_defocus;
    /** Flag whether the phases of the experimental images are flipped already */
    bool phase_flipped;
    /** Matrix with resolution shell at each Fourier pixel */
    MultidimArray<int> Mresol_int;
    /** Vectors with sigma2 (for each defocus group) */
    std::vector<Matrix1D<double> > Vsig, Vctf, Vdec;
    /** Multiplicative factor for SSNR */
    double reduce_snr;
    /** number of defocus groups */
    int nr_focus;
    /** Overall low and high resolution cutoffs for fourier-mode (in Fourier pixels) */
    double lowres_limit, highres_limit, ini_highres_limit;
    /** Do not multiply signal with CTF in the first iteration */
    bool first_iter_noctf;
    /** Divide by CTF (until first zero) instead of wiener filter */
    bool do_divide_ctf;
    /** Include all frequencies in the refinement */
    bool do_include_allfreqs;
    /** Fix high-resolution limit */
    double fix_high;
    /** Pointers to the 2D matrices (in FourierTransformHalf format) */
    std::vector<int> pointer_2d, pointer_i, pointer_j;
    int nr_points_prob, nr_points_2d, dnr_points_2d;
    /** Current highest resolution shell */
    int current_highres_limit;
    /** Random number generator seed */
    int seed;
    /** Matrix2D for mpi passing of docfiledata */
    Matrix2D<double> docfiledata;
    /** Flag for restart */
    bool do_restart;

    /// IN DEVELOPMENT

    /// USe t-distribution instead of normal one
    /** Use t-student distribution instead of normal one */
    bool do_student;
    /** Degrees of freedom for the t-student distribution */
    double df, df2;
    /** Perform sigma-trick for faster convergence (Mclachlan&Peel, p. 228)*/
    bool do_student_sigma_trick;

    /// Re-normalize internally
    /** Flag to refine normalization of each experimental image */
    bool do_norm;
    /** Grey-scale correction values */
    std::vector<double> imgs_scale, refs_avgscale;
    /** Overall average scale (to be forced to one)*/
    double average_scale;

    /// Statistical analysis of the noise distributions
    /** Perform Kolmogorov-Smirnov test on noise distribution */
    bool do_kstest;
    /** Iteration at which to write out histograms */
    int iter_write_histograms;
    /** Average histogram */
    Histogram1D sumhist;
    std::vector<Histogram1D > resolhist;

    /** debug flag */
    int debug;

public:

    /// Read arguments from command line
    void read(int argc, char **argv, bool ML3D = false);

    /// Show
    void show(bool ML3D = false);

    /// Usage
    void usage();

    /// Extended Usage
    void extendedUsage(bool ML3D = false);

    /// Setup lots of stuff
    void produceSideInfo(int rank = 0);

    /// Read reference images in memory & set offset vectors
    /// (This produce_side_info is Selfile-dependent!)
    void produceSideInfo2(int nr_vols = 1, int size = 1, int rank = 0);

    /// Randomize image order
    void randomizeImagesOrder();

    /// Calculate initial sigma2 from average power spectrum of the
    /// experimental images
    void estimateInitialNoiseSpectra();

    /// Calculate Wiener filter for defocus series as defined by Frank
    /// (2nd ed. formula 2.32b on p.60)
    void updateWienerFilters(Matrix1D<double> &spectral_signal,
                             std::vector<double> &sumw_defocus, int iter);

    /// Vary in-plane and translational sampling rates with resolution
    void setCurrentSamplingRates(double current_probres_limit);

    /// Generate initial references from random subset averages
    void generateInitialReferences();


    /// Calculate probability density distribution for in-plane transformations
    void calculateInPlanePDF();

    // Append a FourierTransform (in half format!) to a vector
    void appendFTtoVector(const Matrix2D<std::complex<double> > &Fin,
                          std::vector<double> &out);

    // get a FourierTransform (in half format!) from a vector
    void getFTfromVector(const std::vector<double> &in,
                         const int start_point,
                         Matrix2D<std::complex<double> > &out,
                         bool only_real = false);

    /// Fill vector of matrices with all rotations of reference
    void rotateReference(std::vector<double> &out);

    /// Apply reverse rotations to all matrices in vector and fill new matrix with their sum
    void reverseRotateReference(const std::vector<double> &in,
                                std::vector<Matrix2D<double > > &out);

    /// Calculate which references have projection directions close to
    /// phi and theta
    void preselectDirections(float &phi, float &theta,
                             std::vector<double> &pdf_directions);

    /// Pre-calculate which model and phi have significant probabilities
    /// without taking translations into account!
    void preselect_significant_model_phi(Matrix2D<double> &Mimg, std::vector<double> &offsets,
                                         std::vector <std::vector< Matrix2D<double > > > &Mref,
                                         Matrix2D<int> &Msignificant,
                                         std::vector<double > &pdf_directions);

    // Calculate the FT of a translated matrix using a phase shift in
    // Fourier space
    void fourierTranslate2D(const std::vector<double> &in,
                            Matrix1D<double> &trans,
                            std::vector<double> &out,
                            int point_start = 0);

    // If not determined yet: search optimal offsets using maxCC
    // Then for all optimal translations, calculate all translated FTs
    // for each of the flipped variants
    void calculateFourierOffsets(const Matrix2D<double> &Mimg,
                                 const std::vector<double > &offsets,
                                 std::vector<double>  &out,
                                 Matrix2D<int> &Moffsets,
                                 Matrix2D<int> &Moffsets_mirror);

    /// Perform expectation step for a single image
    void processOneImage(const Matrix2D<double> &Mimg,
                         const int focus, bool apply_ctf,
                         const std::vector<double> &Fref,
                         std::vector<double> &Fwsum_imgs,
                         std::vector<double> &Fwsum_ctfimgs,
                         std::vector<double> &Mwsum_sigma2,
                         double &wsum_sigma_offset,
                         std::vector<double> &sumw,
                         std::vector<double> &sumw2,
                         std::vector<double> &sumwsc,
                         std::vector<double> &sumwsc2,
                         std::vector<double> &sumw_mirror,
                         double &LL, double &fracweight,  double &maxweight2, double &sum_refw2,
                         double &opt_scale, int &opt_refno, double &opt_psi,
                         int &opt_ipsi, int &opt_iflip,
                         Matrix1D<double> &opt_offsets,
                         std::vector<double> &opt_offsets_ref,
                         std::vector<double > &pdf_directions,
                         bool do_kstest, bool write_histograms,
                         FileName fn_img, double &KSprob);

    /// Perform Kolmogorov-Smirnov test
    double performKSTest(Matrix2D<double> &Mimg,  const int focus, bool apply_ctf,
                         FileName &fn_img, bool write_histogram,
                         std::vector <std::vector< Matrix2D<std::complex<double> > > > &Fref,
                         double &opt_scale, int &opt_refno, int &opt_ipsi, int &opt_iflip,
                         Matrix1D<double> &opt_offsets);

    /// Integrate over all experimental images
    void expectation(std::vector< ImageXmippT<double> > &Iref, int iter,
                     double &LL, double &sumcorr,
                     std::vector<Matrix2D<double> > &wsum_Mref,
                     std::vector<Matrix2D<double> > &wsum_ctfMref,
                     std::vector<std::vector<double> > &Mwsum_sigma2,
                     double &wsum_sigma_offset,
                     std::vector<double> &sumw,
                     std::vector<double> &sumw2,
                     std::vector<double> &sumwsc,
                     std::vector<double> &sumwsc2,
                     std::vector<double> &sumw_mirror,
                     std::vector<double> &sumw_defocus);

    /// Update all model parameters (maximization step)
    void maximization(std::vector<Matrix2D<double> > &wsum_Mref,
                      std::vector<Matrix2D<double> > &wsum_ctfMref,
                      std::vector<std::vector<double> > &Mwsum_sigma2,
                      double &wsum_sigma_offset,
                      std::vector<double> &sumw,
                      std::vector<double> &sumw2,
                      std::vector<double> &sumwsc,
                      std::vector<double> &sumwsc2,
                      std::vector<double> &sumw_mirror,
                      std::vector<double> &sumw_defocus,
                      double &sumcorr, double &sumw_allrefs,
                      Matrix1D<double> &spectral_signal,
                      int refs_per_class=1);

    /// check convergence
    bool checkConvergence(std::vector<double> &conv);

    /// Write out reference images, selfile and logfile
    void writeOutputFiles(const int iter,
                          double &sumw_allrefs, double &LL, double &avecorr,
                          std::vector<double> &conv);

    /// Write partial docfile
    void addPartialDocfileData(Matrix2D<double> data, int first, int last);

};
//@}

#endif
