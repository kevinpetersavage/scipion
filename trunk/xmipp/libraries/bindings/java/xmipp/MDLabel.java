package xmipp; 
public class MDLabel {
   public static final int MDL_UNDEFINED = -1;
   public static final int MDL_FIRST_LABEL = 0;  ///< The label MDL_OBJID is special and should not be used
   public static final int MDL_OBJID = MDL_FIRST_LABEL; ///< object id (int); NOTE: This label is special and shouldn't be used
   public static final int MDL_ANGLE_COMPARISON = 1;  ///< Angular comparison (see angular_distance.cpp)
   public static final int MDL_ANGLEPSI2 = 2;  ///< Psi angle of an image (double = 2; degrees)
   public static final int MDL_ANGLEPSI = 3;  ///< Psi angle of an image (double = 3; degrees)
   public static final int MDL_ANGLEROT2 = 4;  ///< Rotation angle of an image (double = 4; degrees)
   public static final int MDL_ANGLEROT = 5;  ///< Rotation angle of an image (double = 5; degrees)
   public static final int MDL_ANGLETILT2 = 6;  ///< Tilting angle of an image (double = 6; degrees)
   public static final int MDL_ANGLETILT = 7;  ///< Tilting angle of an image (double = 7; degrees)
   public static final int MDL_ASSOCIATED_IMAGE1 = 8;  ///< Image associated to this object (std::string)
   public static final int MDL_ASSOCIATED_IMAGE2 = 9;  ///< Image associated to this object (std::string)
   public static final int MDL_ASSOCIATED_IMAGE3 = 10;  ///< Image associated to this object (std::string)
   public static final int MDL_ASSOCIATED_IMAGE4 = 11;  ///< Image associated to this object (std::string)
   public static final int MDL_ASSOCIATED_IMAGE5 = 12;  ///< Image associated to this object (std::string)
   public static final int MDL_AVG = 13;  ///< average value (double)
   public static final int MDL_AZIMUTALANGLE = 14;  ///< ctf definition azimutal angle
   public static final int MDL_BGMEAN = 15;  ///< Mean background value for an image
   public static final int MDL_BLOCK = 16;  ///< Current block number (for incremental EM)
   public static final int MDL_CELLX = 17;  ///< Cell location for crystals
   public static final int MDL_CELLY = 18;  ///< Cell location for crystals
   public static final int MDL_COMMENT = 19;  ///< A comment for this object /*** NOTE THIS IS A SPECIAL CASE AND SO IS TREATED ***/
   public static final int MDL_COST = 20;  ///< Cost for the image (double)
   public static final int MDL_COUNT = 21;  ///< Number of elements of a type (int) [this is a genereic type do not use to transfer information to another program]
   public static final int MDL_CTFINPUTPARAMS = 22;  ///< Parameters file for the CTF Model (std::string)
   public static final int MDL_CTFMODEL = 23;  ///< Name for the CTF Model (std::string)
   public static final int MDL_CTFMODEL2 = 24;  ///< Name for another CTF model (std::string)
   public static final int MDL_CTF_SAMPLING_RATE = 25;  ///< Sampling rate
   public static final int MDL_CTF_SAMPLING_RATE_Z = 26;  ///< Sampling rate in Z direction
   public static final int MDL_CTF_VOLTAGE = 27;  ///< Microscope voltage (kV)
   public static final int MDL_CTF_DEFOCUSA = 28;  ///< aver (Angage defocusstroms)
   public static final int MDL_CTF_DEFOCUSU = 29;  ///< Defocus U (Angstroms)
   public static final int MDL_CTF_DEFOCUSV = 30;  ///< Defocus V (Angstroms)
   public static final int MDL_CTF_DEFOCUS_ANGLE = 31;  ///< Defocus angle (degrees)
   public static final int MDL_CTF_CS = 32;  ///< Spherical aberration
   public static final int MDL_CTF_CA = 33;  ///< Chromatic aberration
   public static final int MDL_CTF_ENERGY_LOSS = 34;  ///< Energy loss
   public static final int MDL_CTF_LENS_STABILITY = 35;  ///< Lens stability
   public static final int MDL_CTF_CONVERGENCE_CONE = 36;  ///< Convergence cone
   public static final int MDL_CTF_LONGITUDINAL_DISPLACEMENT = 37;  ///< Longitudinal displacement
   public static final int MDL_CTF_TRANSVERSAL_DISPLACEMENT = 38;  ///< Transversal displacemente
   public static final int MDL_CTF_Q0 = 39;  ///< Inelastic absorption
   public static final int MDL_CTF_K = 40;  ///< CTF gain
   public static final int MDL_CTFBG_GAUSSIAN_K = 41;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN_SIGMAU = 42;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN_SIGMAV = 43;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN_CU = 44;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN_CV = 45;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN_ANGLE = 46;  ///< CTF Background parameter
   public static final int MDL_CTFBG_SQRT_K = 47;  ///< CTF Background parameter
   public static final int MDL_CTFBG_SQRT_U = 48;  ///< CTF Background parameter
   public static final int MDL_CTFBG_SQRT_V = 49;  ///< CTF Background parameter
   public static final int MDL_CTFBG_SQRT_ANGLE = 50;  ///< CTF Background parameter
   public static final int MDL_CTFBG_BASELINE = 51;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN2_K = 52;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN2_SIGMAU = 53;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN2_SIGMAV = 54;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN2_CU = 55;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN2_CV = 56;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN2_ANGLE = 57;  ///< CTF Background parameter
   public static final int MDL_CTF_CRITERION_PSDCORRELATION90 = 58;  ///< PSD correlation at 90 degrees
   public static final int MDL_CTF_CRITERION_FIRSTZERORATIO = 59;  ///< First zero ratio
   public static final int MDL_CTF_CRITERION_FIRSTZEROAVG = 60;  ///< First zero average (in Angstroms)
   public static final int MDL_CTF_CRITERION_FIRSTZERODISAGREEMENT = 61;  ///< First zero disagreement with second model (in Angstroms)
   public static final int MDL_CTF_CRITERION_DAMPING = 62;  ///< Minimum damping at border
   public static final int MDL_CTF_CRITERION_PSDRADIALINTEGRAL = 63;  ///< Integral of the radial PSD
   public static final int MDL_CTF_CRITERION_FITTINGSCORE = 64;  ///< Score of the fitting
   public static final int MDL_CTF_CRITERION_FITTINGCORR13 = 65;  ///< Correlation between the 1st and 3rd ring of the CTF
   public static final int MDL_CTF_CRITERION_PSDVARIANCE = 66;  ///< PSD variance
   public static final int MDL_CTF_CRITERION_PSDPCA1VARIANCE = 67;  ///< Variance in the first principal component of the PSDs
   public static final int MDL_CTF_CRITERION_PSDPCARUNSTEST = 68;  ///< Runs test on the projection of the PSD on the first principal component
   public static final int MDL_CTF_CRITERION_COMBINED = 69;  ///< Combined criterion formed by several other criteria
   public static final int MDL_CTF_CRITERION_NORMALITY = 70;  ///< Normality test between histogram of micrography and gaussian distribution
   public static final int MDL_CTF_XRAY_DIMENSIONS = 71;  // Size in pixels of the 3D PSF to be created (Xdim = 71;  Ydim = 71;  Zdim)
   public static final int MDL_CTF_XRAY_LAMBDA = 72;  /// X-ray wavelength (nm)
   public static final int MDL_CTF_XRAY_LENS_TYPE = 73;  ///Algorithm used to generate Xray PSF
   public static final int MDL_CTF_XRAY_MAGNIFICATION = 74;  /// Magnification of the X-ray microscope
   public static final int MDL_CTF_XRAY_OUTER_ZONE_WIDTH = 75;  /// Outermost zone width of the X-ray Fresnel lens (nm)
   public static final int MDL_CTF_XRAY_ZONES_NUMBER = 76;  // Number of zones of the X-ray Fresnel lens
   public static final int MDL_DATATYPE = 77;  ///< if read from file original image datatype = 77;  this is an struct defined in image
   public static final int MDL_DEFGROUP = 78;  ///< Defocus group
   public static final int MDL_DM3_IDTAG = 79; 
   public static final int MDL_DM3_NODEID = 80; 
   public static final int MDL_DM3_NUMBER_TYPE = 81; 
   public static final int MDL_DM3_PARENTID = 82; 
   public static final int MDL_DM3_TAGCLASS = 83; 
   public static final int MDL_DM3_TAGNAME = 84; 
   public static final int MDL_DM3_SIZE = 85; 
   public static final int MDL_DM3_VALUE = 86; 
   public static final int MDL_ENABLED = 87;  ///< Is this image enabled? (int [-1 or 1])
   public static final int MDL_FLIP = 88;  ///< Flip the image? (bool)
   public static final int MDL_IMAGE_CLASS_COUNT = 89;  ///< Number of images assigned to the same class as this image
   public static final int MDL_IMAGE_CLASS_GROUP = 90;  ///< Name of the class group for this image (metadata with all the images assigned to that class)
   public static final int MDL_IMAGE_CLASS = 91;  ///< Name of the class representative for this image
   public static final int MDL_IMAGE = 92;  ///< Name of an image (std::string)
   public static final int MDL_IMAGE_ORIGINAL = 93;  ///< Name of an image from which MDL_IMAGE is coming from
   public static final int MDL_IMAGE_TILTED = 94;  ///< Name of the tilted images associated to MDL_IMAGE
   public static final int MDL_IMGMD = 95;  ///< Name of Metadata file for all images (string)
   public static final int MDL_INTSCALE = 96;  ///< Intensity scale for an image
   public static final int MDL_ITER = 97;  ///< Current iteration number (int)
   public static final int MDL_K = 98;  ///< //ctf definition K
   public static final int MDL_KSTEST = 99;  ///<KS-test statistics
   public static final int MDL_LL = 100;  ///< contribution of an image to log-likelihood value
   public static final int MDL_MASK = 101;  ///< Name of a mask associated to image
   public static final int MDL_MAXCC = 102;  ///< Cross-correlation for the image (double)
   public static final int MDL_MAX = 103;  ///<maximum value (double)
   public static final int MDL_MICROGRAPH = 104;  ///< Name of a micrograph (std::string)
   public static final int MDL_MICROGRAPH_TILTED = 105;  ///< Name of the corresponding tilted micrograph (std::string)
   public static final int MDL_MIN = 106;  ///<minimum value (double)
   public static final int MDL_MIRRORFRAC = 107;  ///< Mirror fraction for a Maximum Likelihood model
   public static final int MDL_MISSINGREGION_NR = 108;  ///< Number of missing region in subtomogram
   public static final int MDL_MISSINGREGION_TYPE = 109;  ///< Type of missing region in subtomogram
   public static final int MDL_MISSINGREGION_THY0 = 110;  ///< Initial tilt angle in Y for missing region in subtomogram
   public static final int MDL_MISSINGREGION_THYF = 111;  ///< Final tilt angle in Y for missing region in subtomogram
   public static final int MDL_MISSINGREGION_THX0 = 112;  ///< Initial tilt angle in X for missing region in subtomogram
   public static final int MDL_MISSINGREGION_THXF = 113;  ///< Final tilt angle in X for missing region in subtomogram
   public static final int MDL_MODELFRAC = 114;  ///< Model fraction (alpha_k) for a Maximum Likelihood model
   public static final int MDL_NMA = 115;  ///< Normal mode displacements (vector double)
   public static final int MDL_NOISE_ANGLES = 116;  ///< Noise description for projected angles
   public static final int MDL_NOISE_PARTICLE_COORD = 117;  ///< Noise description for particle's center coordenates (when projecting)
   public static final int MDL_NOISE_PIXEL_LEVEL = 118;  ///< Noise description for pixels' gray level (when projecting)
   public static final int MDL_ORDER = 119;  /// auxiliary label to be used as an index
   public static final int MDL_ORIGINX = 120;  ///< Origin for the image in the X axis (double)
   public static final int MDL_ORIGINY = 121;  ///< Origin for the image in the Y axis (double)
   public static final int MDL_ORIGINZ = 122;  ///< Origin for the image in the Z axis (double)
   public static final int MDL_PMAX = 123;  ///< Maximum value of normalized probability function (now called "Pmax/sumP") (double)
   public static final int MDL_PRJ_DIMENSIONS = 124;  // X = 124; Y dimensions for the generated projections
   public static final int MDL_PRJ_TILT_RANGE = 125;  // Vector with the initial and final tilt angle values = 125;  and step size
   public static final int MDL_PRJ_VOL = 126;         // Volume file name to generate projections from
   public static final int MDL_PSD = 127;  ///< A Power Spectrum Density file name (std::string)
   public static final int MDL_RANDOMSEED = 128;  ///< Seed for random number generator
   public static final int MDL_REF3D = 129;  ///< 3D Class to which the image belongs (int)
   public static final int MDL_REF = 130;  ///< Class to which the image belongs (int)
   public static final int MDL_REFMD = 131;  ///< Name of Metadata file for all references(string)
   public static final int MDL_RESOLUTION_DPR = 132;  ///<differential phase residual (double)
   public static final int MDL_RESOLUTION_ERRORL2 = 133;  ///<Error in l2 (double)
   public static final int MDL_RESOLUTION_FRC = 134;  ///<Fourier shell correlation (double)
   public static final int MDL_RESOLUTION_FRCRANDOMNOISE = 135;  ///<Fourier shell correlation noise (double)
   public static final int MDL_RESOLUTION_FREQ = 136;  ///<Frequency in 1/A (double)
   public static final int MDL_RESOLUTION_FREQREAL = 137;  ///< Frequency in A (double)
   public static final int MDL_SAMPLINGRATE = 138;  ///< sampling rate in A/pixel (double)
   public static final int MDL_SAMPLINGRATEX = 139;  ///< sampling rate in A/pixel (double)
   public static final int MDL_SAMPLINGRATEY = 140;  ///< sampling rate in A/pixel (double)
   public static final int MDL_SAMPLINGRATEZ = 141;  ///< sampling rate in A/pixel (double)
   public static final int MDL_SCALE = 142;  ///< scaling factor for an image or volume (double)
   public static final int MDL_SELFILE = 143;  ///< Name of an image (std::string)
   public static final int MDL_SERIE = 144;  ///< A collection of micrographs = 144;  e.g. a tilt serie (std::string)
   public static final int MDL_SHIFTX = 145;  ///< Shift for the image in the X axis (double)
   public static final int MDL_SHIFTY = 146;  ///< Shift for the image in the Y axis (double)
   public static final int MDL_SHIFTZ = 147;  ///< Shift for the image in the Z axis (double)
   public static final int MDL_SHIFT_CRYSTALX = 148;  ///< Shift for the image in the X axis (double) for crystals
   public static final int MDL_SHIFT_CRYSTALY = 149;  ///< Shift for the image in the Y axis (double) for crystals
   public static final int MDL_SHIFT_CRYSTALZ = 150;  ///< Shift for the image in the Z axis (double) for crystals
   public static final int MDL_SIGMANOISE = 151;  ///< Standard deviation of the noise in ML model
   public static final int MDL_SIGMAOFFSET = 152;  ///< Standard deviation of the offsets in ML model
   public static final int MDL_SIGNALCHANGE = 153;  ///< Signal change for an image
   public static final int MDL_SPHERICALABERRATION = 154;  ///<ctf definition azimutal angle
   public static final int MDL_STDDEV = 155;  ///<stdandard deviation value (double)
   public static final int MDL_SUM = 156;  ///< Sum of elements of a given type (double) [this is a genereic type do not use to transfer information to another program]
   public static final int MDL_SUMWEIGHT = 157;  ///< Sum of all weights in ML model
   public static final int MDL_SYMNO = 158;  ///< Symmetry number for a projection (used in ART)
   public static final int MDL_TRANSFORMATIONMTRIX = 159;  ///< transformation matrix(vector double)
   public static final int MDL_VOLTAGE = 160;  ///< microscope voltage (double)
   public static final int MDL_WEIGHT = 161;  ///< Weight assigned to the image (double)
   public static final int MDL_WROBUST = 162;  ///< Weight of t-student distribution in robust Maximum likelihood
   public static final int MDL_XINT = 163;  ///< X component (int)
   public static final int MDL_XINTTILT = 164;  ///< X component in tilted micrograph (int)
   public static final int MDL_X = 165;  ///< X component (double)
   public static final int MDL_YINT = 166;  ///< Y component (int)
   public static final int MDL_YINTTILT = 167;  ///< Y component in tilted micrograph (int)
   public static final int MDL_Y = 168;  ///< Y component (double)
   public static final int MDL_ZINT = 169;  ///< Z component (int)
   public static final int MDL_Z = 170;  ///< Z component (double)
   public static final int MDL_ZSCORE = 171;  ///< Z Score (double)
   public static final int MDL_LAST_LABEL = 172;                       // **** NOTE ****: Do keep this label always at the end
}
