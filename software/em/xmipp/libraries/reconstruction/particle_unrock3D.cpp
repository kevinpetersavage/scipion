/***************************************************************************
 *
 * Authors:   Vahid Abrishami         vabrishami@cnb.csic.es (2015)
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

#include "particle_unrock3D.h"
#include <data/numerical_tools.h>
#include <data/mask.h>
#include <math.h>
#include <time.h>
#include <data/xmipp_filename.h>

#include "data/xmipp_fftw.h"
#include <data/ctf.h>
#include <reconstruction/fourier_filter.h>

using namespace std;

// This function used by Powell to compute a cost function
double L1costImage(double *x, void *_prm)
{
    MultidimArray<double> transParticle;
    ProgParticlePolishing *prm=(ProgParticlePolishing *)_prm;

    // Generate transformation matrix
    double cosx3,sinx3;

    //[s*(cos(rot) + a*sin(rot) + a*b*cos(rot), -s*(sin(rot)-a*cos(rot)+a*b*sin(rot)), dx]
    //[s*(sin(rot) + b*cos(rot)), s*(cos(rot) - b*sin(rot)), dy ]
    //[0                        ,0                         , 1  ]
    //dx(shiftx)~x[1] dy(shifty)~x[2] rot~x[3] a(shearx)~x[4] b(sheary)~x[5] s(scale)~x[6]
    sincos(x[3],&sinx3,&cosx3);
    MAT_ELEM(prm->A,0,0)=x[6]*(cosx3+x[4]*sinx3+x[4]*x[5]*cosx3);
    MAT_ELEM(prm->A,0,1)=(-1.0*x[6])*(sinx3-x[4]*cosx3+x[4]*x[5]*sinx3);
    MAT_ELEM(prm->A,0,2)=x[1];
    MAT_ELEM(prm->A,1,0)=x[6]*(sinx3+x[5]*cosx3);
    MAT_ELEM(prm->A,1,1)=x[6]*(cosx3-x[5]*sinx3);
    MAT_ELEM(prm->A,1,2)=x[2];
    MAT_ELEM(prm->A,2,0)=0;
    MAT_ELEM(prm->A,2,1)=0;
    MAT_ELEM(prm->A,2,2)=1;

    //Check the constraints on each parameter
    double scale=x[6]-1;
    if (fabs(scale)>prm->scaleLimit
        || fabs(x[1])>prm->shiftLimit || fabs(x[2])>prm->shiftLimit
        || fabs(x[4])>prm->shearLimit || fabs(x[5])>prm->shearLimit
        || fabs(x[3]*(180.0/PI))>prm->rotLimit)
        return 1e38;
    // Apply the transformation to the current particle
    applyGeometry(LINEAR, transParticle, prm->avgParticle, prm->A, IS_NOT_INV, WRAP);

    //**********************************************************************************
    // Remove the comments to compute cost function using L1
    /*double cost=0, E;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(transParticle)
	{
        E = DIRECT_A2D_ELEM(transParticle,i,j) - DIRECT_A2D_ELEM(prm->avgParticle,i,j);
        cost+=fabs(E);
	}
    cost/=MULTIDIM_SIZE(transParticle);
    return cost;*/
    //***********************************************************************************

    //Computation of cost function using correlation
    double cost=(-1.0)*correlationIndex(prm->currentParticle, transParticle);
    return cost;
}

// Some useful procedures

// To convert integer to string (this function may not be required)
std::string int2Str(int num)
{
    std::ostringstream ostr;
    ostr << (num);
    return ostr.str();
}

// Merge two optical maps for X (flowX) and Y (flowY) into one map
void mergeFlows(cv::Mat flowX, cv::Mat flowY, cv::Mat& flow)
{
    float flowXVal,flowYVal;
    int cols = flowX.cols;
    int rows = flowX.rows;

    flow = cv::Mat::zeros(rows, cols, CV_32FC2);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
        	cv::Vec2f flow_at_point;
            flow_at_point[0] = flowX.at<float>(i, j);
            flow_at_point[1] = flowY.at<float>(i, j);
            flow.at<cv::Vec2f>(i, j) = flow_at_point;
        }
    }
}

// To save Mats of OpenCv on disk
int saveMat(const string& filename, const cv::Mat& M)
{
    if (M.empty())
    {
        return 0;
    }
    ofstream out(filename.c_str(), ios::out|ios::binary);
    if (!out)
        return 0;

    int cols = M.cols;
    int rows = M.rows;
    int chan = M.channels();
    int eSiz = (M.dataend-M.datastart)/(cols*rows*chan);

    // Write header
    out.write((char*)&cols,sizeof(cols));
    out.write((char*)&rows,sizeof(rows));
    out.write((char*)&chan,sizeof(chan));
    out.write((char*)&eSiz,sizeof(eSiz));

    // Write data.
    if (M.isContinuous())
    {
        out.write((char *)M.data,cols*rows*chan*eSiz);
    }
    else
    {
        return 0;
    }
    out.close();
    return 1;
}

// Load a matrix which is generated by saveMat
int readMat(const string& filename, cv::Mat& M)
{
    ifstream in(filename.c_str(), ios::in|ios::binary);
    if (!in)
    {
        //M = NULL_MATRIX;
        return 0;
    }
    int cols;
    int rows;
    int chan;
    int eSiz;

    // Read header
    in.read((char*)&cols,sizeof(cols));
    in.read((char*)&rows,sizeof(rows));
    in.read((char*)&chan,sizeof(chan));
    in.read((char*)&eSiz,sizeof(eSiz));

    // Determine type of the matrix
    int type = 0;
    switch (eSiz)
    {
    case sizeof(char):
                    type = CV_8UC(chan);
        break;
    case sizeof(float):
                    type = CV_32FC(chan);
        break;
    case sizeof(double):
                    type = CV_64FC(chan);
        break;
    }

    // Alocate Matrix.
    M = cv::Mat(rows,cols,type,cv::Scalar(1));

    // Read data.
    if (M.isContinuous())
{
        in.read((char *)M.data,cols*rows*chan*eSiz);
    }
    else
    {
        return 0;
    }
    in.close();
    return 1;
}

// Converts a XMIPP MultidimArray to OpenCV matrix
void xmipp2Opencv(const MultidimArray<double> &xmippArray, cv::Mat &opencvMat)
{
    int h = YSIZE(xmippArray);
    int w = XSIZE(xmippArray);
    opencvMat = cv::Mat::zeros(h, w,CV_32FC1);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
    opencvMat.at<float>(i,j) = DIRECT_A2D_ELEM(xmippArray,i,j);
}

// Converts an OpenCV matrix to XMIPP MultidimArray
void opencv2Xmipp(const cv::Mat &opencvMat, MultidimArray<double> &xmippArray)
{
    int h = opencvMat.rows;
    int w = opencvMat.cols;
    xmippArray.initZeros(h, w);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
    DIRECT_A2D_ELEM(xmippArray,i,j) = opencvMat.at<float>(i,j);
}

// Converts an OpenCV float matrix to an OpenCV Uint8 matrix
void convert2Uint8(cv::Mat opencvDoubleMat, cv::Mat &opencvUintMat)
{
    cv::Point minLoc,maxLoc;
    double min,max;
    cv::minMaxLoc(opencvDoubleMat, &min, &max, &minLoc, &maxLoc, cv::noArray());
    opencvDoubleMat.convertTo(opencvUintMat, CV_8U, 255.0/(max - min), -min * 255.0/(max - min));
}

// Implementation of the particle polishing class
// Constructor of the class
ProgParticlePolishing::ProgParticlePolishing()
{}

// Destructor of the class
ProgParticlePolishing::~ProgParticlePolishing()
{}

void ProgParticlePolishing::readParams()
{
    fnMovie = getParam("-i");
    fnPolishedParticles  = getParam("-o");
    particleSize = getIntParam("--particleSize");
    fnparticleCoords = getParam("--coords");
    w1 = getDoubleParam("--cutFre");
    fnInitFlow = getParam("--flows");
    shiftLimit = getDoubleParam("--shiftLim");
    scaleLimit = getDoubleParam("--scaleLim");
    shearLimit = getDoubleParam("--shearLim");
    rotLimit = getDoubleParam("--shearLim");
    iterationNum = getIntParam("--rotlimit");
}

void ProgParticlePolishing::defineParams()
{
    addUsageLine(
        "Polishes the selected particles from the corrected average using frames");
    addParamsLine("  -i <movie>               : Input movie");
    addParamsLine("  -o <outputstack>       : Polished particles");
    addParamsLine("  [--cutFre <cut_frequency=0.05>]       : Cut frequency for low pass filter");
    addParamsLine("  [--scaleLim <scale_limit=0.03>]       : Scale limit for transformation matrix");
    addParamsLine("  [--shearLim <shear_limit=0.01>]       : Shear limit for transformation matrix");
    addParamsLine("  [--shiftLim <shift_limit=1.5>]        : Shift limit for transformation matrix");
    addParamsLine("  [--rotlimit <rot_limit=1>]            : Rotation limit for transformation matrix");
    addParamsLine("  --particleSize <size>         : Particle size in pixels");
    addParamsLine("  --coords <coordsfile>         : Particles' coordinates");
    addParamsLine("  [--iterations <num_iteration=2>]      : Number of iteration for polishing particles");
    addParamsLine("  --flows <flow_rootname>      : Root name for the initial flows");
}

void ProgParticlePolishing::show()
{
    std::cerr<<"=========================================================="<<std::endl;
    std::cerr<<"Input movie: "<<fnMovie<<std::endl;
    std::cerr<<"Corrected average: "<<fnAvg<<std::endl;
    std::cerr<<"Particle coordinates: "<<fnparticleCoords<<std::endl;
    std::cerr<<"Particle size: "<<particleSize<<std::endl;
    std::cerr<<"=========================================================="<<std::endl;
}
void ProgParticlePolishing::produceSideInfo()
{
    ArrayDim aDim;

    // ***************** Variables are configured for phantom data ***************
    fnParticleStack = fnMovie + "_frame"; //e.g. Movie000053_frame10.stk
    fnAvg =  fnMovie + "_avg.stk"; //e.g. Movie000053_avg.stk (Reference frame)
    fnPolishedParticles = fnMovie + "_aligned.stk"; //e.g. Movie000053_frame10.stk
    frameNum = 16;
    particleNum = 50;
    // ***************************************************************************

    // ***************** Variables are configured for real data ******************
    //fnParticleStack = fnMovie.removeLastExtension() + "_frame";
    //fnAvg =  fnparticleCoords.removeLastExtension() + "_aligned.mrc";
    //fnPolishedParticles = fnparticleCoords.removeLastExtension() + "_aligned.stk";
    // Obtain the number of frames
    //movieStack.read(fnMovie+":mrcs",HEADER);
    //movieStack.read(fnMovie,HEADER);
    //movieStack.getDimensions(aDim);
    //frameNum = aDim.ndim;
    // Obtain the number of particles
    //particleCoords.read(fnparticleCoords);
    //particleNum = particleCoords.size();
    //******************************************************************************
    // Calculate the radius of the particle
    particleRadius=(size_t)(particleSize*0.5);
    extParticles.initZeros(particleNum, 1, particleSize, particleSize);
    particleAvgStack.initZeros(particleNum, 1, particleSize, particleSize);
    // Initialization of the transformation matrix
    A.initIdentity(3);
    transformationVec.resize(frameNum*particleNum);
}

// Extract a particle from a micrograph centered at x and y
void ProgParticlePolishing::micExtractParticle(const int x, const int y,
        const MultidimArray<double> frameImage,
        MultidimArray<double> &particleImage)
{
    int startX, startY, endX, endY;
    startX = x - particleRadius;
    startY = y - particleRadius;
    endX = x + particleRadius;
    endY = y + particleRadius;
    frameImage.window(particleImage, startY, startX, endY-1, endX-1);
}

// Compute the average of particle inside a group
void ProgParticlePolishing::computeAvgStack()
{
    FileName fnCurrentStack;
    MultidimArray<double> particleAvg, particleImg;
    ImageGeneric particleStack;

    for (size_t i=0;i<particleNum;i++)
    {
        // Set average to zero at the beginning of each iteration
        particleAvg.aliasImageInStack(particleAvgStack, i);
        particleAvg.initZeros(particleSize, particleSize);
        for (size_t j=0;j<frameNum;j++)
        {
            fnCurrentStack = fnParticleStack + int2Str(j+1) + ".stk";
            particleStack.readMapped(fnCurrentStack, i+1);
            particleStack().getImage(particleImg);
            particleAvg += particleImg;
        }
        particleAvg /= double(frameNum);
    }
}

/*
 * The goal of this function is to estimate the position
 * of instances of each particle within each frame respect
 * to the position of the corrected average (from OF) by computing
 * optical flow between the average and each frame.
 */
void ProgParticlePolishing::extractParticels()
{
    int x, y, shiftedX, shiftedY;
    // For reading each frame from the movie
    MultidimArray<double> frameImg, particleImg, avgImage;
    cv::Mat frameImgOpenCV, avgImgOpenCV;
    // Stack for particles of a frame
    MultidimArray<double> particlesStack, avgParticleImage;
    Image<double> Iaux;
    CTFDescription ctf;
    // Stack filename for particles in a frame
    FileName fnCurrentStack, fnCTFDescr;
    FileName fnFlowX, fnFlowY;

    // OpenCV structures
    // Matrices for computing optical flow
    cv::Mat flowX, flowY, flow;
    // Equal 8 bit images for average and current frame for optical flow
    cv::Mat frameImg8, avgImg8;

    std::cerr<<"Step 1: Extracting particles from the frames"<<std::endl;
    // Read the OF corrected average
    Iaux.read(fnAvg);
    avgImage = Iaux();
    // Convert the XMIPP array format to OpenCV matrix format
    xmipp2Opencv(avgImage, avgImgOpenCV);
    convert2Uint8(avgImgOpenCV, avgImg8);
    // Read each frame of the movie
    for (size_t i=0;i<frameNum;i++)
    {
        std::cerr<<"Extracting particles from frame "<<i+1<<std::endl;
        // Allocate memory for the particles in the current frame
        particlesStack.resize(particleNum, 1, particleSize, particleSize);
        //movieStack.readMapped(fnMovie+":mrcs",i+1);
        movieStack.readMapped(fnMovie,i+1);
        movieStack().getImage(frameImg);
        xmipp2Opencv(frameImg, frameImgOpenCV);
        convert2Uint8(frameImgOpenCV, frameImg8);
        // Compute the shifts between the corrected average and current frame
        cv::calcOpticalFlowFarneback(avgImg8, frameImg8, flow, 0.5, 6, 150, 1, 5, 1.1, 0);
        // Read particles of each frame
        size_t particleIdx = 0;
        // Prepare the stack file name to put corrected particles
        fnCurrentStack = fnParticleStack + int2Str(i+1) + ".stk";

        /****************************************************************************

         * Remove the comments of you want to apply phase flip
        //fnCTFDescr = fnMovie.removeLastExtension() + "-aligned/xmipp_ctf.ctfparam";
        //ctf.clear();
        //ctf.read(fnCTFDescr);
        //ctf.produceSideInfo();
        actualPhaseFlip(frameImg,ctf);
        ******************************************************************************/
        // For the all particles in the current frame
        FOR_ALL_OBJECTS_IN_METADATA(particleCoords)
        {
            // The position of the particle in the average image
            particleCoords.getValue(MDL_XCOOR, x, __iter.objId);
            particleCoords.getValue(MDL_YCOOR, y, __iter.objId);
            cv::Vec2f flow_at_point = flow.at<cv::Vec2f>(y, x);
            // Compute the real position of the particle in the frame
            shiftedX = round(x + flow_at_point[0]);
            shiftedY = round(y + flow_at_point[1]);
            // Extract particle from the current frame at shiftedX, shiftedY
            particleImg.aliasImageInStack(particlesStack, particleIdx);
            micExtractParticle(shiftedX, shiftedY, frameImg, particleImg);
            particleIdx++;
        }
        Iaux() = particlesStack;
        Iaux.write(fnCurrentStack,ALL_IMAGES,true,WRITE_APPEND);
    }
    std::cerr<<"Extracting particles from the frames has been done"<<std::endl;
}

void ProgParticlePolishing::computeGroupAvg(size_t start, size_t end, size_t particleIdx, MultidimArray<double> &groupAverage)
{
    FileName fnCurrentStack;
    ImageGeneric particleStack;
    MultidimArray<double> particleOF;
    size_t particleNum;

    groupAverage.initZeros(particleSize, particleSize);
    particleNum=end-start;
    for (size_t frameNum=start;frameNum<end;frameNum++)
    {

        cout<<"Frame number is: "<<frameNum+1<<endl;
        fnCurrentStack = fnParticleStack + int2Str(frameNum+1) + ".stk";
        particleStack.readMapped(fnCurrentStack, particleIdx+1);
        particleStack().getImage(particleOF);
        groupAverage+=particleOF;
    }
    groupAverage/=particleNum;
}

void ProgParticlePolishing::particlePolishing3D()
{
    FileName fnCurrentStack, fnCurrentStackAl;
    ImageGeneric particleStack, avgStack;
    Matrix1D<double> p, steps;
    MultidimArray<double> particleImg, transparticleImg, particleAvg, frameImg;
    MultidimArray<double> groupAvg;
    MultidimArray<int> mask;
    Image<double> avgImg, Iaux;
    size_t particleIdx;
    FourierFilter filter;
    int iter, x, y, divCoeff=2;
    int groupSize, cntFrame;
    double xScale, yScale, ccSum;
    double xShear, yShear;
    double cosp2,sinp2;

    cv::Mat avgParticleCV, currentParticleCV, warp_matrix, groupAvgCV;
    cv::Mat currentParticleAlignedCV;
    cv::Mat_<cv::Point2f> flow;
    // Initialization
    steps.initZeros(6);
    steps.initConstant(1.);
    avgImg.read(fnAvg);
    filter.FilterShape=RAISED_COSINE;
    filter.FilterBand=LOWPASS;
    filter.w1=w1;
    particleIdx = 0;
    mask.resize(particleSize,particleSize);
    mask.setXmippOrigin();
    BinaryCircularMask(mask, particleRadius);
    MultidimArray<double> Maux(particleSize, particleSize);
    for (size_t pt=0;pt<particleNum;pt++)
    {
        cout<<"We are processing the particle: "<<particleIdx+1<<endl;
        avgStack.readMapped(fnAvg, pt+1);
        //avgParticle has the current ref for unrock3D
        avgStack().getImage(avgParticle);
        //ParticleAvg will have the aligned particle
        particleAvg.aliasImageInStack(extParticles, pt);
        particleAvg.initZeros(particleSize, particleSize);
        avgParticle.setXmippOrigin();
        //We just mask the reference
        apply_binary_mask(mask, avgParticle, avgParticle);
        xmipp2Opencv(avgParticle, avgParticleCV);
        cntFrame=0;
        p.initZeros(6);
        p(0) = p(4) = 1.;
        for (size_t i=0;i<frameNum;i++)
        {
            cout<<"Frame number is: "<<i+1<<endl;
            fnCurrentStack = fnParticleStack + int2Str(i+1) + ".stk";
            particleStack.readMapped(fnCurrentStack, pt+1);
            particleStack().getImage(currentParticle);
            //To have the real particle after filtering
            groupAvg=currentParticle;
            currentParticle.setXmippOrigin();
            apply_binary_mask(mask, currentParticle, currentParticle);
            filter.applyMaskSpace(currentParticle);
            apply_binary_mask(mask, currentParticle, currentParticle);
            p.initZeros(6);
            p(5) = 1.;
            double cost=1e38;
            //Compute free transformation
            powellOptimizer(p, 1, 6, &L1costImage, this, 0.01, cost, iter, steps, true);
            if (cost == 1e38)
            	A.initIdentity(3);
            else
            {
            	// Transformation matrix
            	sincos(p(2),&sinp2,&cosp2);
            	MAT_ELEM(A,0,0)=p(5)*(cosp2+p(3)*sinp2+p(3)*p(4)*cosp2);
                MAT_ELEM(A,0,1)=(-1.0*p(5))*(sinp2-p(3)*cosp2+p(3)*p(4)*sinp2);
                MAT_ELEM(A,0,2)=p(0);
                MAT_ELEM(A,1,0)=p(5)*(sinp2+p(4)*cosp2);
                MAT_ELEM(A,1,1)=p(5)*(cosp2-p(4)*sinp2);
                MAT_ELEM(A,1,2)=p(1);
				MAT_ELEM(A,2,0)=0;
				MAT_ELEM(A,2,1)=0;
				MAT_ELEM(A,2,2)=1;
            }
            cntFrame++;
            std::cout<<"The transformatj¡ion matrix is :"<<A<<std::endl;
            applyGeometry(BSPLINE3, transparticleImg, groupAvg, A, IS_INV, WRAP);
            particleAvg+=transparticleImg;
        }
        particleAvg/=cntFrame;
        particleIdx++;
    }
    // Write the aligned particles
    Iaux() = extParticles;
    Iaux.write(fnPolishedParticles, ALL_IMAGES, true, WRITE_APPEND);
}

void ProgParticlePolishing::extractAverageParticle()
{
    int x,y;
    Image<double> avgImg, Iaux;
    avgImg.read(fnAvg);

    FOR_ALL_OBJECTS_IN_METADATA(particleCoords)
    {
        particleCoords.getValue(MDL_XCOOR, x, __iter.objId);
        particleCoords.getValue(MDL_YCOOR, y, __iter.objId);
        micExtractParticle(x, y, avgImg(), avgParticle);
        Iaux()=avgParticle;
        Iaux.write("testrealavg.stk",ALL_IMAGES,true,WRITE_APPEND);
    }
}

void ProgParticlePolishing::computeRadialAvg()
{
    ImageGeneric particleStack;
    FileName fnCurrentStack;
    MultidimArray<std::complex<double> > Faux;
    MultidimArray<double> Maux, rmean_signal, SMaux;
    MultidimArray<int> radial_count;
    Matrix1D<int> center(2);

    center.initZeros();
    for (size_t i=0;i<frameNum;i++)
    {
        cout<<"*********************************"<<endl;
        cout<<"Frame under processing is"<<std::endl;
        fnCurrentStack = fnParticleStack + int2Str(i+1) + ".stk";
        SMaux.initZeros(particleSize, particleSize);
        SMaux.setXmippOrigin();
        for (size_t j=0;j<particleNum;j++)
        {
            particleStack.readMapped(fnCurrentStack, j+1);
            particleStack().getImage(currentParticle);
            FourierTransform(currentParticle, Faux);
            FFT_magnitude(Faux, Maux);
            CenterFFT(Maux, true);
            Maux *= Maux;
            Maux.setXmippOrigin();
            SMaux += Maux;
        }
        SMaux /= particleNum;
        rmean_signal.initZeros();
        radialAverage(SMaux, center, rmean_signal, radial_count, true);
        cout<<rmean_signal<<endl;
    }
}

void ProgParticlePolishing::run()
{
    clock_t startTime;
    //Image<double> avgImages;
    produceSideInfo();
    show();
    //extractAverageParticle();
    //startTime = clock();
    //extractParticels();
    //printf("Time for extracting particles: %.2fs\n", (double)(clock() - startTime)/CLOCKS_PER_SEC);
    //computeAvgStack();
    //avgImages()=particleAvgStack;
    //avgImages.write("testavgimag.mrc");
    startTime = clock();
    particlePolishing3D();
    printf("Time for polishing particles: %.2fs\n", (double)(clock() - startTime)/CLOCKS_PER_SEC);
    //    std::vector<Mat> warpGround;
    //    std::vector< Mat > simulParticles;
    //    Image<double> II;
    //    MultidimArray<double> refImage, transImage,noisyImage;
    //    MultidimArray<double> refImageNormalize, tranBackImage, avgImage;
    //    Mat refImageCV, transImageCV, warp_matrix, invert_warp_matrix;
    //    Mat refImageNormalizeCV, tranBackImageCV, noisyImageCV, avgImageCV;
    //    warpGround.resize(10);
    //    simulParticles.resize(10);
    //
    //    warpGround[0]=(Mat_<float>(3,3) << 0.99517697, -0.002530592, 1.2480022, 0.0019662017, 0.99426794, 0.05981759, 0, 0, 1);
    //    warpGround[1]=(Mat_<float>(3,3) << 0.99560261, 0.0023737659, 0.44687042, -0.00074027939, 0.99518883, 1.108709, 0, 0, 1);
    //    warpGround[2]=(Mat_<float>(3,3) << 1.0042983,  0.0040052789, -0.8620984, 0.0015993031, 1.0078131, -1.1137034, 0, 0, 1);
    //    warpGround[3]=(Mat_<float>(3,3) << 0.99628067, -0.0016642481,  0.86237586, -0.0031945514, 0.99755478, 1.0432063, 0, 0, 1);
    //    warpGround[4]=(Mat_<float>(3,3) << 0.9977138, 0.0011662564, 0.17499633, 0.0011602666, 0.99887437, -0.56903493, 0, 0, 1);
    //    warpGround[5]=(Mat_<float>(3,3) << 1.0017791, -0.0049244128, 0.1899711,0.0023116327, 1.0010291, -0.2773093, 0, 0, 1);
    //    warpGround[6]=(Mat_<float>(3,3) << 1.0073326, 0.0040010153, -0.86902905, -0.00066119432, 0.99547476, 0.26132536, 0, 0, 1);
    //    warpGround[7]=(Mat_<float>(3,3) << 1.0005574, 0.0056234887, -1.1801665, -0.0047680577, 1.0032096, 0.14632595, 0, 0, 1);
    //    warpGround[8]=(Mat_<float>(3,3) << 1.0031086, 0.0025546292, -0.96174616, 0.0073704273, 0.99284703, 0.51307178, 0, 0, 1);
    //    warpGround[9]=(Mat_<float>(3,3) << 1.0039386, 0.0014224659, -0.83403122, 0.0049594352, 0.99990237,  -1.1123345, 0, 0, 1);
    //
    //    // Initializing the variables
    //    A.initIdentity(3);
    //    II.read("refimage.xmp");
    //    refImage=II();
    //    xmipp2Opencv(refImage, refImageCV);
    //    cv::Mat noise = Mat(refImageCV.size(),CV_32F);
    //    normalize(refImageCV, refImageNormalizeCV, 0.0, 1.0, CV_MINMAX, CV_32F);
    //    opencv2Xmipp(refImageNormalizeCV, refImageNormalize);
    //    noisyImage=refImageNormalize;
    //    tranBackImage=refImageNormalize;
    //    avgImage.initZeros(refImageNormalize);
    //
    //    // Generating frames
    //    for (size_t i=0;i<10;i++)
    //    {
    //        MAT_ELEM(A,0,0)=warpGround[i].at<float>(0,0);
    //        MAT_ELEM(A,0,1)=warpGround[i].at<float>(0,1);
    //        MAT_ELEM(A,0,2)=warpGround[i].at<float>(0,2);
    //        MAT_ELEM(A,1,0)=warpGround[i].at<float>(1,0);
    //        MAT_ELEM(A,1,1)=warpGround[i].at<float>(1,1);
    //        MAT_ELEM(A,1,2)=warpGround[i].at<float>(1,2);
    //        MAT_ELEM(A,2,0)=0;
    //        MAT_ELEM(A,2,1)=0;
    //        MAT_ELEM(A,2,2)=1;
    //
    //        cv::randn(noise, 0, 0.97);
    //        applyGeometry(BSPLINE3, noisyImage, refImageNormalize, A, IS_INV, WRAP);
    //        xmipp2Opencv(noisyImage, noisyImageCV);
    //        noisyImageCV=noisyImageCV+noise;
    //        normalize(noisyImageCV, noisyImageCV, 0.0, 1.0, CV_MINMAX, CV_32F);
    //        opencv2Xmipp(noisyImageCV,noisyImage);
    //        simulParticles[i].push_back(noisyImageCV);
    //        avgImage+=noisyImage;
    //        II()=noisyImage;
    //        II.write("noisyframes.stk",ALL_IMAGES,true,WRITE_APPEND);
    //    }
    //    avgImage/=10.0;
    //    xmipp2Opencv(avgImage, avgImageCV);
    //    normalize(avgImageCV, avgImageCV, 0.0, 1.0, CV_MINMAX, CV_32F);
    //    opencv2Xmipp(avgImageCV, avgImage);
    //    II()=avgImage;
    //    II.write("avgframe.xmp");
    //    // End of generating frames
    //
    //    // Aligning of the frames
    //    for (size_t i=0;i<10;i++)
    //    {
    //        warp_matrix=Mat::eye(2, 3, CV_32FC1);
    //        //cv::Scalar avgPixelIntensity = cv::mean(abs(simulParticles[i]));
    //        //cout<<"The average of the matrix is "<<avgPixelIntensity<<std::endl;
    //        double cc = findTransformECC(avgImageCV,  simulParticles[i], warp_matrix, MOTION_AFFINE);
    //        //cout<<"The correlation value is "<<cc<<std::endl;
    //
    //        MAT_ELEM(A,0,0)=warp_matrix.at<float>(0,0);
    //        MAT_ELEM(A,0,1)=warp_matrix.at<float>(0,1);
    //        MAT_ELEM(A,0,2)=warp_matrix.at<float>(0,2);
    //        MAT_ELEM(A,1,0)=warp_matrix.at<float>(1,0);
    //        MAT_ELEM(A,1,1)=warp_matrix.at<float>(1,1);
    //        MAT_ELEM(A,1,2)=warp_matrix.at<float>(1,2);
    //        MAT_ELEM(A,2,0)=0;
    //        MAT_ELEM(A,2,1)=0;
    //        MAT_ELEM(A,2,2)=1;
    //
    //        opencv2Xmipp(simulParticles[i], transImage);
    //        applyGeometry(BSPLINE3, tranBackImage, transImage, A, IS_INV, WRAP);
    //        II()=tranBackImage;
    //        II.write("transnBacbackimage.stk",ALL_IMAGES,true,WRITE_APPEND);
    //        //invertAffineTransform(warp_matrix, invert_warp_matrix);
    //        invertAffineTransform(warp_matrix, invert_warp_matrix);
    //        Mat warp_matrix_tmp=Mat::eye(3, 3, CV_32FC1);
    //        warp_matrix_tmp.at<float>(0,0)=invert_warp_matrix.at<float>(0,0);
    //        warp_matrix_tmp.at<float>(0,1)=invert_warp_matrix.at<float>(0,1);
    //        warp_matrix_tmp.at<float>(0,2)=invert_warp_matrix.at<float>(0,2);
    //        warp_matrix_tmp.at<float>(1,0)=invert_warp_matrix.at<float>(1,0);
    //        warp_matrix_tmp.at<float>(1,1)=invert_warp_matrix.at<float>(1,1);
    //        warp_matrix_tmp.at<float>(1,2)=invert_warp_matrix.at<float>(1,2);
    //        warp_matrix_tmp.at<float>(2,0)=0;
    //        warp_matrix_tmp.at<float>(2,1)=0;
    //        warp_matrix_tmp.at<float>(2,2)=1;
    //        cout<<"original transformation is "<<std::endl<<warpGround[i]<<std::endl;
    //        cout<<"warp matrix is "<<std::endl<<invert_warp_matrix<<std::endl;
    //        Mat AA;
    //        gemm(warpGround[i], warp_matrix_tmp, 1.0, NULL, 0, AA);
    //        cout<<"final error is "<<AA<<std::endl;
    //    }

    // End of aligning frames

    //
    //cv::randn(noise, 0, 0.97);
    //    normalize(refImageCV, refImageNormalizeCV, 0.0, 1.0, CV_MINMAX, CV_32F);
    //    refImageNormalizeCV=refImageNormalizeCV+noise;
    //    normalize(refImageNormalizeCV, refImageNormalizeCV, 0.0, 1.0, CV_MINMAX, CV_32F);
    //    opencv2Xmipp(refImageNormalizeCV, refImageNormalize);
    //    II()=refImageNormalize;
    //    II.write("noisyref.xmp");
    //    transImage=refImageNormalize;
    //    tranBackImage=refImageNormalize;
    //    A.initIdentity(3);
    //    for (size_t i=0;i<10;i++)
    //    {
    //        cout<<"Initializing Matrix"<<std::endl;
    //        warp_matrix = Mat::eye(2, 3, CV_32FC1)
    //                      ;
    //        MAT_ELEM(A,0,0)=warpGround[i].at<float>(0,0);
    //        MAT_ELEM(A,0,1)=warpGround[i].at<float>(0,1);
    //        MAT_ELEM(A,0,2)=warpGround[i].at<float>(0,2);
    //        MAT_ELEM(A,1,0)=warpGround[i].at<float>(1,0);
    //        MAT_ELEM(A,1,1)=warpGround[i].at<float>(1,1);
    //        MAT_ELEM(A,1,2)=warpGround[i].at<float>(1,2);
    //        MAT_ELEM(A,2,0)=0;
    //        MAT_ELEM(A,2,1)=0;
    //        MAT_ELEM(A,2,2)=1;
    //        cout<<"Initializing Matrix"<<std::endl;
    //        applyGeometry(BSPLINE3, transImage, refImageNormalize, A, IS_INV, WRAP);
    //        cout<<"Matrix applied"<<std::endl;
    //warpAffine(refImageNormalizeCV, transImageCV, warpGround[i], Size(YSIZE(refImage),XSIZE(refImage)), INTER_CUBIC + WARP_INVERSE_MAP);
    //opencv2Xmipp(transImageCV,transImage);
    //        II()=transImage;
    //        II.write("transimage.stk",ALL_IMAGES,true,WRITE_APPEND);
    //        xmipp2Opencv(transImage, transImageCV);
    //        double cc = findTransformECC(refImageNormalizeCV, transImageCV, warp_matrix, MOTION_AFFINE);
    //warpAffine(transImageCV, tranBackImageCV, warp_matrix, Size(YSIZE(refImage),XSIZE(refImage)), INTER_CUBIC + WARP_INVERSE_MAP);
    //        MAT_ELEM(A,0,0)=warp_matrix.at<float>(0,0);
    //        MAT_ELEM(A,0,1)=warp_matrix.at<float>(0,1);
    //        MAT_ELEM(A,0,2)=warp_matrix.at<float>(0,2);
    //        MAT_ELEM(A,1,0)=warp_matrix.at<float>(1,0);
    //        MAT_ELEM(A,1,1)=warp_matrix.at<float>(1,1);
    //        MAT_ELEM(A,1,2)=warp_matrix.at<float>(1,2);
    //        MAT_ELEM(A,2,0)=0;
    //        MAT_ELEM(A,2,1)=0;
    //        MAT_ELEM(A,2,2)=1;
    //        MAT_ELEM(A,0,0)=1;
    //        MAT_ELEM(A,0,1)=0;
    //        MAT_ELEM(A,0,2)=0;
    //        MAT_ELEM(A,1,0)=0;
    //        MAT_ELEM(A,1,1)=1;
    //        MAT_ELEM(A,1,2)=0;
    //        MAT_ELEM(A,2,0)=0;
    //        MAT_ELEM(A,2,1)=0;
    //        MAT_ELEM(A,2,2)=1;
    //        applyGeometry(BSPLINE3, tranBackImage, transImage, A, IS_INV, WRAP);
    //        II()=tranBackImage;
    //        II.write("transnBacbackimage.stk",ALL_IMAGES,true,WRITE_APPEND);
    //        xmipp2Opencv(tranBackImage, tranBackImageCV);
    //        Mat errorImage;
    //        subtract(refImageNormalizeCV, tranBackImageCV, errorImage);
    //        double max_of_error;
    //        cv::Scalar avgPixelIntensity = cv::mean(abs(errorImage));
    //        std::cout<<avgPixelIntensity<<std::endl;
    //
    //        invertAffineTransform(warp_matrix, invert_warp_matrix);
    //        cout<<"Transformation matrix is "<<endl<<invert_warp_matrix<<endl;
    //  }

}
