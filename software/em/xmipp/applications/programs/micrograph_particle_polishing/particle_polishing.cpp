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

#include "particle_polishing.h"
#include <data/numerical_tools.h>
#include <math.h>
#include <time.h>
#include <data/xmipp_filename.h>
#include <reconstruction/ctf_phase_flip.h>
#include "data/xmipp_fftw.h"

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/video/video.hpp"

using namespace std;
using namespace cv;

double L1costImage(double *x, void *_prm)
{
    MultidimArray<double> transParticle;
    ProgParticlePolishing *prm=(ProgParticlePolishing *)_prm;
    // Transformation matrix
    MAT_ELEM(prm->A,0,0)=x[1];
    MAT_ELEM(prm->A,0,1)=x[2];
    MAT_ELEM(prm->A,0,2)=x[3];
    MAT_ELEM(prm->A,1,0)=x[4];
    MAT_ELEM(prm->A,1,1)=x[5];
    MAT_ELEM(prm->A,1,2)=x[6];
    MAT_ELEM(prm->A,2,0)=0;
    MAT_ELEM(prm->A,2,1)=0;
    MAT_ELEM(prm->A,2,2)=1;

    // Check the limits
    double xScale, yScale;
    double xShear, yShear;
    xScale=x[1]-1;
    yScale=x[5]-1;
    xShear=x[2];
    yShear=x[4];
    if (fabs(xScale)>prm->scaleLimit || fabs(yScale)>prm->scaleLimit
        || fabs(x[3])>prm->shiftLimit || fabs(x[6])>prm->shiftLimit
        || fabs(xShear)>prm->shearLimit || fabs(yShear)>prm->shearLimit)
        return 1e38;
    // Apply the transformation to the current particle
    applyGeometry(LINEAR, transParticle, prm->currentParticle, prm->A, IS_NOT_INV, WRAP);

    // Compute the cost
    double cost=0, E;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(transParticle)
    {
        E = DIRECT_A2D_ELEM(transParticle,i,j) - DIRECT_A2D_ELEM(prm->avgParticle,i,j);
        cost+=fabs(E);
    }
    cost/=MULTIDIM_SIZE(transParticle);
    //cout<<"Transformation matrix in L1cost is :"<<prm->A<<endl;
    //cout<<"Cost is: "<<cost<<endl;
    return cost;
}

// Some useful procedures
std::string int2Str(int num)
{
    std::ostringstream ostr;
    ostr << (num);
    return ostr.str();
}

void mergeFlows(Mat flowX, Mat flowY, Mat& flow)
{
    float flowXVal,flowYVal;
    int cols = flowX.cols;
    int rows = flowX.rows;

    flow = Mat::zeros(rows, cols, CV_32FC2);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            Vec2f flow_at_point;
            flow_at_point[0] = flowX.at<float>(i, j);
            flow_at_point[1] = flowY.at<float>(i, j);
            flow.at<Vec2f>(i, j) = flow_at_point;
        }
    }
}

int saveMat(const string& filename, const Mat& M)
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
int readMat(const string& filename, Mat& M)
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
    M = Mat(rows,cols,type,Scalar(1));

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
void xmipp2Opencv(const MultidimArray<double> &xmippArray, Mat &opencvMat)
{
    int h = YSIZE(xmippArray);
    int w = XSIZE(xmippArray);
    opencvMat = cv::Mat::zeros(h, w,CV_32FC1);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
    opencvMat.at<float>(i,j) = DIRECT_A2D_ELEM(xmippArray,i,j);
}

// Converts an OpenCV matrix to XMIPP MultidimArray
void opencv2Xmipp(const Mat &opencvMat, MultidimArray<double> &xmippArray)
{
    int h = opencvMat.rows;
    int w = opencvMat.cols;
    xmippArray.initZeros(h, w);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
    DIRECT_A2D_ELEM(xmippArray,i,j) = opencvMat.at<float>(i,j);
}

// Converts an OpenCV float matrix to an OpenCV Uint8 matrix
void convert2Uint8(Mat opencvDoubleMat, Mat &opencvUintMat)
{
    cv::Point minLoc,maxLoc;
    double min,max;
    cv::minMaxLoc(opencvDoubleMat, &min, &max, &minLoc, &maxLoc, noArray());
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
    fnInitFlow = getParam("--flows");
    shiftLimit = getDoubleParam("--shiftLim");
    scaleLimit = getDoubleParam("--scaleLim");
    shearLimit = getDoubleParam("--shearLim");
    iterationNum = getIntParam("--iterations");
}

void ProgParticlePolishing::defineParams()
{
    addUsageLine(
        "Polishes the selected particles from the corrected average using frames");
    addParamsLine("  -i <movie>               : Input movie");
    addParamsLine("  -o <outputstack>       : Polished particles");
    addParamsLine("  [--scaleLim <scale_limit=0.03>]       : Scale limit for transformation matrix");
    addParamsLine("  [--shearLim <shear_limit=0.03>]       : Shear limit for transformation matrix");
    addParamsLine("  [--shiftLim <shift_limit=3.0>]       : Shift limit for transformation matrix");
    addParamsLine("  --particleSize <size>         : Particle size in pixels");
    addParamsLine("  --coords <coordsfile>         : Particles' coordinates");
    addParamsLine("  [--iterations <num_iteration=2>]      : Number of iteration for polishing particles");
    addParamsLine("  --flows <flow_rootname>      : Root name for the initial flows");
    //    addParamsLine("  [--thr <p=1>]                 : Number of threads for automatic picking");
    //    addParamsLine("  [--fast]                      : Perform a fast preprocessing of the micrograph (Fourier filter instead of Wavelet filter)");
    //    addParamsLine("  [--in_core]                   : Read the micrograph in memory");
    //    addParamsLine("  [--filter_num <n=6>]          : The number of filters in filter bank");
    //    addParamsLine("  [--NPCA <n=4>]               : The number of PCA components");
    //    addParamsLine("  [--NCORR <n=2>]               : The number of PCA components");
    //    addParamsLine("  [--autoPercent <n=90>]               : The number of PCA components");
    //    addExampleLine("Automatically select particles during training:", false);
    //    addExampleLine("xmipp_micrograph_automatic_picking -i micrograph.tif --particleSize 100 --model model --thr 4 --outputRoot micrograph --mode try ");
    //    addExampleLine("Training:", false);
    //    addExampleLine("xmipp_micrograph_automatic_picking -i micrograph.tif --particleSize 100 --model model --thr 4 --outputRoot micrograph --mode train manual.pos");
    //    addExampleLine("Automatically select particles after training:", false);
    //    addExampleLine("xmipp_micrograph_automatic_picking -i micrograph.tif --particleSize 100 --model model --thr 4 --outputRoot micrograph --mode autoselect");
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

    fnParticleStack = fnMovie.removeLastExtension() + "_frame";
    fnAvg =  fnparticleCoords.removeLastExtension() + "-aligned.mrc";
    // Obtain the number of frames
    movieStack.read(fnMovie,HEADER);
    movieStack.getDimensions(aDim);
    frameNum = aDim.ndim;
    // Obtain the number of particles
    particleCoords.read(fnparticleCoords);
    particleNum = particleCoords.size();
    // Calculate the radius of the particle
    particleRadius=(size_t)(particleSize*0.5);
    extParticles.initZeros(particleNum, 1, particleSize, particleSize);
    particleAvgStack.initZeros(particleNum, 1, particleSize, particleSize);
    // Initialization of the transformation matrix
    A.initIdentity(3);
    transformationVec.resize(frameNum*particleNum);
}

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

void ProgParticlePolishing::extractParticels()
{
    int x, y, shiftedX, shiftedY;
    // For reading each frame from the movie
    MultidimArray<double> frameImg, particleImg, avgImage;
    Mat frameImgOpenCV, avgImgOpenCV;
    // Stack for particles of a frame
    MultidimArray<double> particlesStack, avgParticleImage;
    Image<double> Iaux;
    CTFDescription ctf;
    // Stack filename for particles in a frame
    FileName fnCurrentStack, fnCTFDescr;
    FileName fnFlowX, fnFlowY;

    // OpenCV structures
    // Matrices for computing optical flow
    Mat flowX, flowY, flow;
    // Equal 8 bit images for average and current frame for optical flow
    Mat frameImg8, avgImg8;

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
        movieStack.readMapped(fnMovie,i+1);
        movieStack().getImage(frameImg);
        xmipp2Opencv(frameImg, frameImgOpenCV);
        convert2Uint8(frameImgOpenCV, frameImg8);

        /* Reload the initial shifts between the current frame and the average
        fnFlowX = fnInitFlow + "x" + int2Str(i) + ".mat";
        fnFlowY = fnInitFlow + "y" + int2Str(i) + ".mat";
        readMat(fnFlowX.c_str(), flowX);
        readMat(fnFlowY.c_str(), flowY);
        mergeFlows(flowX, flowY, flow);*/

        // Compute the shifts between the corrected average and current frame
        calcOpticalFlowFarneback(avgImg8, frameImg8, flow, 0.5, 6, 150, 1, 5, 1.1, 0);
        // Read particles of each frame
        size_t particleIdx = 0;
        fnCurrentStack = fnParticleStack + int2Str(i+1) + ".stk";
        fnCTFDescr = fnMovie.removeLastExtension() + "-aligned/xmipp_ctf.ctfparam";
        ctf.clear();
        ctf.read(fnCTFDescr);
        ctf.produceSideInfo();
        actualPhaseFlip(frameImg,ctf);
        FOR_ALL_OBJECTS_IN_METADATA(particleCoords)
        {
            // The position of the particle in the average image
            particleCoords.getValue(MDL_XCOOR, x, __iter.objId);
            particleCoords.getValue(MDL_YCOOR, y, __iter.objId);
            cout<<x<<" "<<y<<endl;
            // Shifts in x and y at the particle position
            //Point2f u = flow(y, x);
            Vec2f flow_at_point = flow.at<Vec2f>(y, x);
            shiftedX = round(x + flow_at_point[0]);
            shiftedY = round(y + flow_at_point[1]);
            //std::cerr<<"The value of x is "<<x<<" and the value of y is "<<y<<std::endl;
            //std::cerr<<"The value of x after shift is "<<shiftedX<<" and the value of y is "<<shiftedY<<std::endl;
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

void ProgParticlePolishing::particlePolishing3D()
{
    FileName fnCurrentStack, fnCurrentStackAl;
    ImageGeneric particleStack;
    Matrix1D<double> p, steps;
    MultidimArray<double> particleImg, transparticleImg, particleAvg;
    Image<double> avgImg, Iaux;
    size_t particleIdx;
    int iter, x, y;

    // Initialization
    steps.initZeros(6);
    steps.initConstant(1.);
    avgImg.read(fnAvg);
    for (size_t it=0;it<iterationNum;it++)
    {
        cout<<"Iteration number: "<<it<<endl;
        particleIdx = 0;
        //if (fnCurrentStackAl.exists())
        //fnCurrentStackAl.deleteFile();
        FOR_ALL_OBJECTS_IN_METADATA(particleCoords)
        {
            cout<<"****************************************"<<endl;
            cout<<"We are processing the particle: "<<particleIdx+1<<endl;
            if (it == 0)
            {
                particleCoords.getValue(MDL_XCOOR, x, __iter.objId);
                particleCoords.getValue(MDL_YCOOR, y, __iter.objId);
                micExtractParticle(x, y, avgImg(), avgParticle);
            }
            else
                extParticles.getImage(particleIdx, avgParticle);
            particleAvg.aliasImageInStack(extParticles, particleIdx);
            particleAvg.initZeros(particleSize, particleSize);
            for (size_t i=0;i<frameNum;i++)
            {
                cout<<"Frame number is: "<<i+1<<endl;
                fnCurrentStack = fnParticleStack + int2Str(i+1) + ".stk";
                fnCurrentStackAl = "stack" + int2Str(i+1) + ".stk";
                particleStack.readMapped(fnCurrentStack, particleIdx+1);
                particleStack().getImage(currentParticle);
                if (it == 0)
                {
                    p.initZeros(6);
                    p(0) = p(4) = 1.;
                }
                else
                    p = transformationVec[particleIdx*frameNum+i];
                double cost=1e38;
                // Compute free transformation
                powellOptimizer(p, 1, 6, &L1costImage, this, 0.01, cost, iter, steps, true);
                if (cost == 1e38)
                    A.initIdentity(3);
                else
                    // Transformation matrix
                    MAT_ELEM(A,0,0)=p(0);
                MAT_ELEM(A,0,1)=p(1);
                MAT_ELEM(A,0,2)=p(2);
                MAT_ELEM(A,1,0)=p(3);
                MAT_ELEM(A,1,1)=p(4);
                MAT_ELEM(A,1,2)=p(5);
                MAT_ELEM(A,2,0)=0;
                MAT_ELEM(A,2,1)=0;
                MAT_ELEM(A,2,2)=1;
                transformationVec[particleIdx*frameNum+i]=p;
                cout<<"Transformation matrix is "<<endl<<A<<endl;
                // Apply the best transformation
                applyGeometry(BSPLINE3, transparticleImg, currentParticle, A, IS_NOT_INV, WRAP);
                particleAvg += transparticleImg;
                Iaux() = transparticleImg;
                Iaux.write(fnCurrentStackAl, ALL_IMAGES, true, WRITE_APPEND);
            }
            particleAvg /= double(frameNum);
            particleIdx++;
        }
    }
    // Write the polished particles
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
    Image<double> avgImages;
    produceSideInfo();
    show();
    extractAverageParticle();
    startTime = clock();
    extractParticels();
    printf("Time for extracting particles: %.2fs\n", (double)(clock() - startTime)/CLOCKS_PER_SEC);
    //computeAvgStack();
    //avgImages()=particleAvgStack;
    //avgImages.write("testavgimag.mrc");
    startTime = clock();
    particlePolishing3D();
    printf("Time for polishing particles: %.2fs\n", (double)(clock() - startTime)/CLOCKS_PER_SEC);

    //    avgImages() = particleAvgStack;
    //    avgImages.write("test2.stk",ALL_IMAGES,true,WRITE_APPEND);
    //computeRadialAvg();

    // Part for testing the program with phantom
    //    Image<double> II;
    //    MultidimArray<double> particleStack, avgImage, transparticleImg, particleAvg;
    //    Matrix1D<double> p, steps;
    //    std::vector< Matrix2D<double> > trMatrixVec;
    //    double cost=1e38;
    //    int iter;
    //
    //    p.initZeros(6);
    //    steps.initZeros(6);
    //    steps.initConstant(1.);
    //    p(0) = p(4) = 1.0;
    //
    //    A.initIdentity(3);
    //    // Reading the stack and the average
    //    II.read("ribosomeprojset30.stk");
    //    particleStack = II();
    //    II.read("refimage.xmp");
    //    avgParticle = II();
    //    transparticleImg.initZeros(particleSize, particleSize);
    //    currentParticle.initZeros(particleSize, particleSize);
    //    //currentParticle.initZeros(particleSize, particleSize);
    //    for (size_t it=0;it<4;it++)
    //    {
    //     particleAvg.initZeros(particleSize, particleSize);
    //        for (size_t i=0;i<30;i++)
    //        {
    //            //currentParticle.aliasImageInStack(particleStack, i);
    //         if (it>0)
    //         {
    //          p(0) = MAT_ELEM(trMatrixVec[i],0,0);
    //                p(1) = MAT_ELEM(trMatrixVec[i],0,1);
    //                p(2) = MAT_ELEM(trMatrixVec[i],0,2);
    //                p(3) = MAT_ELEM(trMatrixVec[i],1,0);
    //                p(4) = MAT_ELEM(trMatrixVec[i],1,1);
    //                p(5) = MAT_ELEM(trMatrixVec[i],1,2);
    //         }
    //            particleStack.getImage(i, currentParticle, 0);
    //            powellOptimizer(p, 1, 6, &L1costImage, this, 0.01, cost, iter, steps, false);
    //            MAT_ELEM(A,0,0)=p(0);
    //            MAT_ELEM(A,0,1)=p(1);
    //            MAT_ELEM(A,0,2)=p(2);
    //            MAT_ELEM(A,1,0)=p(3);
    //            MAT_ELEM(A,1,1)=p(4);
    //            MAT_ELEM(A,1,2)=p(5);
    //            MAT_ELEM(A,2,0)=0;
    //            MAT_ELEM(A,2,1)=0;
    //            MAT_ELEM(A,2,2)=1;
    //            trMatrixVec.push_back(A);
    //            cout<<A<<endl;
    //            applyGeometry(LINEAR, transparticleImg, currentParticle, A, IS_NOT_INV, DONT_WRAP, 0.);
    //            if (it==0)
    //            {
    //             II()=transparticleImg;
    //             II.write("testpowell.stk",ALL_IMAGES,true,WRITE_APPEND);
    //            }
    //            particleAvg += transparticleImg;
    //        }
    //        avgParticle = particleAvg / 30.0;
    //    }
}

