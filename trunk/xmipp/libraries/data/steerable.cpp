/***************************************************************************
 *
 * Authors:     Manuel Sanchez Pau 
 *              Carlos Oscar Sanchez Sorzano
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "steerable.h"
#include "fft.h"
#include "fftw.h"
#include "volume.h"
#include "histogram.h"
#include "filters.h"
#include "morphology.h"

// Remove wedge ------------------------------------------------------------
void MissingWedge::removeWedge(Matrix3D<double> &V) const
{
    Matrix2D<double> Epos, Eneg;
    Euler_angles2matrix(rotPos,tiltPos,0,Epos);
    Euler_angles2matrix(rotNeg,tiltNeg,0,Eneg);
    
    Matrix1D<double> freq(3), freqPos, freqNeg;
    Matrix1D<int> idx(3);

    XmippFftw transformer;
    Matrix3D< std::complex<double> > Vfft;
    transformer.FourierTransform(V,Vfft,false);

    FOR_ALL_ELEMENTS_IN_MATRIX3D(Vfft)
    {
        // Frequency in the coordinate system of the volume
        VECTOR_R3(idx,j,i,k);
        FFT_idx2digfreq(V,idx,freq);

        // Frequency in the coordinate system of the plane
        freqPos=Epos*freq;
        freqNeg=Eneg*freq;
        if (ZZ(freqPos)<0 || ZZ(freqNeg)>0)
            Vfft(k,i,j)=0;
    }
    transformer.inverseFourierTransform();
}

// Constructor -------------------------------------------------------------
Steerable::Steerable(double sigma, Matrix3D<double> &Vtomograph, 
    double deltaAng, const std::string &filterType, const MissingWedge *_MW) 
{
    MW=_MW;
    buildBasis(Vtomograph,sigma);
 
    // Choose a,b,c parameters as a function of the filterType
    double a; 
    double b; 
    double c;
    if (filterType == "wall") {
        // for wall structures
        a = -(1.0/4.0);
        b = 5.0/4.0;
        c = 5.0/2.0;
    }
    else{      
        // for filament structures
        a = 1.0;
        b = -(5.0/3.0);
        c = 10.0/3.0;
    }

    // Filter along tilt=0 and 180 => u=(1,0,0)
    double u0=1;
    double u1=0;
    double u2=0;                
    FOR_ALL_ELEMENTS_IN_MATRIX3D(Vtomograph){
        Vtomograph(k,i,j) = basis[0](k,i,j) * (a+b*u0*u0) +
                            basis[1](k,i,j) * (a+b*u1*u1) +
                            basis[2](k,i,j) * (a+b*u2*u2) +
                         c*(basis[3](k,i,j) * u0*u1 + 
                            basis[4](k,i,j) * u0*u2 + 
                            basis[5](k,i,j) * u1*u2);
    }

    // Filter the rest of directions and keep the maximum
    double Ntilt = round(180.0/deltaAng);  
    for (int i=1;i<Ntilt;i++){
        double tilt = deltaAng*i;
        double deltaRoti = deltaAng/SIND(tilt);
        double NrotP = round(360.0/deltaRoti);        
        for (int j=0;j<NrotP;j++){
            double rot = j*deltaRoti;
            double u0 = SIND(rot)*COSD(tilt);
            double u1 = SIND(rot)*SIND(tilt);
            double u2 = COSD(rot);
            FOR_ALL_ELEMENTS_IN_MATRIX3D(Vtomograph)
            {
                double filterval =
                    basis[0](k,i,j) * (a+b*u0*u0) +
                    basis[1](k,i,j) * (a+b*u1*u1) +
                    basis[2](k,i,j) * (a+b*u2*u2) +
                 c*(basis[3](k,i,j) * u0*u1 + 
                    basis[4](k,i,j) * u0*u2 + 
                    basis[5](k,i,j) * u1*u2);

                if(filterval>Vtomograph(k,i,j))
                    Vtomograph(k,i,j) = filterval;
            }
        }
    }
}

/* Build basis ------------------------------------------------------------- */
void Steerable::buildBasis(const Matrix3D<double> &Vtomograph, double sigma)
{
    std::vector< Matrix1D<double> > hx, hy, hz;
    generate1DFilters(sigma, Vtomograph, hx, hy, hz);
    for (int n=0; n<6; n++)
    {
        Matrix3D<double> aux;
        singleFilter(Vtomograph,hx[n],hy[n],hz[n],aux);
        basis.push_back(aux);
    }    
}        

void Steerable::singleFilter(const Matrix3D<double>& Vin,
    const Matrix1D<double> &hx, const Matrix1D<double> &hy, const Matrix1D<double> &hz,
    Matrix3D<double> &Vout){

    Matrix1D< std::complex<double> > H, Aux;
    Vout.initZeros(Vin);

    // Filter in X
    #define MINUS_ONE_POWER(n) (((n)%2==0)? 1:-1)
    FourierTransform(hx,H);    
    FOR_ALL_ELEMENTS_IN_MATRIX1D(H)
       if (i<(XSIZE(H)/2))
          H(i)*= MINUS_ONE_POWER(i);
       else
          H(i)*= MINUS_ONE_POWER(XSIZE(H)-i);

    Matrix1D<double> aux(XSIZE(Vin));    
    for (int k=0; k<ZSIZE(Vin); k++)
        for (int i=0; i<YSIZE(Vin); i++)
        {
            for (int j=0; j<XSIZE(Vin); j++)
                DIRECT_VEC_ELEM(aux,j)=DIRECT_VOL_ELEM(Vin,k,i,j);
            FourierTransform(aux,Aux);
            Aux*=H;
            InverseFourierTransform(Aux,aux);
            for (int j=0; j<XSIZE(Vin); j++)
                DIRECT_VOL_ELEM(Vout,k,i,j)=XSIZE(Aux)*DIRECT_VEC_ELEM(aux,j);
        }

    // Filter in Y
    FourierTransform(hy,H);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(H)
       if (i<(XSIZE(H)/2))
          H(i)*= MINUS_ONE_POWER(i);
       else
          H(i)*= MINUS_ONE_POWER(XSIZE(H)-i);

    aux.initZeros(YSIZE(Vin));
    for (int k=0; k<ZSIZE(Vin); k++)
        for (int j=0; j<XSIZE(Vin); j++)
        {
            for (int i=0; i<YSIZE(Vin); i++)
                DIRECT_VEC_ELEM(aux,i)=DIRECT_VOL_ELEM(Vout,k,i,j);
            FourierTransform(aux,Aux);
            Aux*=H;
            InverseFourierTransform(Aux,aux);
            for (int i=0; i<YSIZE(Vin); i++)
                DIRECT_VOL_ELEM(Vout,k,i,j)=XSIZE(Aux)*DIRECT_VEC_ELEM(aux,i);
        }

    // Filter in Z
    FourierTransform(hz,H);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(H)
       if (i<(XSIZE(H)/2))
          H(i)*= MINUS_ONE_POWER(i);
       else
          H(i)*= MINUS_ONE_POWER(XSIZE(H)-i);

    aux.initZeros(ZSIZE(Vin));    
    for (int i=0; i<YSIZE(Vin); i++)
        for (int j=0; j<XSIZE(Vin); j++)
        {
            for (int k=0; k<ZSIZE(Vin); k++)
                DIRECT_VEC_ELEM(aux,k)=DIRECT_VOL_ELEM(Vout,k,i,j);
            FourierTransform(aux,Aux);
            Aux*=H;
            InverseFourierTransform(Aux,aux);
            for (int k=0; k<ZSIZE(Vin); k++)
                DIRECT_VOL_ELEM(Vout,k,i,j)=XSIZE(Aux)*DIRECT_VEC_ELEM(aux,k);
        }
    
    // If Missing wedge
    if (MW!=NULL)
        MW->removeWedge(Vout);
}

/* Filter generation ------------------------------------------------------- */
void Steerable::generate1DFilters(double sigma,
    const Matrix3D<double> &Vtomograph,
    std::vector< Matrix1D<double> > &hx,
    std::vector< Matrix1D<double> > &hy,
    std::vector< Matrix1D<double> > &hz){

    // Initialization 
    Matrix1D<double> aux;    
    aux.initZeros(XSIZE(Vtomograph));
    aux.setXmippOrigin();
    for (int i=0; i<6; i++) hx.push_back(aux);
    
    aux.initZeros(YSIZE(Vtomograph));
    aux.setXmippOrigin();
    for (int i=0; i<6; i++) hy.push_back(aux);

    aux.initZeros(ZSIZE(Vtomograph));
    aux.setXmippOrigin();
    for (int i=0; i<6; i++) hz.push_back(aux);

    double sigma2=sigma*sigma;       
    double k1 =  1.0/pow((2.0*PI*sigma),(3.0/2.0));
    double k2 = -1.0/(sigma2);
    
    FOR_ALL_ELEMENTS_IN_MATRIX1D(hx[0])
    {        
        double i2=i*i;
        double g = -exp(-i2/(2.0*sigma2));
	hx[0](i) = k1*k2*g*(1.0-(i2/sigma2));
	hx[1](i) = k1*k2*g;
	hx[2](i) = k1*k2*g;
	hx[3](i) = k1*k2*k2*g*i;
	hx[4](i) = k1*k2*k2*g*i;
	hx[5](i) = k1*k2*k2*g;
    }    
    FOR_ALL_ELEMENTS_IN_MATRIX1D(hy[0])
    {
        double i2=i*i;
        double g = -exp(-i2/(2.0*sigma2));
        hy[0](i) = g;
        hy[1](i) = g*(1.0-(i2/sigma2));
        hy[2](i) = g;
        hy[3](i) = g*i;
        hy[4](i) = g;
        hy[5](i) = g*i;
    }
    FOR_ALL_ELEMENTS_IN_MATRIX1D(hz[0])
    {
        double i2=i*i;
        double g = -exp(-i2/(2.0*sigma2));
	hz[0](i) = g;
	hz[1](i) = g;
	hz[2](i) = g*(1.0-(i2/sigma2));
	hz[3](i) = g;
	hz[4](i) = g*i;
	hz[5](i) = g*i;
    }
}

void Steerable::generate3DFilter(Matrix3D<double>& h3D,
    std::vector< Matrix1D<double> > &hx,
    std::vector< Matrix1D<double> > &hy,
    std::vector< Matrix1D<double> > &hz)
{
    h3D.initZeros(XSIZE(hz[0]),XSIZE(hy[0]),XSIZE(hx[0]));
    h3D.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX3D(h3D)
        for (int n=0; n<6; n++)
            h3D(k,i,j)+=(hz[n](k)*hy[n](i)*hx[n](j));    
}

