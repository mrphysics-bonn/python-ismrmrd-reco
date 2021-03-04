#include "vdspiral.h"
#include <iomanip>
#include <stdio.h>
#include <fstream>

// copy in binary mode
bool copyFile(const char *SRC, const char* DEST) {
    std::ifstream src(SRC, std::ios::binary);
    std::ofstream dest(DEST, std::ios::binary);
    dest << src.rdbuf();
    return src && dest;
}

vdspiral::vdspiral(void) // Constructor
    : NonCartesianTraj()
    , m_eSpiralType(SpiralOut)
    , m_Nitlv(1)
    , m_dResolution(5.) //mm
{}


vdspiral::~vdspiral(void) // Destructor
{}


bool vdspiral::prep(int Nitlv, double res, std::vector<double> fov, std::vector<double> radius, double dMaxAmplitude, double dMinRiseTime, eSpiralType spiralType, double dLarmorConst, double dGradRasterTime) {
    
    m_dLarmorConst  = dLarmorConst;
    m_dResolution   = res;
    m_Nitlv         = Nitlv;
    m_fov           = fov;
    m_radius        = radius;
    m_dMaxAmplitude = dMaxAmplitude;
    m_dMinRiseTime  = dMinRiseTime;
    m_dGradRasterTime = dGradRasterTime;
    m_eSpiralType   = spiralType;
    
    if (fov.size()!=radius.size())
        return false;
    
    double kmax = 5./m_dResolution;  // kmax = 1/(2*res) BUT: kmax in 1/cm, res in mm
    double Gmax = m_dMaxAmplitude/10.;   // mT/m    -> G/cm
    double Smax = 100./m_dMinRiseTime;   // -> G/cm/ms
    double T = m_dGradRasterTime/1000.;  // us -> ms

    // unit transformation of fov and radius
    double fovmax = 0.;
    for (size_t k=0; k<fov.size();++k) {
        fov[k] /= 10.; // mm->cm
        if (fov[k]>fovmax)
            fovmax = fov[k];
        radius[k] *= kmax;
    }

    if (!vdspiral::vdSpiralDesign(m_Nitlv, fovmax, kmax, Gmax, Smax, fov, radius, spiralType, T))
        return false;
    
    #ifdef BUILD_SEQU
    if (!vdspiral::prepGradients())
        return false;
    #endif
    
    return true;
}


bool vdspiral::vdSpiralDesign(int Nitlv, double fovmax, double kmax, double Gmax, double Smax, std::vector<double> fov, std::vector<double> radius, eSpiralType spiralType, double T) {

    size_t k; // loop index

    //double dr = 1./500. * 1./(fovmax/Nitlv);
    double dr = 1./100. * 1./(fovmax/Nitlv); // a little faster
    long   nr = long(kmax/dr) + 1;
    
    std::vector<double> x, y, z;
    x.resize(nr, 0.);
    y.resize(nr, 0.);
    z.resize(nr, 0.);
    
    // calculate parametrized curve
    double theta = 0.;
    for (k=0; k<nr; k++) {
        double r = k*dr;
        double cFoV = fov.back();
        for (int l=0; l<fov.size(); ++l) {
            if (r < radius[l]) {
                if (l==0 || l==fov.size()-1)
                    cFoV = fov[l];
                else {// linearer übergang
                    double step = (r-radius[l-1])/(radius[l]-radius[l-1]);
                    cFoV = step * fov[l] + (1.-step)*fov[l-1];
                }
                break;
            }
        }        
        x[k] = r * cos(theta);
        y[k] = r * sin(theta);
        if (spiralType == DoubleSpiral) {
            theta += M_PI*dr*cFoV/Nitlv;
        } else {
            theta += 2.*M_PI*dr*cFoV/Nitlv;
        }
    }
            
    int n;
    double g0   = 0.; // to simplify sequence development, our gradient will start at 0.
    double gfin = 0.; // and end at 0.
    double *gx; double *gy; double *gz;
    //clock_t start = clock();
    minTimeGradientRIV(&x[0], &y[0], &z[0], nr, g0, gfin, Gmax, Smax, T, gx, gy, gz, n, -0.5, m_dLarmorConst/10.);
    //clock_t end = clock();
    //cout << "start = " << start << "  end = " << end << "  end-start = " << end-start << "  CLOCKS_PER_SEC = " << CLOCKS_PER_SEC << endl;

    // determine max gradient amplitudes
    m_dAx  = 0.;
    m_dAy  = 0.;
    m_dAmp = 0.;
    for (k=0; k<n; ++k) {
        if (fabs(gx[k])>m_dAx)
            m_dAx  = fabs(gx[k]);
        if (fabs(gy[k])>m_dAy)
            m_dAy  = fabs(gy[k]);
        double dAmp = sqrt(gx[k]*gx[k]+gy[k]*gy[k]);
        if (dAmp>m_dAmp)
            m_dAmp = dAmp;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////
    // calculate final gradient arrays
    ///////////////////////////////////////////////////////////////////////////////////
    m_vfGx.clear();
    m_vfGy.clear();
    m_vfGx.resize(n, 0.);
    m_vfGy.resize(n, 0.);
    // copy & scale gradient to interval -1...+1
    for (k=0; k<n; ++k) {
        m_vfGx[k] = (float) (gx[k] / (m_dAx>0?m_dAx:1.));
        m_vfGy[k] = (float) (gy[k] / (m_dAy>0?m_dAy:1.));
    }
    delete[] gx; delete[] gy; delete[] gz;
    
    // mirror & copy gradient
    if (spiralType == DoubleSpiral) {
        m_vfGx.resize(2*n, 0.);
        m_vfGy.resize(2*n, 0.);
        for (k=0; k<n; ++k) {
            m_vfGx[n+k] = m_vfGx[k];
            m_vfGy[n+k] = m_vfGy[k];
        }
        for (k=0; k<n/2; ++k) {
            std::swap(m_vfGx[k], m_vfGx[n-1-k]);
            std::swap(m_vfGy[k], m_vfGy[n-1-k]);
        }
        n *= 2;
    } else if (spiralType == ROI) {
        m_vfGx.resize(2*n+4, 0.);
        m_vfGy.resize(2*n+4, 0.);
        for (k=0; k<n; ++k) {
            m_vfGx[n+k+4] = -m_vfGx[k];
            m_vfGy[n+k+4] = -m_vfGy[k];
        }
        for (k=0; k<n/2; ++k) {
            std::swap(m_vfGx[n+k+4], m_vfGx[2*n+3-k]);
            std::swap(m_vfGy[n+k+4], m_vfGy[2*n+3-k]);
        }
        m_vfGx[n] = m_vfGx[n-1]/2.;
        m_vfGy[n] = m_vfGy[n-1]/2.;
        m_vfGx[n+1] = m_vfGx[n-1]/4.;
        m_vfGy[n+1] = m_vfGy[n-1]/4.;
        m_vfGx[n+2] = m_vfGx[n+4]/4.;
        m_vfGy[n+2] = m_vfGy[n+4]/4.;
        m_vfGx[n+3] = m_vfGx[n+4]/2.;
        m_vfGy[n+3] = m_vfGy[n+4]/2.;
        n = 2*n+2;
    } else if (spiralType == RIO) {
        m_vfGx.resize(2*n+2, 0.);
        m_vfGy.resize(2*n+2, 0.);
        for (k=0; k<n; ++k) {
            m_vfGx[n+k+2] = -m_vfGx[k];
            m_vfGy[n+k+2] = -m_vfGy[k];
        }
        for (k=0; k<n/2; ++k) {
            std::swap(m_vfGx[k], m_vfGx[n-1-k]);
            std::swap(m_vfGy[k], m_vfGy[n-1-k]);
        }
        m_vfGx[n] = m_vfGx[n-1]/2.;
        m_vfGy[n] = m_vfGy[n-1]/2.;
        m_vfGx[n+1] = m_vfGx[n+2]/2.;
        m_vfGy[n+1] = m_vfGy[n+2]/2.;
        n = 2*n+2;
    }

    // add points at start and end of gradient to avoid slewrate overflow
    m_vfGx.insert(m_vfGx.begin(), m_vfGx[0]/2.);
    m_vfGy.insert(m_vfGy.begin(), m_vfGy[0]/2.);
    m_vfGx.push_back(m_vfGx.back()/2.);
    m_vfGy.push_back(m_vfGy.back()/2.);
    n+=2;
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    // G/cm -> mT/m   
    m_dAx  *= 10.;
    m_dAy  *= 10.;
    m_dAmp *= 10.;

    // now calculate the gradient moments
    m_dMomX = 0.; m_dMomY = 0.; m_dMomZ = 0.;
    for (k=0; k < (int)m_vfGx.size(); ++k) {
        m_dMomX += m_dAx * m_vfGx[k] * m_dGradRasterTime;
        m_dMomY += m_dAy * m_vfGy[k] * m_dGradRasterTime;
    }
    
    m_dPreMomX  = 0.; m_dPreMomY = 0.; m_dPreMomZ = 0.;
    m_dPostMomX = 0.; m_dPostMomY = 0.; m_dPostMomZ = 0.;
    if (spiralType==SpiralIn) {
        m_dPreMomX = m_dMomX;
        m_dPreMomY = m_dMomY;
        // we have to time-reverse the trajectory!
        for (long k=0; k<(int)m_vfGx.size()/2; ++k) {
            std::swap(m_vfGx[k],m_vfGx[m_vfGx.size()-1-k]);
            std::swap(m_vfGy[k],m_vfGy[m_vfGy.size()-1-k]);
        }
    } else if (spiralType==SpiralOut) {
        m_dPostMomX = m_dMomX;
        m_dPostMomY = m_dMomY;
    } else if (spiralType==DoubleSpiral) {
        m_dPreMomX = m_dMomX/2.;
        m_dPreMomY = m_dMomY/2.;
        m_dPostMomX = m_dMomX/2.;
        m_dPostMomY = m_dMomY/2.;
    } else if (spiralType==RIO) {
        for (k=0; k < (int)m_vfGx.size()/2; ++k) {
            m_dPreMomX += m_dAx * m_vfGx[k] * m_dGradRasterTime;
            m_dPreMomY += m_dAy * m_vfGy[k] * m_dGradRasterTime;
            m_dPostMomX += m_dAx * m_vfGx[k+m_vfGx.size()/2] * m_dGradRasterTime;
            m_dPostMomY += m_dAy * m_vfGy[k+m_vfGx.size()/2] * m_dGradRasterTime;
        }
    }

    return true;
}


bool vdspiral::setSpiralType(eSpiralType spiralType) {
    bool bStatus = true;
    if (spiralType != m_eSpiralType) {
        // if ((m_eSpiralType != DoubleSpiral) && (spiralType != DoubleSpiral)) {
        //     // we have to time-reverse the trajectory!
        //     for (long k=0; k<(int)m_vfGx.size()/2; ++k) {
        //         std::swap(m_vfGx[k],m_vfGx[m_vfGx.size()-1-k]);
        //         std::swap(m_vfGy[k],m_vfGy[m_vfGy.size()-1-k]);
        //     }
        //     std::swap(m_dPreMomX, m_dPostMomX);
        //     std::swap(m_dPreMomY, m_dPostMomY);
        //     std::swap(m_dPreMomZ, m_dPostMomZ);
        //     #ifdef BUILD_SEQU
        //     bStatus = vdspiral::prepGradients();
        //     #endif
        //     return bStatus;
        // } else { // we need to recalculate the trajectory
            m_eSpiralType = spiralType;
            return this->prep(m_Nitlv, m_dResolution, m_fov, m_radius, m_dMaxAmplitude, m_dMinRiseTime, m_eSpiralType, m_dLarmorConst, m_dGradRasterTime);
        // }
    }
    return true;
}

bool vdspiral::calcTrajectory(std::vector<float> &vfKx, std::vector<float> &vfKy, std::vector<float> &vfDcf, long lADCSamples, int gridsize, double dADCshift, double dGradDelay) {
    //only spiral out supported for now!
    
    if (m_vfGx.size()==0 || m_vfGx.size()!=m_vfGy.size())
        return false;
    
    long lGradSamples = m_vfGx.size();
    long k,l;
    
    double dwelltime = (lGradSamples*m_dGradRasterTime + dADCshift)/lADCSamples;
    // iflag used in spline method to signal error 
    int *iflag;
    int iflagp;
    iflag = &iflagp;
    int *last;
    int lastp = 0;
    last = &lastp;
    
    // kgradx,kgrady: k-space trajectory on gradient raster
    int nFillpre  = 2;
    int nFillpost = 2;
    if (m_eSpiralType == SpiralOut)
        nFillpre  += int(dADCshift/m_dGradRasterTime);
    else
        nFillpost += int(dADCshift/m_dGradRasterTime);
    
    long lFilledSamples = lGradSamples+nFillpre+nFillpost;
    double *kgradx  = new double[lFilledSamples];
    double *kgrady  = new double[lFilledSamples];
    for (k=0;k<nFillpre+1;++k) {
        kgradx[k] = 0.;
        kgrady[k] = 0.;
    }
    double cumsumx=0.,cumsumy=0.;
    for (k=1;k<lGradSamples;++k) {
        cumsumx += m_dAx * (m_vfGx[k]+m_vfGx[k-1])/2.;
        kgradx[k+nFillpre]  = cumsumx * m_dGradRasterTime * m_dLarmorConst/1e5;
        cumsumy += m_dAy * (m_vfGy[k]+m_vfGy[k-1])/2.;
        kgrady[k+nFillpre]  = cumsumy * m_dGradRasterTime * m_dLarmorConst/1e5;
    }
    for (k=0;k<nFillpost;++k) {
        cumsumx += m_dAx * m_vfGx[lGradSamples-1]/2.;
        cumsumy += m_dAy * m_vfGy[lGradSamples-1]/2.;
        kgradx[k+nFillpre+lGradSamples] = cumsumx * m_dGradRasterTime * m_dLarmorConst/1e5;
        kgrady[k+nFillpre+lGradSamples] = cumsumy * m_dGradRasterTime * m_dLarmorConst/1e5;
    }
    
    double *coeff1x = new double[lFilledSamples];
    double *coeff2x = new double[lFilledSamples];
    double *coeff3x = new double[lFilledSamples];
    double *coeff1y = new double[lFilledSamples];
    double *coeff2y = new double[lFilledSamples];
    double *coeff3y = new double[lFilledSamples];
    double *t_grad  = new double[lFilledSamples];
    // Initalize gradient raster time
    for (k=0;k<lFilledSamples;++k)
        t_grad[k] = (k-nFillpre+1)*m_dGradRasterTime;
    
    spline(lFilledSamples, 0, 0, 0, 0, t_grad, kgradx, coeff1x, coeff2x, coeff3x, iflag);
    spline(lFilledSamples, 0, 0, 0, 0, t_grad, kgrady, coeff1y, coeff2y, coeff3y, iflag);

    // --------------------------------------------------------------
    // Interpolated curve
    // --------------------------------------------------------------
    vfKx.resize(lADCSamples*m_Nitlv,0.);
    vfKy.resize(lADCSamples*m_Nitlv,0.);
    for (k=0; k<lADCSamples; k++) {
        // Time for ACD sampling point
        double t_relativeToGrad = (k+0.5)*dwelltime + dGradDelay;
        if (m_eSpiralType == SpiralOut)
            t_relativeToGrad -= dADCshift;
        else
            t_relativeToGrad += dADCshift;
        vfKx[k] = (float) seval(lFilledSamples, t_relativeToGrad, t_grad, kgradx, coeff1x, coeff2x, coeff3x, last);
        vfKy[k] = (float) seval(lFilledSamples, t_relativeToGrad, t_grad, kgrady, coeff1y, coeff2y, coeff3y, last);
    }
    delete[] kgradx;
    delete[] kgrady;
    delete[] t_grad;
    delete[] coeff1x;
    delete[] coeff2x;
    delete[] coeff3x;
    delete[] coeff1y;
    delete[] coeff2y;
    delete[] coeff3y;
    
    if (m_eSpiralType == SpiralIn) {
        //for spiral in: make sure that trajectory ends in the center of k-space
        for (k=0; k<lADCSamples; k++) {
            vfKx[k] -= vfKx[lADCSamples-1];
            vfKy[k] -= vfKy[lADCSamples-1];
        }
    } else if (m_eSpiralType == DoubleSpiral || m_eSpiralType == RIO) {
        //for double spiral: make sure that middle of trajectory is in the center of k-space
        float midX = vfKx[lADCSamples/2];
        float midY = vfKy[lADCSamples/2];
        for (k=0; k<lADCSamples; k++) {
            vfKx[k] -= midX;
            vfKy[k] -= midY;
        }
    }
    
    //now calculate trajectory for other interleaves by rotation
    for (k=1;k<m_Nitlv;++k) {
        float phi = (float) (2.* M_PI * k / m_Nitlv);
        if (m_eSpiralType == DoubleSpiral && m_Nitlv%2==0) {
            // we only need to distribute the spirals over M_PI 
            phi /= 2.f;
        }
        for (l=0; l<lADCSamples; ++l) {
//             vfKx[l+k*lADCSamples] = (float) (cos(phi) * vfKx[l] + sin(phi) * vfKy[l]);
//             vfKy[l+k*lADCSamples] = (float) (-sin(phi) * vfKx[l] + cos(phi) * vfKy[l]);
            vfKx[l+k*lADCSamples] = (float) (cos(phi) * vfKx[l] - sin(phi) * vfKy[l]);
            vfKy[l+k*lADCSamples] = (float) (sin(phi) * vfKx[l] + cos(phi) * vfKy[l]);
        }
    }
    
    // Calculate density compensation function
    vfDcf = jacksonDCF(vfKx, vfKy, gridsize, 1.f);
    
    return true;
}


std::vector<float> vdspiral::jacksonDCF(std::vector<float> &vfKx, std::vector<float> &vfKy, int gridsize, float zeta) {
    
    int k,l,m;
    long nsamples = vfKx.size()/m_Nitlv;
    
    //scale zeta:
    //find maxk:
    float kmax = 0.f;
    for (k=0;k<nsamples;++k) {
        float tmp = vfKx[k]*vfKx[k]+vfKy[k]*vfKy[k];
        if (tmp>kmax)
            kmax = fabs(vfKx[k]);
    }
    kmax  = sqrt(kmax);
    zeta *= 2.f * kmax * 2.f/gridsize;
    
    //cut dcf at 0.85 of kmax
    for (k=0;k<nsamples;++k) {
        float tmp = sqrt(vfKx[k]*vfKx[k]+vfKy[k]*vfKy[k]);
        if (tmp>0.85f*kmax)
            break;
    }
    // wi is cutoff and normalized to average wi between cutoff_ix1 and cutoff_ix2
    int cutoff_ix1 = k;
    int cutoff_ix2 = cutoff_ix1 + (nsamples-cutoff_ix1)/4+1; 
    
    float zeta_sq = zeta*zeta;
    std::vector<float> wi(nsamples*m_Nitlv,1.);
    for (k=0;k<nsamples;++k) {
        float goal = 0.;
        float kxk = vfKx[k];
        float kyk = vfKy[k];
        // vorsicht: Skript nimmt gleichverteilte Interleaves an (Winkel)
        for (l=0; l<m_Nitlv;l++) {
            for (m=0;m<nsamples;++m) {
                float dx = vfKx[m+l*nsamples] - kxk;
                dx = dx*dx;
                if (dx < zeta_sq) {
                    float dr = vfKy[m+l*nsamples] - kyk;
                    dr = dx + dr*dr;
                    if (dr < zeta_sq) {
                        dr = (float) sqrt(dr);
                        //simple hann filter
                        float kern = 0.5f - 0.5f * (float)cos(2.*M_PI*(1.+dr/zeta)/2.);
                        goal += kern;
                    }
                }
            }
        }
        wi[k] = 1.f/goal;
    }
    
    // determine cutoff value for wi
    float cutoff_val = 0.;
    for (k=cutoff_ix1;k<cutoff_ix2;++k) 
        cutoff_val += wi[k];
    cutoff_val /= MAX(1,cutoff_ix2-cutoff_ix1);
    // normalize wi by cutoff value
    for (k=0;k<nsamples;++k)
        wi[k] = MIN(1.f,wi[k]/cutoff_val);
            
    //now copy wi from interleave 0 to other interleaves;
    for (k=1; k<m_Nitlv;++k) {
        for (l=0; l<nsamples; ++l)
            wi[l+k*nsamples] = wi[l];
    }
    return wi;
}


void vdspiral::saveTrajectory(long lADCSamples, int gridsize, double dADCshift, double dGradDelay) {
#ifdef WIN32
    std::vector<float> vfKx, vfKy, vfDcf;
    this->calcTrajectory(vfKx, vfKy, vfDcf, lADCSamples, gridsize, dADCshift, dGradDelay);

    std::ofstream kxfile, kyfile, wifile;
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );  
    strftime (buffer,80,"%Y%m%d_%H%M%S",timeinfo);

    std::string s;
    std::stringstream ss;
    ss << buffer;
    ss >> s;

    #if defined(__linux__)
    std::string pathname   = "/tmp/";
    #else
    std::string pathname   = "C:\\Temp\\";
    #endif
    
    std::string kxfilename = pathname + "spiral_kx_" + s + ".csv";
    std::string kyfilename = pathname + "spiral_ky_" + s + ".csv";
    std::string wifilename = pathname + "spiral_wi_" + s + ".csv";

    kxfile.open(kxfilename.c_str());
    kyfile.open(kyfilename.c_str());
    wifile.open(wifilename.c_str());
    for (size_t k=0; k<vfKx.size(); ++k) {
        if (k != 0) {
            kxfile << ", ";
            kyfile << ", ";
            wifile << ", ";
        }
        kxfile << vfKx[k];
        kyfile << vfKy[k];
        wifile << vfDcf[k];
    }
    kxfile.close();
    kyfile.close();
    wifile.close();
#endif
}


#ifdef BUILD_SEQU
void vdspiral::saveGradientShapes(sGRAD_PULSE* pGradPreX, sGRAD_PULSE* pGradPreY, sGRAD_PULSE* pGradPostX, sGRAD_PULSE* pGradPostY) { 
#ifdef WIN32
    if (m_vfGx.size()*m_vfGy.size()==0)
        return;

    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );  
    strftime (buffer,80,"%Y%m%d_%H%M%S",timeinfo);

    std::string s;
    std::stringstream ss;
    ss << buffer;
    ss >> s;

    #if defined(__linux__)
    std::string pathname   = "/tmp/";
    #else
    // std::string pathname   = "C:\\Temp\\";
    std::string pathname   = "gradshapes\\";
    #endif
    
    std::string filename_gr = pathname + "spiral_dephaser_Gr.csv";
    std::string filename_gp = pathname + "spiral_dephaser_Gp.csv";
    // try to delete dephaser gradient files (so that they only exists if written below)
    remove(filename_gr.c_str());
    remove(filename_gp.c_str());
    std::ofstream gpfile, grfile;
    if ((m_eSpiralType != SpiralOut) && (pGradPreX!=NULL) && (pGradPreY!=NULL)) {
        filename_gr = pathname + "spiral_dephaser_Gr.csv";
        filename_gp = pathname + "spiral_dephaser_Gp.csv";
        gpfile.open(filename_gr.c_str(), std::ofstream::out | std::ofstream::trunc);
        grfile.open(filename_gp.c_str(), std::ofstream::out | std::ofstream::trunc);
        long prephase_time = MAX(pGradPreX->getTotalTime(), pGradPreY->getTotalTime());
        for (size_t t = 0; t < prephase_time; t+=GRAD_RASTER_TIME) {
            // account for distant possibility that gx and gy gradients do not have same length
            long tx = t - prephase_time + pGradPreX->getTotalTime();
            long ty = t - prephase_time + pGradPreY->getTotalTime();
            if (tx != 0) {
                gpfile << ", ";
            }
            if (tx != 0) {
                grfile << ", ";
            }
            gpfile << std::setprecision(12) << pGradPreX->getCurrentAmplitude(tx);
            grfile << std::setprecision(12) << pGradPreY->getCurrentAmplitude(ty);
        }
        gpfile.close();
        grfile.close();

        copyFile(filename_gr.c_str(), (pathname + "spiral_dephaser_Gr_" + s + ".csv").c_str());
        copyFile(filename_gp.c_str(), (pathname + "spiral_dephaser_Gp_" + s + ".csv").c_str());
    }
    
    // if ((m_eSpiralType != SpiralIn) && (pGradPostX!=NULL) && (pGradPostY!=NULL)) {
    //     filename = pathname + "spiral_rephaser_read_" + s + ".csv";
    //     gpfile.open(filename.c_str());
    //     filename = pathname + "spiral_rephaser_phase_" + s + ".csv";
    //     grfile.open(filename.c_str());
    //     long prephase_time = MAX(pGradPostX->getTotalTime(), pGradPostY->getTotalTime());
    //     for (size_t t = 0; t < prephase_time; t+=GRAD_RASTER_TIME) {
    //         // account for distant possibility that gx and gy gradients do not have same length
    //         long tx = t - prephase_time + pGradPostX->getTotalTime();
    //         long ty = t - prephase_time + pGradPostY->getTotalTime();
    //         if (tx != 0) {
    //             gpfile << ", ";
    //         }
    //         if (tx != 0) {
    //             grfile << ", ";
    //         }
    //         gpfile << std::setprecision(12) << pGradPostX->getCurrentAmplitude(tx);
    //         grfile << std::setprecision(12) << pGradPostY->getCurrentAmplitude(ty);
    //     }
    //     gpfile.close();
    //     grfile.close();
    // }

    filename_gr = pathname + "spiral_Gr.csv";
    filename_gp = pathname + "spiral_Gp.csv";
    gpfile.open(filename_gr.c_str(), std::ofstream::out | std::ofstream::trunc);
    grfile.open(filename_gp.c_str(), std::ofstream::out | std::ofstream::trunc);
    for (size_t t=0; t<m_GSpiralX.getTotalTime(); t+=GRAD_RASTER_TIME) {
        if (t != 0) {
            gpfile << ", ";
            grfile << ", ";
        }
        gpfile << std::setprecision(12) << m_GSpiralX.getCurrentAmplitude(t);
        grfile << std::setprecision(12) << m_GSpiralY.getCurrentAmplitude(t);
    }
    gpfile.close();
    grfile.close();
    copyFile(filename_gr.c_str(), (pathname + "spiral_Gr_" + s + ".csv").c_str());
    copyFile(filename_gp.c_str(), (pathname + "spiral_Gp_" + s + ".csv").c_str());
#endif
}

bool vdspiral::prepGradients() {
    m_GSpiralX.setRampShape(&m_vfGx[0], m_vfGx.size(), 0, true);
    m_GSpiralX.set(GRAD_RASTER_TIME*m_vfGx.size(), GRAD_RASTER_TIME*m_vfGx.size(), 0, m_dAx);
    if (!m_GSpiralX.prep()) {
        return false;
    }
    
    m_GSpiralY.setRampShape(&m_vfGy[0], m_vfGy.size(), 0, true);
    m_GSpiralY.set(GRAD_RASTER_TIME*m_vfGy.size(), GRAD_RASTER_TIME*m_vfGy.size(), 0, m_dAy);
    if (!m_GSpiralY.prep()) {
        return false;
    }
    return true;
}
#endif
