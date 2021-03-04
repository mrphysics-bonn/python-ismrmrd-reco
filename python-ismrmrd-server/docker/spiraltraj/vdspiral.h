#ifndef VDSPIRAL
#define VDSPIRAL 1

#include "nonCartesianTraj.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>

#ifdef BUILD_SEQU
    #include "MrServers/MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"
    #include "MrServers/MrMeasSrv/SeqIF/libRT/libRT.h"
#endif

#if defined (BUILD_SEQU) || defined (BUILD_ICE)
    #include "MrServers/MrMeasSrv/SeqIF/libRT/libRTDefines.h"
#else
    #define GRAD_RASTER_TIME 10
#endif

#ifndef MAX
    #define MAX(a,b)  ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
    #define MIN(a,b)    ((a)>(b)?(b):(a))
#endif

class vdspiral: public NonCartesianTraj {
public:
    enum eSpiralType {
        SpiralOut = 1,
        SpiralIn = 2,
        DoubleSpiral = 3,
        ROI = 4,
        RIO = 5
    };

    vdspiral  (void);            // Constructor
    ~vdspiral (void);            // Destructor

    bool prep(int Nitlv, double res, std::vector<double> fov, std::vector<double> radius, double dMaxAmplitude, double dMinRiseTime, eSpiralType spiralType = SpiralOut, double dLarmorConst = 42.5756, double dGradRasterTime = GRAD_RASTER_TIME);
    
    inline long getlPoints ()  {return m_vfGx.size();}
    inline long getTotalTime() {return long(m_dGradRasterTime * getlPoints());}
    
    bool setSpiralType(eSpiralType spiralType = SpiralOut);
    virtual bool calcTrajectory(std::vector<float> &vfKx, std::vector<float> &vfKy, std::vector<float> &vfDcf, long lADCSamples, int gridsize=128, double dADCshift=0., double dGradDelay=0.);
    
    void saveTrajectory(long lADCSamples, int gridsize=128, double dADCshift=0., double dGradDelay=0.);
    #ifdef BUILD_SEQU
    void saveGradientShapes(sGRAD_PULSE* pGradPreX=NULL, sGRAD_PULSE* pGradPreY=NULL, sGRAD_PULSE* pGradPostX=NULL, sGRAD_PULSE* pGradPostY=NULL);
    inline void run (long lT) {
        // fRTEI(lT, 0,  0, 0, &m_GSpiralY, &m_GSpiralX, 0, 0);
        fRTEI(lT, 0,  0, 0, &m_GSpiralX, &m_GSpiralY, 0, 0);
    }
    inline sGRAD_PULSE_ARB* getPointerToX() {return &(m_GSpiralX);}
    inline sGRAD_PULSE_ARB* getPointerToY() {return &(m_GSpiralY);}
    #endif

    std::vector<float>& getGradX() { return m_vfGx; }
    std::vector<float>& getGradY() { return m_vfGy; }
protected:
    eSpiralType m_eSpiralType;

    double m_dResolution, m_dMaxAmplitude, m_dMinRiseTime;
    int m_Nitlv;
    std::vector<double> m_fov, m_radius;

    bool vdSpiralDesign(int Nitlv, double fovmax, double kmax, double Gmax, double Smax, std::vector<double> fov, std::vector<double> radius, eSpiralType spiralType = SpiralOut, double T = GRAD_RASTER_TIME/1000.);

    std::vector<float> jacksonDCF(std::vector<float> &vfKx, std::vector<float> &vfKy, int gridsize=128, float zeta=1.);
    
    #ifdef BUILD_SEQU
    sGRAD_PULSE_ARB m_GSpiralX, m_GSpiralY;
    bool prepGradients(void);
    #endif
};

#endif
