#ifndef NONCARTESIANTRAJ_H
#define NONCARTESIANTRAJ_H 1

#ifdef WIN32
#pragma warning(disable: 4514 4710)
#endif

#include "mtg_functions.h"
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class NonCartesianTraj { 
public:
    NonCartesianTraj(void);
    ~NonCartesianTraj(void);

    inline double getMaxAmplitudeX()   {return m_dAx;}
    inline double getMaxAmplitudeY()   {return m_dAy;}
    inline double getMaxAmplitudeZ()   {return m_dAz;}
    inline double getMaxAbsAmplitude() {return m_dAmp;}
    inline double getMomentumX()       {return m_dMomX;}
    inline double getMomentumY()       {return m_dMomY;}
    inline double getMomentumZ()       {return m_dMomZ;}
    inline double getPreMomentumX()    {return m_dPreMomX;}
    inline double getPreMomentumY()    {return m_dPreMomY;}
    inline double getPreMomentumZ()    {return m_dPreMomZ;}
    inline double getPostMomentumX()   {return m_dPostMomX;}
    inline double getPostMomentumY()   {return m_dPostMomY;}
    inline double getPostMomentumZ()   {return m_dPostMomZ;}
    inline double getMaxAbsMomentum()  {return sqrt(m_dMomX*m_dMomX + m_dMomY*m_dMomY + m_dMomZ*m_dMomZ);}
    inline double getMaxAbsPreMomentum()   {return sqrt(m_dPreMomX*m_dPreMomX + m_dPreMomY*m_dPreMomY + m_dPreMomZ*m_dPreMomZ);}
    inline double getMaxAbsPostMomentum()  {return sqrt(m_dPostMomX*m_dPostMomX + m_dPostMomY*m_dPostMomY + m_dPostMomZ*m_dPostMomZ);}
    
    inline double getMomentumX(double phi, double theta=0.) {return (cos(theta)*cos(phi)*m_dMomX - sin(phi)*m_dMomY + sin(theta)*cos(phi)*m_dMomZ);}
    inline double getMomentumY(double phi, double theta=0.) {return (cos(theta)*sin(phi)*m_dMomX + cos(phi)*m_dMomY + sin(theta)*sin(phi)*m_dMomZ);}
    inline double getMomentumZ(double phi, double theta=0.) {(void) phi; return (-sin(theta)*m_dMomX + cos(theta)*m_dMomZ);}
    
    inline double getPreMomentumX(double phi, double theta=0.) {return (cos(theta)*cos(phi)*m_dPreMomX - sin(phi)*m_dPreMomY + sin(theta)*cos(phi)*m_dPreMomZ);}
    inline double getPreMomentumY(double phi, double theta=0.) {return (cos(theta)*sin(phi)*m_dPreMomX + cos(phi)*m_dPreMomY + sin(theta)*sin(phi)*m_dPreMomZ);}
    inline double getPreMomentumZ(double phi, double theta=0.) {(void) phi; return (-sin(theta)*m_dPreMomX + cos(theta)*m_dPreMomZ);}
    
    inline double getPostMomentumX(double phi, double theta=0.) {return (cos(theta)*cos(phi)*m_dPostMomX - sin(phi)*m_dPostMomY + sin(theta)*cos(phi)*m_dPostMomZ);}
    inline double getPostMomentumY(double phi, double theta=0.) {return (cos(theta)*sin(phi)*m_dPostMomX + cos(phi)*m_dPostMomY + sin(theta)*sin(phi)*m_dPostMomZ);}
    inline double getPostMomentumZ(double phi, double theta=0.) {(void) phi; return (-sin(theta)*m_dPostMomX + cos(theta)*m_dPostMomZ);}
    
protected:
    std::vector<float> m_vfGx, m_vfGy, m_vfGz;
    double m_dAx, m_dAy, m_dAz, m_dAmp, m_dMomX, m_dMomY, m_dMomZ, m_dPreMomX, m_dPreMomY, m_dPreMomZ, m_dPostMomX, m_dPostMomY, m_dPostMomZ;
    double m_dGradRasterTime;
    double m_dLarmorConst;
};

#endif
