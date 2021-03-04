#include "nonCartesianTraj.h"

NonCartesianTraj::NonCartesianTraj(void) // Constructor
    : m_dAx(0.)
    , m_dAy(0.)
    , m_dAz(0.)
    , m_dAmp(0.)
    , m_dMomX(0.)
    , m_dMomY(0.)
    , m_dMomZ(0.)
    , m_dPreMomX(0.)
    , m_dPreMomY(0.)
    , m_dPreMomZ(0.)
    , m_dPostMomX(0.)
    , m_dPostMomY(0.)
    , m_dPostMomZ(0.)
    , m_dLarmorConst(0.)
    , m_dGradRasterTime(0)
{}

NonCartesianTraj::~NonCartesianTraj(void) { // Destructor
}
