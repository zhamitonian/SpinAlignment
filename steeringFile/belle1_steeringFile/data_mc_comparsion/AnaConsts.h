#ifndef ANACONST_H
#define ANACONST_H

#include "CLHEP/Vector/LorentzVector.h"
#include "particle/Particle.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

    extern double xmass[5];
    extern int leafParticleID[7];
    extern double PI;

    namespace cuts
    {
        extern const double vertexZ;  //cuts on the vertex position
        extern const double vertexR;
        extern const double trkPt;
        extern const double trkP;
        extern const double trk_cosTheta_min;
        extern const double trk_cosTheta_max;
        extern const double Evis;
        extern const double PrimeVr;
        extern const double PrimeVz;
        extern const double HJM_min;
        extern const double HJM_Evis_ratio;
    }

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif