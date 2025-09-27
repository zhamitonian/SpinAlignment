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
    extern double PI;

    namespace cuts
    {
        extern const double vertexZ;  //cuts on the vertex position
        extern const double vertexR;
        extern const double minTrackPt;
        extern const double minPhotonE;
        extern const int minNGoodTracks;
        extern const double minPt;
    }

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif