#ifndef ANACONST_H
#define ANACONST_H

#include "CLHEP/Vector/LorentzVector.h"
#include "particle/Particle.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

    extern float xmass[5];

    namespace cuts
    {
        extern const double vertexZ;  //cuts on the vertex position
        extern const double vertexR;
        extern const double minTrackPt;
        extern const double minPhotonE;
        extern const int minNGoodTracks;
        extern const double minPt;
    }

    namespace kinematics
    {
        extern double ler_e;
        extern double her_e;
        extern HepLorentzVector firstElectronCM;
        extern HepLorentzVector secondElectronCM;
        extern Hep3Vector CMBoost;
        extern HepLorentzVector cm;
        extern double Q;
        extern double thrust;
        extern double thrust_theta;
        extern double thrust_phi;
    }

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif