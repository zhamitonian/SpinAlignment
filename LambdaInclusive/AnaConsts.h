#ifndef ANACONST_H
#define ANACONST_H

#include "CLHEP/Vector/LorentzVector.h"
#include "particle/Particle.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  //#if defined(BELLE_NAMESPACE)
  //using namespace Belle;
  //#endif

  extern double D0mass_lo;
  extern double D0mass_up;
  extern double Dpmass_lo;
  extern double Dpmass_up;
  extern double Ximass_lo;
  extern double Ximass_up;
  extern double lcmass_lo;
  extern double lcmass_up;
  extern double pi0mass_lo;
  extern double pi0mass_up;
  extern float m_pi;
  extern float m_k;
  extern float m_ks;
  extern float m_pr;
  extern float m_lmbda;
  extern float m_muon;
  extern float m_e;//doesn't matter
  extern float xmass[5];//doesn't matter
  extern int LambdaId;
  extern int anti_LambdaId;
  extern int leafParticleID[7];
  extern double m_inputPolarizationArray[5];
  extern double m_inputPolPtArray[5];
  extern double zbinRange[6];
  extern double ptbinRange[6];
  extern double ptbinRange_href[6];
  extern double m_inputCombinedPolArray[5][5];
  extern double m_inputCombinedPolArray_B[5][5];
  extern double m_inputZPtcombinedArray[5][5];
  extern double HadronZRange[6];
  extern double weight_zpt[5][5];

namespace cuts
{
  extern const float pi0SigLower;
  extern const float pi0SigUpper;
  extern const float pi0BGLower;
  extern const float pi0BGUpper;

  extern const float minZThrust; //min Z so that this particle goes into thrust computation
  extern const float minPtThrust;
  extern const float minZ;//min z so that this particle is used in asymmetry extraction
  extern const float minThrust;
  extern const float minGammaEBarrel;
  extern const float minGammaEEndcapFwd;
  extern const float minGammaEEndcapBkwd;
  extern const float minPi0GammaE;
  extern const float minThrustProj;//too big lets theta around thrust peak at around 1 and cos decay theta at zero
  extern const float maxThrustZ;///?????
  extern const float min2H_Z;
  extern const float minCosTheta; //barrell region?
  extern const float maxCosTheta;
  extern const float maxLabThrustCosTheta;
  extern const float minLabThrustCosTheta;

  extern const float barrelThetaMin;
  extern const float barrelThetaMax;

  extern const float minVisEnergy;//here only hadrons and gammas so far
  extern const float maxPi0GAsym;
  extern const float maxQt;
  extern const int minNTracks;//in charged pion prod, this doesn't make sense...
  extern const float vertexZ;  //cuts on the vertex position
  extern const float vertexR;
}
namespace kinematics
{

  extern int D0Tag;
  extern int DStarTag;
  //not ideal way to tag the fast pion
  extern float fastPionZ;

  extern int evtNr;
  extern int runNr;
  extern float E_miss;
  extern const double eler;//energies of l, h beam
  extern const double eher;
  extern const double theta;
  extern HepLorentzVector firstElectronCM;
  extern HepLorentzVector secondElectronCM;
  extern Hep3Vector CMBoost;
  extern HepLorentzVector cm;
  extern double Q;
  extern Hep3Vector thrustDirCM;
  extern Hep3Vector thrustDirLab;
  extern Hep3Vector thrustDirCMS_MC;
  extern Hep3Vector thrustDirLab_MC;
  extern Hep3Vector qqDirCMS_MC;
  extern Hep3Vector qqDirLab_MC;
  extern float thrustMag;
  extern double thrust;
  extern bool thrustZReverted;
  extern float thetaEThrust;
  extern Hep3Vector jet1;
  extern Hep3Vector jet2;
  extern Hep3Vector dijet[2];
  extern float jetE1;
  extern float jetE2;

  extern int jetNumParts1;
  extern int jetNumParts2;

  extern float jetFluffiness1;
  extern float jetFluffiness2;
  extern float R;

}
extern const double pi;

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
