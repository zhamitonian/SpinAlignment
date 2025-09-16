#include "event/BelleEvent.h"
#include "particle/Particle.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
//#include "TF1.h"
//#include "TLorentzVector.h"
//#include "TVector3.h"
//#include "TFile.h"
//#include "TH1F.h"
//#include "LorentzVector.h" //clhep
#include "belle.h"
#include <kid/atc_pid.h>
#include <eid/eid.h>
#include <mdst/Muid_mdst.h>
#include "math.h"


#include "CLHEP/Vector/LorentzVector.h"
#include <vector>
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


  //#if defined(BELLE_NAMESPACE)
  //using namespace Belle;
  //#endif
  double D0mass_lo=1.82;
  double D0mass_up=1.90;
  double Dpmass_lo=1.82;
  double Dpmass_up=1.90;
  double Dsmass_lo=1.90;
  double Dsmass_up=2.02;
  double Ximass_lo =1.29;
  double Ximass_up =1.35;
  double lcmass_lo=2.23;
  double lcmass_up=2.34;
  double pi0mass_lo=0.10;
  double pi0mass_up=0.16;
  float m_pi=0.13957;
  float m_k=0.4937;
  float m_ks=0.4977;
  float m_pr=0.938272;
  float m_lmbda=1.1156;
  float m_muon=0.105658;
  float m_e=0;//doesn't matter
  float xmass[5]={0,0.105658,0.13957,0.4937,0.938272};
  int   LambdaId=3122;
  int   anti_LambdaId=-3122;
  //iteration first time
  //double weight_zpt[5][5]={{0.812795, 0.911199, 1.12967, 1.78332, 2.78214},{ 0.932453, 1.04155, 1.15736, 1.42422, 2.24349}, {0.964413, 1.07055, 1.18032, 1.2892, 1.56963}, {0.861968, 0.962688, 1.0955, 1.24486, 1.45266},{0.861968, 0.962688, 1.0955, 1.24486, 1.45266}};
  //iteration second time
  //double weight_zpt[5][5]= {{0.793364, 0.909597, 1.17436, 1.98032, 3.31144}, {0.942514, 1.08807, 1.13048, 1.25215, 2.66497}, {0.962147, 1.10134, 1.21964, 1.28562, 1.82804}, {0.78636, 0.896728, 1.04836, 1.26317, 1.67445}, {0.78636, 0.896728, 1.04836, 1.26317, 1.67445}};
 double weight_zpt[5][5]= {  
    {0.771298, 0.847986, 1.03685, 1.59515, 2.46909},
    {0.866569, 0.983894, 1.05231, 1.15568, 2.03338},
    {0.877056, 0.999202, 1.09419, 1.07028, 1.49217},
    {0.756853, 0.836067, 0.968921, 1.0377, 1.34338},
    {0.756853, 0.836067, 0.968921, 1.0377, 1.34338}};
  //int   leafParticleID[5]={11,13,211,321,2212,12,14,22,2112}
  int   leafParticleID[7]={11,13,211,321,2212,22,2112};
  double m_inputPolarizationArray[5]={-0.005,-0.01, -0.03, -0.040, -0.050}; //in z bins
  double m_inputPolPtArray[5] = {-0.01, -0.03, -0.025, -0.02, -0.02};
  double zbinRange[6]={0.0,0.3, 0.4, 0.5, 0.7, 100};
  double ptbinRange[6]={0.0,0.3,0.5,0.8,1.0, 100};//thrust frame
  double ptbinRange_href[6]={0.0,0.5,1.0,1.5,2.5,100.0};//thrust frame
  //double m_inputZPtcombinedArray[5][5] = {{-0.01, -0.015,-0.01,-0.01,-0.01,},{-0.01,-0.03,-0.02,-0.02,-0.02},{-0.03,-0.05,-0.04,-0.03,-0.03},{-0.04,-0.05,-0.07,-0.06,-0.06}, {-0.04,-0.05,-0.07,-0.06,-0.06}};
  //reconstructed  from data + smearing correction 
  //double m_inputZPtcombinedArray[5][5] = {{-0.00125074 , -0.00963629 , -0.00934841 , -0.00785457 , -0.0499679},{-0.0245057 , -0.0137538 , -0.00410763 , -0.00978075 , -0.00834879},{-0.0519334 , -0.0497287 , -0.0420356 , -0.0319212 , -0.0323314},{-0.0447815 , -0.059117 , -0.0663114 , -0.0600209 , -0.0802691}, {-0.0447815 , -0.059117 , -0.0663114 , -0.0600209 , -0.0802691}};
  double m_inputZPtcombinedArray[5][5] = {
   { 0.000529477 , -0.0107711 , -0.00937779 , -0.00619273 , -0.0600498},
   {-0.0291174 , -0.0136605 , 0.000325202 , -0.0138386 , -0.00706579},
   {-0.0531582 , -0.053586 , -0.0410579 , -0.0229692 , -0.0334279},
   //{-0.0307192 , -0.0661024 , -0.0775013 , -0.023404 , -0.118468},
   //{-0.0307192 , -0.0661024 , -0.0775013 , -0.023404 , -0.118468}
   {-0.0307192 , -0.0661024 , -0.0775013 , -0.09, -0.118468},
   {-0.0307192 , -0.0661024 , -0.0775013 , -0.09, -0.118468}
   };
  double m_inputZPtcombinedArray_href[5][5] = {{-0.01, -0.015,-0.01,-0.01,-0.01},{-0.01,-0.03,-0.02,-0.02,-0.02},{-0.03,-0.05,-0.04,-0.03,-0.03},{-0.04,-0.05,-0.07,-0.06,-0.06},{-0.04,-0.05,-0.07,-0.06,-0.06}};
  //double m_inputCombinedPolArray[5][5]={{-0.03, -0.05, -0.07, -0.08, -0.08},{-0.03, -0.05, -0.07, -0.08, -0.08}, {-0.03, -0.05, -0.07, -0.08, -0.08}, {-0.03, -0.05, -0.07, -0.08, -0.08}, {-0.03, -0.05, -0.07, -0.08, -0.08}};
  //double m_inputCombinedPolArray_B[5][5]= {{0.02, 0.01, 0.005, 0.005, 0.005}, {0.00,  0.00, 0.00, 0.00, 0.00}, {-0.005, -0.005, -0.005, -0.01, -0.01}, {-0.02, -0.02, -0.03, -0.05, -0.05}, {-0.02, -0.02, -0.03, -0.05, -0.05}};
  double HadronZRange[6]={0.0,0.3, 0.4, 0.5, 0.7, 1.0};

namespace cuts
{
  //pi0 cuts for signal and background
  float pi0SigLower=0.12;
  float pi0SigUpper=0.15;
  float pi0BGLower=0.21;
  float pi0BGUpper=0.3;


  //0.02 should corresnpond to 0.1 gev
  //  float minZThrust=0.02; //min Z so that this particle goes into thrust computation
    float minZThrust=0.0; //min Z so that this particle goes into thrust computation
    float minZ=0.1;//min z so that this particle is used in asymmetry extraction
  //float minZ=0.0;//min z so that this particle is used in asymmetry extraction
    float minThrust=0.8;
  //    float minThrust=0.0;

  //this is for the barrel, 
  float minGammaEBarrel=0.05;
  float minGammaEEndcapFwd=0.075;
  float minGammaEEndcapBkwd=0.1;
  float minPi0GammaE=0.05;
  //  float minThrustProj=0.1;//too big lets theta around thrust peak at around 1 and cos decay theta at zero
  // float minThrustProj=0.8;//too big lets theta around thrust peak at around 1 and cos decay theta at zero
     float minThrustProj=0.4;//too big lets theta around thrust peak at around 1 and cos decay theta at zero

    float maxThrustZ=0.75;//make sure that thrust doesn't point in the direction of the endcaps
  //  float maxThrustZ=0.878; //this should be symmetric in CMS (obviously) and at least 0.3 in eta away from cdc border--> this is for +-1.37, so +-0.1
  //  float maxThrustZ=0.878; //this should be symmetric in CMS (obviously) and at least 0.3 in eta away from cdc border
  //  float maxThrustZ=0.79; //this is for eta+-1.07, allowing +- 0.4 on both sides
  //  float maxThrustZ=0.88; //this is for eta+-1.37, just 0.1 away from the endcaps, just used for studies how the resolution looks like. Can always cut on later...
  //   float maxThrustZ=1.0;//don't care about the endcaps
  float min2H_Z=0.2;
  //float min2H_Z=0.0;
  //    float minCosTheta=-0.6; //barrell region?
  //    float maxCosTheta=0.9;  //barrel region
  //if we use jets we should leave this open
  float minCosTheta=-100;
  float maxCosTheta=2.0;
  //these angles correspond to eta of -1.32 to 1.9 in the lab system. Since eta is linear under boosts we can use that to have reasonable cuts on the thrust theta
  //if we assume a jet cone of 2*0.3 (very generous, half of that should be sufficient): 23 deg / 140 deg
  // 2* 0.2 ( so 0.2 from the edge: 21/143  
  //  float minCosTheta=-0.86; //wih endcaps (17 - 150 degrees, this is the CDC, the EMC goes from 12.4 to 155.1, central is 32.2 - 128.7 //endcaps 12.4-31.4, 130.7 - 155.1
  //  float maxCosTheta=0.95;  //

  //used to determine which energy cut to use for photons
  float barrelThetaMin=32.2;
  float barrelThetaMax=128.7;

  //from the above 0.2 jet cone, does this make sense in addition to z projection?, The values here are 0.936-->20 deg, -0.8==>143.13 deg not even symmetric?
  //change to accept all: 1.0, -1.0
  float maxLabThrustCosTheta=1.0;
  float minLabThrustCosTheta=-1.0;

  float minVisEnergy=7.5;//here only hadrons and gammas so far
  float maxPi0GAsym=0.8;//looked at distribution
  //  float maxQt=3.5;
  float maxQt=1000.0;
  int minNTracks=3;
  float vertexZ=4.0;
  float vertexR=2.0;
  float minPtThrust=0.0;
// min tranverse momentum so that it is used in thrust computation, put to 0.0 for compatibility with Ami
}
namespace kinematics
{
  int evtNr;
  int runNr;

  int D0Tag;
  int DStarTag;
  float fastPionZ;

  float E_miss;
  double eler=3.499218;//energies of l, h beam
  double eher(7.998213);
  double theta(0.022);
  HepLorentzVector firstElectronCM;
  HepLorentzVector secondElectronCM;
  Hep3Vector CMBoost;
  HepLorentzVector cm;
  double Q;
  bool thrustZReverted;
  Hep3Vector thrustDirCM;
  Hep3Vector thrustDirLab;
  Hep3Vector thrustDirCMS_MC;
  Hep3Vector thrustDirLab_MC;
  Hep3Vector qqDirCMS_MC;
  Hep3Vector qqDirLab_MC;
  float thrustMag;
  double thrust;
  float thetaEThrust;

  Hep3Vector jet1;
  Hep3Vector jet2;
  Hep3Vector dijet[2];


 float jetE1;
 float jetE2;

 int jetNumParts1;
 int jetNumParts2;

 float jetFluffiness1;
 float jetFluffiness2;
  float R;

}
double pi=3.14159265;
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
