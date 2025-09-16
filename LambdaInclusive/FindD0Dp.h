#include "belle.h"
#include <cmath>
#include "particle/Particle.h"
#include "particle/utility.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "tables/evtcls.h"
#include "basf/module_descr.h"
#include "TFile.h"
#include "TH1F.h"
#include <TTree.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TRandom3.h>

#include <panther/panther.h>
#include "tables/evtcls.h"
//#include "toolbox/Thrust.h"
#include "toolbox/FuncPtr.h"
//#include "TMath.h"

#include BELLETDF_H
#include HEPEVT_H
#include MDST_H
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<float> Vfloat;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<Particle> VParticle;

class FindD0Dp{

public:
  FindD0Dp(void) {};
  ~FindD0Dp(void) {};

 //void find(Mdst_vee2_Manager::iterator, Vint ipion );
 int find(Vint ipion, Vint ikoan , HepPoint3D runIP);

  
  Particle getDmeson(int index){return m_Dmeson[index];} 
  int getCharge(int index){return m_charge[index];} 
  int getMode(int index){return m_Dmode[index];} 

  

private:

   VParticle m_Dmeson; 
   Vint m_Dmode; //0: D0: 0 ,1 ,2; Dpm: 3,4,5,6,7,8 
   Vint m_charge; //0: D0, +1 Dp, -1, Dm
   int  m_nDmeson;///how many candidates we found


};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
