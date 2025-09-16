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

class FindLambdaC{

public:
  FindLambdaC(void) {};
  ~FindLambdaC(void) {};

 //void find(Mdst_vee2_Manager::iterator, Vint ipion );
 int find(HepLorentzVector pLambda, Vint ipion );

  
  Particle getLambdac(int index){return m_lambdac[index];} 
  int getCharge(int index){return m_charge[index];} 
  int getMode(int index){return m_lcmode[index];} 

  

private:

   VParticle m_lambdac; 
   Vint m_lcmode; //0: lambda + pi ; 1 lambda + pi + pi0 ; 2 lambda + 3pi
   Vint m_charge; 
   int  m_nLambdac;///how many candidates we found


};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
