#include "event/BelleEvent.h"
#include "particle/Particle.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
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

double xmass[5]={0, 0.105658, 0.13956995, 0.4937, 0.938272};
int   leafParticleID[7]={11,13,211,321,2212,22,2112};
double PI = 3.141592653589793;

namespace cuts{
    double vertexZ = 4.0;
    double vertexR = 2.0;
    double trkPt= 0.1;
    double trkP = 0.5;
    double trk_cosTheta_min = -0.511;
    double trk_cosTheta_max = 0.842;
    double Evis = 7;
    double PrimeVr = 1.5;
    double PrimeVz = 3.5;
    double HJM_min = 1.8;
    double HJM_Evis_ratio = 0.25;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif