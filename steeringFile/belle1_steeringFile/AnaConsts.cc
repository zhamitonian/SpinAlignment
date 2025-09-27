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

double xmass[5]={0, 0.105658, 0.13957, 0.4937, 0.938272};
double PI = 3.141592653589793;

namespace cuts{
    double vertexZ = 3.0;
    double vertexR = 1.0;
    double minTrackPt = 0.1;
    double minPhotonE = 0.1;
    double minNGoodTracks = 3;
    double minPt = 0.1;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif