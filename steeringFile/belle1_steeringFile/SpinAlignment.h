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
#include "toolbox/FuncPtr.h"

#include BELLETDF_H
#include HEPEVT_H
#include MDST_H

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<Hep3Vector> Vp3;
typedef std::vector<Particle> VParticle;

class SpinAlignment : public Module{
    public:
        SpinAlignment();
        ~SpinAlignment() {}
        void init(int *dummy);
        void begin_run(BelleEvent* evptr, int* status);
        void end_run(BelleEvent* evptr, int* status);
        void hist_def();
        void event(BelleEvent* evpter, int* status);
        void disp_stat(const char*);
        void term();
        void getDrDz(Mdst_charged_Manager::iterator chr_it, int masshyp, double& dr, double& dz, double& refitPx, double& refitPy, double& refitPz);
        void other(int*, BelleEvent*, int*);

    public: // BASF parameter
        char* output_filename;

    private:
        int countEvt;
        TRandom3 r1;
        TRandom3 r5;
        TRandom3 r6;
        TRandom3 r7;
        TRandom3 r8;
        TRandom3 r9;
        double PI;
        TFile * output_file;
        TTree * tree;
        double ler_e;
        double her_e;

        struct var_collection{
            int evtNo;
            int runNo;
            int expNo;
            double Q;
            double visEnergy;
            double Evis;
            int nGood;
            int nPip;
            int nPim;
            int nCluster;
            double thrust;
            double thrust_theta;
            double thrust_phi;
            double z;
            double pt;
            double cms_vecP[4];
        } m_info;

        int nPhoton;
        int nCluster;
        double photon_p[50];
        double photon_theta[50];
        double photon_phi[50];

        int nPip;
        int nPim;
        double pip_vecP[50][4];
        double pim_vecP[50][4];
        double pip_theta[50];
        double pim_theta[50];
        double pip_phi[50];
        double pim_phi[50];
};


extern "C" Module_descr *mdcl_SpinAlignment() {
    SpinAlignment *module = new SpinAlignment;
    Module_descr *dscr = new Module_descr( "SpinAlignment", module );
    dscr->define_param( "output_filename", "output root file name",255 ,(char*)module->output_filename );
    return dscr;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif