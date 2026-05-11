#include "belle.h"
#include <cmath>
#include "particle/Particle.h"
#include "particle/utility.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "tables/evtcls.h"
#include "basf/module_descr.h"

#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TRandom3.h>

#include <panther/panther.h>
#include "tables/evtcls.h"
#include "toolbox/FuncPtr.h"

#include BELLETDF_H
#include MDST_H

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

typedef std::vector<int> Vint;
typedef std::vector<bool> Vbool;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<Hep3Vector> Vp3;
typedef std::vector<Particle> VParticle;

class KsSpinAlignment_NullTest : public Module{
    public:
        KsSpinAlignment_NullTest();
        ~KsSpinAlignment_NullTest() {}
        void init(int *dummy);
        void begin_run(BelleEvent* evptr, int* status);
        void end_run(BelleEvent* evptr, int* status);
        void hist_def();
        void event(BelleEvent* evpter, int* status);
        void disp_stat(const char*);
        void term();
        void getDrDz(Mdst_charged_Manager::iterator chr_it, int masshyp, double& dr, double& dz, double& refitPx, double& refitPy, double& refitPz);
        void other(int*, BelleEvent*, int*);
        std::pair<double, double> calculateHeavyJetMassEnergy(const Vp4 &particles, const Hep3Vector &thrustAxis);
        std::pair<double, double> calculateSphericityAplanarity(const std::vector<HepLorentzVector>& particles);

    public: // BASF parameter
        char* output_filename;

    private:
        int countEvt;
        int countEvt_fail;
        TRandom3 r1;
        TRandom3 r5;
        TRandom3 r6;
        TRandom3 r7;
        TRandom3 r8;
        TRandom3 r9;
        TFile * output_file;
        TTree * tree;
        double ler_e;
        double her_e;
        double x_angle;

        int total_ks_candidates;
        int total_after_pt;
        int total_after_p;
        int total_after_cosTheta;
        int total_after_lepton;
        int total_after_pid;

        struct var_collection{
            int evtNo;
            int runNo;
            int expNo;
            // event var. used in hadronic selection
            double Evis;
            double Esum;
            double Psum;
            double Pz;
            double R2;
            double HeavyJetMass;
            double thrust[3];
            int Ntrk;
            int Ncls;

            double sqrts;
        } m_info;
        /*
        struct mc_var_collection{
            double thrust_truth[3];
        } m_mc_info; */

        struct KinematicsVars {
            double ler_e;
            double her_e;
            double x_angle;
            HepLorentzVector firstElectronCM;
            HepLorentzVector secondElectronCM;
            Hep3Vector CMBoost;
            HepLorentzVector cm;
            double thrust;
            double thrust_theta;
            double thrust_phi;
            double sqrts;
        } kinematics;

        Vdouble pip_E_cms;
        Vdouble pip_px_cms;
        Vdouble pip_py_cms;
        Vdouble pip_pz_cms;
        Vdouble pim_E_cms;
        Vdouble pim_px_cms;
        Vdouble pim_py_cms;
        Vdouble pim_pz_cms;
        Vdouble pip_p;
        Vdouble pim_p;
        Vdouble pip_costheta;
        Vdouble pim_costheta;

        // pion PID variables (aligned with pip/pim vectors)
        Vdouble pip_atcKPi;
        Vdouble pim_atcKPi;

        Vdouble ks_m_combine;
        Vdouble ks_m_read;
};


extern "C" Module_descr *mdcl_KsSpinAlignment_NullTest() {
    KsSpinAlignment_NullTest *module = new KsSpinAlignment_NullTest;
    Module_descr *dscr = new Module_descr( "KsSpinAlignment_NullTest", module );
    dscr->define_param( "output_filename", "output root file name",255 ,(char*)module->output_filename );
    return dscr;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif