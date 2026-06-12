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
#include HEPEVT_H
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
        void readMC();
        void other(int*, BelleEvent*, int*);
        std::pair<double, double> calculateHeavyJetMassEnergy(const Vp4 &particles, const Hep3Vector &thrustAxis);
        std::pair<double, double> calculateSphericityAplanarity(const std::vector<HepLorentzVector>& particles);
        void fillMCTruthForTopo();

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
        TTree * mctree;
        double ler_e;
        double her_e;
        double x_angle;

        double n_kstar_truth; // added for truth K*(892) count

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
            
            double thrust_truth[3];
            double qqbar_axis[2];

        } m_info;

        Vdouble thrust_truth;
        Vdouble qqbar_axis;
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

        Vdouble kp_E_cms;
        Vdouble kp_px_cms;
        Vdouble kp_py_cms;
        Vdouble kp_pz_cms;
        Vdouble km_E_cms;
        Vdouble km_px_cms;
        Vdouble km_py_cms;
        Vdouble km_pz_cms;
        Vdouble kp_p;
        Vdouble km_p;
        Vdouble kp_costheta;
        Vdouble km_costheta;

        Vbool kp_isSignal;
        Vbool km_isSignal;

        Vdouble kp_E_cms_gen;
        Vdouble kp_px_cms_gen;
        Vdouble kp_py_cms_gen;
        Vdouble kp_pz_cms_gen;
        Vdouble km_E_cms_gen;
        Vdouble km_px_cms_gen;
        Vdouble km_py_cms_gen;
        Vdouble km_pz_cms_gen;

        Vdouble kp_E_cms_truth;
        Vdouble kp_px_cms_truth;
        Vdouble kp_py_cms_truth;
        Vdouble kp_pz_cms_truth;
        Vdouble km_E_cms_truth;
        Vdouble km_px_cms_truth;
        Vdouble km_py_cms_truth;
        Vdouble km_pz_cms_truth;

        Vdouble kp_theta;
        Vdouble km_theta;
        Vdouble kp_p_truth;
        Vdouble km_p_truth;
        Vdouble kp_pt_truth;
        Vdouble km_pt_truth;

        // pion daughters of K*(892): pi- from K*0, pi+ from anti-K*0
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

        Vbool pip_isSignal;
        Vbool pim_isSignal;

        Vdouble pip_E_cms_gen;
        Vdouble pip_px_cms_gen;
        Vdouble pip_py_cms_gen;
        Vdouble pip_pz_cms_gen;
        Vdouble pim_E_cms_gen;
        Vdouble pim_px_cms_gen;
        Vdouble pim_py_cms_gen;
        Vdouble pim_pz_cms_gen;

        Vdouble pip_E_cms_truth;
        Vdouble pip_px_cms_truth;
        Vdouble pip_py_cms_truth;
        Vdouble pip_pz_cms_truth;
        Vdouble pim_E_cms_truth;
        Vdouble pim_px_cms_truth;
        Vdouble pim_py_cms_truth;
        Vdouble pim_pz_cms_truth;

        Vdouble pip_theta;
        Vdouble pim_theta;
        Vdouble pip_p_truth;
        Vdouble pim_p_truth;
        Vdouble pip_pt_truth;
        Vdouble pim_pt_truth;

        bool isMC;

        static const int MAX_MC_PARTICLES = 80;  // Adjust as needed
        int m_nMCGen;
        double m_MCGenPDG[MAX_MC_PARTICLES];
        double m_MCGenMothIndex[MAX_MC_PARTICLES];
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