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
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<Hep3Vector> Vp3;

class HadronB_dataMC_comp : public Module{
    public:
        HadronB_dataMC_comp();
        ~HadronB_dataMC_comp() {}
        void init(int *dummy);
        void begin_run(BelleEvent* evptr, int* status);
        void end_run(BelleEvent* evptr, int* status);
        void hist_def();
        void event(BelleEvent* evpter, int* status);
        void disp_stat(const char*);
        void term();
        void other(int*, BelleEvent*, int*);
        std::pair<double, double> calculateHeavyJetMassEnergy(const Vp4 &particles, const Hep3Vector &thrustAxis);
        std::pair<double, double> calculateSphericityAplanarity(const std::vector<HepLorentzVector>& particles);

    public: // BASF parameter
        char* output_filename;

    private:
        int countEvt;
        TFile * output_file;
        
        // Histograms for event variables
        TH1F *h_Evis;
        TH1F *h_HeavyJetMass_over_Evis;
        TH1F *h_Esum;
        TH1F *h_Psum;
        TH1F *h_Pz;
        TH1F *h_R2;
        TH1F *h_HeavyJetMass;
        TH1F *h_Ntrk;
        TH1F *h_Ncls;
        TH1F *h_thrust;
        TH1F *h_thrust_costheta;
        TH1F *h_thrust_phi;
        TH1F *h_HeavyJetEnergy;
        TH1F *h_Sphericity;
        TH1F *h_Aplanarity;
        
        // Histograms for all charged tracks (lab frame)
        TH1F *h_trk_dr;
        TH1F *h_trk_dz;
        TH1F *h_trk_p;
        TH1F *h_trk_costheta;
        TH1F *h_trk_phi;
        
        // Histograms for event-level vertex info
        TH1F *h_PrimeR;
        TH1F *h_PrimeZ;
        
        // Histograms for photons
        TH1F *h_photon_p;
        TH1F *h_photon_costheta;
        TH1F *h_photon_phi;

        struct KinematicsVars {
            double ler_e;
            double her_e;
            double x_angle;
            Hep3Vector CMBoost;
            HepLorentzVector cm;
            double sqrts;
        } kinematics;

        // Track variables (lab frame)
        Vdouble trk_dr;
        Vdouble trk_dz;
        Vdouble trk_p;
        Vdouble trk_costheta;
        Vdouble trk_phi;
        Vint trk_charge;
        
        // Primary vertex variables
        double PrimeR;
        double PrimeZ;
        
        // Photon variables (lab frame)
        Vdouble photon_p;
        Vdouble photon_costheta;
        Vdouble photon_phi;
};

extern "C" Module_descr *mdcl_HadronB_dataMC_comp() {
    HadronB_dataMC_comp *module = new HadronB_dataMC_comp;
    Module_descr *dscr = new Module_descr( "HadronB_dataMC_comp", module );
    dscr->define_param( "output_filename", "output root file name",255 ,(char*)module->output_filename );
    return dscr;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif