#include "belle.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <map>

#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "particle/Particle.h"
#include "particle/utility.h"
#include "particle/combination.h"

#include "./KsSpinAlignment_NullTest.h"
#include "./AnaConsts.h"
#include <eid/eid.h>
#include "TTree.h"

#include <panther/panther.h>
#include "tables/evtcls.h"

#include "ip/IpProfile.h"
#include <mdst/findLambda.h>
#include <mdst/Muid_mdst.h>
#include "mdst/mdst.h"
#include <kid/atc_pid.h>
#include "toolbox/Thrust.h"
#include "toolbox/FoxWolfr.h"
#include "benergy/BeamEnergy.h"
#include "math.h"
#include "toolbox/FuncPtr.h"
#include "TRandom.h"
#include <mdst/findKs.h> 

#include MDST_H
#include BELLETDF_H

using namespace std;

/*
-----------------------------------------
qqbar(udsc) -> K_S + anything                             
                \-> pi+ pi-                                 
measure K_S's spin alignment as null test

K_S is spin-0, so spin alignment should be zero
This serves as a validation of the analysis method

using pid efficiency from control sample
-----------------------------------------

version : v1.0.0
Date    : 2026.05.11
Author  : Zheng Wang
*/

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

Hep3Vector& retSelf(Hep3Vector& vec)
{
	return vec;
};


KsSpinAlignment_NullTest::KsSpinAlignment_NullTest(){
    output_filename=new char[256];

    r1.SetSeed(100);
    r5.SetSeed(500);
    r6.SetSeed(600);
    r7.SetSeed(700);
    r8.SetSeed(800);
    r9.SetSeed(900);

    countEvt = 0;
    countEvt_fail = 0;
    
    // Initialize cut flow statistics
    total_ks_candidates = 0;
    total_after_pt = 0;
    total_after_p = 0;
    total_after_cosTheta = 0;
    total_after_lepton = 0;
    total_after_pid = 0;
    
    return;
}

void KsSpinAlignment_NullTest::init(int *status){
    cout << "in KsSpinAlignment_NullTest::init()" << endl;
    hist_def();
    *status=0;
    return;
}

void KsSpinAlignment_NullTest::begin_run(BelleEvent* evptr, int* status){
    // loading calibration constants and lookup tables needed for particle ID ?
    eid::init_data();
    (void)evptr; (void)status;

    IpProfile::begin_run();
    BeamEnergy::begin_run();

    
    kinematics.ler_e = BeamEnergy::E_LER();
    kinematics.her_e = BeamEnergy::E_HER();
    kinematics.x_angle = BeamEnergy::Cross_angle();
    /*     
    kinematics.ler_e = 3.5;
    kinematics.her_e = 7.9965;
    kinematics.x_angle = 0.022;
    */
    if(kinematics.ler_e < 3.0 || kinematics.her_e < 7.0 || kinematics.ler_e > 5.0 || kinematics.her_e > 9.0){
        cout<<"someting is wrong ler_e is "<< kinematics.ler_e << " her_e is "<< kinematics.her_e <<endl;
        return;
    }

    kinematics.cm = HepLorentzVector(-kinematics.her_e * sin(kinematics.x_angle), 0., -kinematics.her_e * cos(kinematics.x_angle) + kinematics.ler_e, kinematics.her_e + kinematics.ler_e);
    kinematics.CMBoost = kinematics.cm.boostVector();
    kinematics.firstElectronCM = HepLorentzVector(kinematics.her_e * sin(kinematics.x_angle), 0., kinematics.her_e * cos(kinematics.x_angle), kinematics.her_e);
    kinematics.secondElectronCM = HepLorentzVector(0., 0., -kinematics.ler_e, kinematics.ler_e);
    kinematics.sqrts = kinematics.cm.m();

    return;
}


void KsSpinAlignment_NullTest::end_run(BelleEvent* evptr, int* status){
    (void)evptr; (void)status;
    return; 
}


void KsSpinAlignment_NullTest::hist_def(){
    output_file = new TFile(output_filename, "RECREATE");
    tree = new TTree("event", "event"); 
                                        
    // macro to add branches            
    #define ADDBRANCH__(name, var, type );   tree->Branch(#name, &m_info.var, #type)
    #define ADDBRANCH(x, type)             ADDBRANCH__(x, x, x/type)
    #define ADDBARRAY__(name, var, n, type) ADDBRANCH__(name, var, name[n]/type)
    #define ADDBARRAY(x, n, type)          ADDBARRAY__(x, x, n, type)

    ADDBRANCH(evtNo, I);
    ADDBRANCH(runNo, I);
    ADDBRANCH(expNo, I);
    // event var. used in hadronic selection
    ADDBRANCH(sqrts, D);

    // vector<double> type branches
    tree->Branch("pip_E_cms", &pip_E_cms);
    tree->Branch("pip_px_cms", &pip_px_cms);
    tree->Branch("pip_py_cms", &pip_py_cms);
    tree->Branch("pip_pz_cms", &pip_pz_cms);
    tree->Branch("pim_E_cms", &pim_E_cms);
    tree->Branch("pim_px_cms", &pim_px_cms);
    tree->Branch("pim_py_cms", &pim_py_cms);
    tree->Branch("pim_pz_cms", &pim_pz_cms);

    tree->Branch("pip_p", &pip_p);
    tree->Branch("pim_p", &pim_p);
    tree->Branch("pip_costheta", &pip_costheta);
    tree->Branch("pim_costheta", &pim_costheta);

    tree->Branch("pip_atcKPi", &pip_atcKPi);
    tree->Branch("pim_atcKPi", &pim_atcKPi);

    return;
}


void KsSpinAlignment_NullTest::event(BelleEvent* evptr, int* status){
    (void)evptr; (void)status;

    int  expNo=Belle_event_Manager::get_manager().begin()->ExpNo();
    int  evtNo=Belle_event_Manager::get_manager().begin()->EvtNo();
    int  runNo=Belle_event_Manager::get_manager().begin()->RunNo();

    //if(!(countEvt%10000)) 
    //cout << "evt " <<countEvt<< " expNo "<< expNo << ", runNo "<< runNo << ", evtNo "<< evtNo <<endl;
    if(!IpProfile::usable()) {
        cout <<" ip not usable ..." << endl;
        return;
    }
    countEvt ++;

    pip_E_cms.clear(); pip_px_cms.clear(); pip_py_cms.clear(); pip_pz_cms.clear();
    pim_E_cms.clear(); pim_px_cms.clear(); pim_py_cms.clear(); pim_pz_cms.clear();
    pip_p.clear(); pim_p.clear(); pip_costheta.clear(); pim_costheta.clear();
    pip_atcKPi.clear();
    pim_atcKPi.clear();

    // ------------------------------------------------------------------
    // Save ALL pi+/pi- tracks in the event with PID cut ONLY.
    // (No pt/p/cosTheta/dr/dz/lepton-veto/event-selection cuts.)
    // ------------------------------------------------------------------
    Mdst_charged_Manager& charged_mgr = Mdst_charged_Manager::get_manager();

    for (Mdst_charged_Manager::iterator ch_it = charged_mgr.begin();
         ch_it != charged_mgr.end(); ++ch_it) {

        const Mdst_charged& chg = *ch_it;
        const int charge = chg.charge();
        if (charge != 1 && charge != -1) continue;

        // PID: pion-like
        atc_pid selKPi(3, 1, 5, 3, 2);  // K/pi separation
        const double atcKPi = selKPi.prob(*ch_it);
        const bool isPion = (atcKPi < 0.4);
        if (!isPion) continue;

        Hep3Vector p3(chg.p(0), chg.p(1), chg.p(2));
        const double cosTheta = cos(p3.theta());

        // cosTheta acceptance cut
        if (cosTheta < cuts::trk_cosTheta_min || cosTheta > cuts::trk_cosTheta_max)
            continue;

        const double m_pi = xmass[2];
        const double E = sqrt(p3.mag2() + m_pi * m_pi);
        HepLorentzVector p4_boosted(p3, E);
        p4_boosted.boost(kinematics.CMBoost);

        if (charge == 1) {
            pip_E_cms.push_back(p4_boosted.e());
            pip_px_cms.push_back(p4_boosted.px());
            pip_py_cms.push_back(p4_boosted.py());
            pip_pz_cms.push_back(p4_boosted.pz());
            pip_p.push_back(p3.mag());
            pip_costheta.push_back(cosTheta);
            pip_atcKPi.push_back(atcKPi);
        } else {
            pim_E_cms.push_back(p4_boosted.e());
            pim_px_cms.push_back(p4_boosted.px());
            pim_py_cms.push_back(p4_boosted.py());
            pim_pz_cms.push_back(p4_boosted.pz());
            pim_p.push_back(p3.mag());
            pim_costheta.push_back(cosTheta);
            pim_atcKPi.push_back(atcKPi);
        }

    }

    // Save event identifiers
    m_info.evtNo = evtNo;
    m_info.runNo = runNo;
    m_info.expNo = expNo;
    m_info.sqrts = kinematics.sqrts;

    tree->Fill();

    return;
}



void KsSpinAlignment_NullTest::disp_stat(const char*){
    return;
} 


void KsSpinAlignment_NullTest::term(){
    output_file->Write();
    output_file->Close();
    return;
}


std::pair<double, double> KsSpinAlignment_NullTest::calculateHeavyJetMassEnergy(const Vp4& particles, const Hep3Vector& thrustAxis) {
    HepLorentzVector jet1, jet2;

    // Split particles into two jets based on their projection onto the thrust axis
    for (size_t i = 0; i < particles.size(); ++i) {
        HepLorentzVector p = particles[i];
        double projection = p.vect().dot(thrustAxis);
        if (projection > 0) {
            jet1 += p;
        } else {
            jet2 += p;
        }
    }

    double mJet1 = jet1.m();
    double mJet2 = jet2.m();
    double eJet1 = jet1.e();
    double eJet2 = jet2.e();

    // Return heavier jet mass and energy
    if (mJet1 > mJet2) {
        return std::make_pair(mJet1, eJet1);
    } else {
        return std::make_pair(mJet2, eJet2);
    }
}


std::pair<double, double> KsSpinAlignment_NullTest::calculateSphericityAplanarity(const std::vector<HepLorentzVector>& particles) {
    // Initialize the momentum tensor
    double M[3][3] = {0};

    // Fill the momentum tensor
    for (size_t i = 0; i < particles.size(); ++i) {
        HepLorentzVector p = particles[i];
        Hep3Vector vec = p.vect(); // Extract the 3-momentum
        M[0][0] += vec.x() * vec.x();
        M[0][1] += vec.x() * vec.y();
        M[0][2] += vec.x() * vec.z();
        M[1][1] += vec.y() * vec.y();
        M[1][2] += vec.y() * vec.z();
        M[2][2] += vec.z() * vec.z();
    }

    // Symmetrize the tensor
    M[1][0] = M[0][1];
    M[2][0] = M[0][2];
    M[2][1] = M[1][2];

    // Normalize the tensor
    double trace = M[0][0] + M[1][1] + M[2][2];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            M[i][j] /= trace;
        }
    }

    // Compute eigenvalues of the tensor
    double eigenvalues[3] = {0};
    // Use a numerical method or library to compute eigenvalues
    // Manually compute eigenvalues of the symmetric 3x3 matrix M
    // Characteristic equation: |M - λI| = 0
    double p1 = M[0][1]*M[0][1] + M[0][2]*M[0][2] + M[1][2]*M[1][2];
    if (p1 == 0) {
        // Diagonal matrix
        eigenvalues[0] = M[0][0];
        eigenvalues[1] = M[1][1];
        eigenvalues[2] = M[2][2];
    } else {
        double q = (M[0][0] + M[1][1] + M[2][2]) / 3.0;
        double p2 = (M[0][0] - q)*(M[0][0] - q) + (M[1][1] - q)*(M[1][1] - q) + (M[2][2] - q)*(M[2][2] - q) + 2.0*p1;
        double p = sqrt(p2 / 6.0);

        // B = (1 / p) * (M - q * I)
        double B[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                B[i][j] = M[i][j] - (i == j ? q : 0.0);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                B[i][j] /= p;

        // Compute determinant of B
        double r = 0.5 * (
            B[0][0]*B[1][1]*B[2][2] +
            B[0][1]*B[1][2]*B[2][0] +
            B[0][2]*B[1][0]*B[2][1] -
            B[0][2]*B[1][1]*B[2][0] -
            B[0][1]*B[1][0]*B[2][2] -
            B[0][0]*B[1][2]*B[2][1]
        );

        // Clamp r to [-1, 1] to avoid numerical errors
        r = std::max(-1.0, std::min(1.0, r));

        // Compute the eigenvalues
        double phi = acos(r) / 3.0;
        eigenvalues[0] = q + 2.0 * p * cos(phi);
        eigenvalues[1] = q + 2.0 * p * cos(phi + (2.0*M_PI/3.0));
        eigenvalues[2] = q + 2.0 * p * cos(phi + (4.0*M_PI/3.0));
    }

    // Sort eigenvalues in ascending order
    std::sort(eigenvalues, eigenvalues + 3);

    // Calculate Sphericity and Aplanarity
    double sphericity = 1.5 * (eigenvalues[0] + eigenvalues[1]);
    double aplanarity = 1.5 * eigenvalues[0];

    return std::make_pair(sphericity, aplanarity);
}

 void KsSpinAlignment_NullTest::getDrDz(Mdst_charged_Manager::iterator chr_it, int masshyp, double& dr, double& dz, double& refitPx, double& refitPy, double& refitPz){
    Mdst_trk& mdsttrk = chr_it->trk();
    Mdst_trk_fit &mdsttrkfit=mdsttrk.mhyp(masshyp);
    HepPoint3D pivot(mdsttrkfit.pivot_x(),mdsttrkfit.pivot_y(),mdsttrkfit.pivot_z());

    HepVector a( 5, 0 );
    a[0] = mdsttrkfit.helix( 0 );
    a[1] = mdsttrkfit.helix( 1 );
    a[2] = mdsttrkfit.helix( 2 );
    a[3] = mdsttrkfit.helix( 3 );
    a[4] = mdsttrkfit.helix( 4 );

    HepSymMatrix Ea( 5, 0 );
    Ea[0][0] = mdsttrkfit.error( 0 );
    Ea[1][0] = mdsttrkfit.error( 1 );
    Ea[1][1] = mdsttrkfit.error( 2 );
    Ea[2][0] = mdsttrkfit.error( 3 );
    Ea[2][1] = mdsttrkfit.error( 4 );
    Ea[2][2] = mdsttrkfit.error( 5 );

  Ea[3][0] = mdsttrkfit.error( 6 );
    Ea[3][1] = mdsttrkfit.error( 7 );
    Ea[3][2] = mdsttrkfit.error( 8 );
    Ea[3][3] = mdsttrkfit.error( 9 );
    Ea[4][0] = mdsttrkfit.error( 10 );
    Ea[4][1] = mdsttrkfit.error( 11 );
    Ea[4][2] = mdsttrkfit.error( 12 );
    Ea[4][3] = mdsttrkfit.error( 13 );
    Ea[4][4] = mdsttrkfit.error( 14 );


    Helix helix( pivot, a, Ea );
    helix.pivot( IpProfile::position(1));
    refitPx=helix.momentum().x();
    refitPy=helix.momentum().y();
    refitPz=helix.momentum().z();

    dr  = helix.dr();
    dz  = helix.dz();
}


void KsSpinAlignment_NullTest::other(int* , BelleEvent*, int* ){
    return;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif


// version log 
// v1.0.0 : 
// - first version
// Dec. 24, 2025
// v1.0.1 :
// add continue statement after warning message for K_S decay daughters number not equal to 2
// the crash is due to gen_hepevt(chg) return nullptr || gen.mother() , add check fix this issue. when apply mc truth matching

// v2.0.0 :
// recontruct Ks using official goodKs selection ; and fixed little bug in save thrust truth in readMC()
// Dec. 30, 2025

// v2.1.0 :
// put additional cuts on Ks daughter tracks: pt, cosTheta, lepton veto, pion PID
// add cumulative and final cut flow statistics printout for Ks selection

// v2.2.0 :
// change Ks reconstruct method from nisKsfind to  KSfinder; 
// and cut track's momentum p > 0.5 GeV; change pt cut to 0.1 GeV to keep it same as in HadronBJ reqiurement
// Jan. 19, 2026

// v2.3.0 :
// fix little error on isSignal saving 
// add topology info fill function
// v2.3.1 :
// change back to nisKsFind for test
// Jan. 20, 2026

// v2.3.2 :
// use FindKs; not cut p > 0.5 GeV ; save pion's costheta ; and change MAX_MC_PARTICLES to 80
// Jan. 22, 2026

// v2.3.3 : 
// fix isSignal definition: 0 not from Ks/K0; 1 pi from Ks/K0 directly; 
// 2 mu from pi with grandmother Ks/K0; and add check for gen.mother() when apply mc truth matching to avoid crash
// Mar. 11, 2026