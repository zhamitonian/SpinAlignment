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
#include HEPEVT_H

using namespace std;

/*
-----------------------------------------
qqbar(udsc) -> K_S + anything                             
                \-> pi+ pi-                                 
measure K_S's spin alignment as null test

K_S is spin-0, so spin alignment should be zero
This serves as a validation of the analysis method
-----------------------------------------

version : v2.3.2
Date    : 2026.01.20
Author  : Zhen Wang
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

    isMC = true;
    
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
    mctree = new TTree("truth", "truth");
                                        
    // macro to add branches            
    #define ADDBRANCH__(name, var, type );   tree->Branch(#name, &m_info.var, #type)
    #define ADDBRANCH(x, type)             ADDBRANCH__(x, x, x/type)
    #define ADDBARRAY__(name, var, n, type) ADDBRANCH__(name, var, name[n]/type)
    #define ADDBARRAY(x, n, type)          ADDBARRAY__(x, x, n, type)

    ADDBRANCH(evtNo, I);
    ADDBRANCH(runNo, I);
    ADDBRANCH(expNo, I);
    // event var. used in hadronic selection
    ADDBRANCH(Evis, D);
    ADDBRANCH(Esum, D);
    ADDBRANCH(Psum, D);
    ADDBRANCH(Pz, D);
    ADDBRANCH(R2, D);
    ADDBRANCH(HeavyJetMass, D);
    ADDBRANCH(Ntrk, I);
    ADDBRANCH(Ncls, I);

    ADDBRANCH(sqrts, D);

    ADDBARRAY(thrust, 3, D);
   
    if (isMC) {
        tree->Branch("nMCGen", &m_nMCGen, "nMCGen/I");
        for (int i = 0; i < MAX_MC_PARTICLES; i++) {
            tree->Branch(Form("MCGenPDG_%d", i), &m_MCGenPDG[i], Form("MCGenPDG_%d/D", i));
            tree->Branch(Form("MCGenMothIndex_%d", i), &m_MCGenMothIndex[i], Form("MCGenMothIndex_%d/D", i));
        }
    }

    if (isMC){
        mctree->Branch("pip_E_cms_truth", &pip_E_cms_truth);
        mctree->Branch("pip_px_cms_truth", &pip_px_cms_truth);
        mctree->Branch("pip_py_cms_truth", &pip_py_cms_truth);
        mctree->Branch("pip_pz_cms_truth", &pip_pz_cms_truth);
        mctree->Branch("pim_E_cms_truth", &pim_E_cms_truth);
        mctree->Branch("pim_px_cms_truth", &pim_px_cms_truth);
        mctree->Branch("pim_py_cms_truth", &pim_py_cms_truth);
        mctree->Branch("pim_pz_cms_truth", &pim_pz_cms_truth);

        mctree->Branch("thrust_truth", &thrust_truth);
        mctree->Branch("qqbar_axis", &qqbar_axis);

        mctree->Branch("pip_theta", &pip_theta);
        mctree->Branch("pim_theta", &pim_theta);
        mctree->Branch("pip_p_truth", &pip_p_truth);
        mctree->Branch("pim_p_truth", &pim_p_truth);
        mctree->Branch("pip_pt_truth", &pip_pt_truth);
        mctree->Branch("pim_pt_truth", &pim_pt_truth);
    }

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

    //tree->Branch("ks_m_combine", &ks_m_combine);
    //tree->Branch("ks_m_read", &ks_m_read);

    if(isMC){
        tree->Branch("pip_E_cms_gen", &pip_E_cms_gen);
        tree->Branch("pip_px_cms_gen", &pip_px_cms_gen);
        tree->Branch("pip_py_cms_gen", &pip_py_cms_gen);
        tree->Branch("pip_pz_cms_gen", &pip_pz_cms_gen);
        tree->Branch("pim_E_cms_gen", &pim_E_cms_gen);
        tree->Branch("pim_px_cms_gen", &pim_px_cms_gen);
        tree->Branch("pim_py_cms_gen", &pim_py_cms_gen);
        tree->Branch("pim_pz_cms_gen", &pim_pz_cms_gen);

        tree->Branch("pip_isSignal", &pip_isSignal);
        tree->Branch("pim_isSignal", &pim_isSignal);

        tree->Branch("n_ks_truth", &n_ks_truth);

        tree->Branch("thrust_truth", &thrust_truth);
        tree->Branch("qqbar_axis", &qqbar_axis);
    }
    return;
}


void KsSpinAlignment_NullTest::event(BelleEvent* evptr, int* status){
    (void)evptr; (void)status;

    int  expNo=Belle_event_Manager::get_manager().begin()->ExpNo();
    int  evtNo=Belle_event_Manager::get_manager().begin()->EvtNo();
    int  runNo=Belle_event_Manager::get_manager().begin()->RunNo();

    const HepPoint3D ip_position = IpProfile::e_position();

    Evtcls_hadronic_flag_Manager& ClassMgr = Evtcls_hadronic_flag_Manager::get_manager();

    int EvtClsFlgB = ClassMgr.begin()->hadronic_flag(2);

    if(!EvtClsFlgB){
        //cout << "event " << evtNo << " in run " << runNo << " exp " << expNo << " is not hadronB skimed data!" << endl;
        //cout << EvtClsFlgB << endl;
        return; // only process hadronB skimed data)
    }

    if (isMC == 1){
        thrust_truth.clear();
        qqbar_axis.clear();
        readMC();
    }

    //if(!(countEvt%10000)) 
    //cout << "evt " <<countEvt<< " expNo "<< expNo << ", runNo "<< runNo << ", evtNo "<< evtNo <<endl;
    if(!IpProfile::usable()) {
        cout <<" ip not usable ..." << endl;
        return;
    }
    countEvt ++;

    Mdst_charged_Manager& charged_mgr = Mdst_charged_Manager::get_manager();
    Mdst_gamma_Manager& gamma_mgr = Mdst_gamma_Manager::get_manager();

    //---------------------------------------------------------
    Vp3 particles_p3_cms;
    Vp4 particles_p4_cms;
    particles_p3_cms.clear();
    particles_p4_cms.clear();

    // ----------------------------------------------------------------------------------


    Evtcls_hadron_charged_Manager& hadron_chg_mgr = Evtcls_hadron_charged_Manager::get_manager();
    for(Evtcls_hadron_charged_Manager::iterator it = hadron_chg_mgr.begin();
        it != hadron_chg_mgr.end(); it++){

        if (!it->charged()) continue;  // Check if charged is valid

        const Mdst_charged& Charged = it->charged();
        Hep3Vector P(Charged.px(), Charged.py(), Charged.pz());
        
        double E = sqrt(P.mag2() + xmass[2]*xmass[2]);
        HepLorentzVector P_cms(P,E);
        P_cms.boost(kinematics.CMBoost);
        
        particles_p3_cms.push_back(P_cms.vect());
        particles_p4_cms.push_back(P_cms);
    }

    // 对 gamma 也类似
    Evtcls_hadron_neutral_Manager& hadron_neu_mgr = Evtcls_hadron_neutral_Manager::get_manager();

    for(Evtcls_hadron_neutral_Manager::iterator it = hadron_neu_mgr.begin();
        it != hadron_neu_mgr.end(); it++){
        
        if (!it->gamma()) continue;  // Check if gamma is valid
        
        const Mdst_gamma& Gamma = it->gamma();
        Hep3Vector P(Gamma.px(), Gamma.py(), Gamma.pz());
        
        double E = P.mag();
        HepLorentzVector vec_p4(P, E);
        vec_p4.boost(kinematics.CMBoost);
        
        particles_p3_cms.push_back(vec_p4.vect());
        particles_p4_cms.push_back(vec_p4);
    }

    //Hep3Vector thr = thrust(particles_p3_cms.begin(), particles_p3_cms.end(), retSelf);
    Hep3Vector thr = thrust(particles_p3_cms.begin(), particles_p3_cms.end(), SelfFunc(Hep3Vector()));
    Hep3Vector tn  = thr.unit();
    double ThrParam = thr.mag();

    // ----------------------------------------------------------------------------------
    // Additional hadron event selection (for HadronB skimed data)

    Evtcls_hadron_info_Manager& hadronInfo_mgr = Evtcls_hadron_info_Manager::get_manager();
    double Evis = hadronInfo_mgr.begin()->Evis();
    double Esum = hadronInfo_mgr.begin()->Esum();
    double Psum = hadronInfo_mgr.begin()->Psum();
    double Pz   = hadronInfo_mgr.begin()->Pz();
    double R2   = hadronInfo_mgr.begin()->R2();
    double HeavyJetMass = hadronInfo_mgr.begin()->HeavyJetMass();
    int Ntrk = hadronInfo_mgr.begin()->Ntrk();
    int Ncls = hadronInfo_mgr.begin()->Ncls();
    double Thrust = hadronInfo_mgr.begin()->Thrust();
    if(Evis < cuts::Evis)
        return;
    //if(abs(Thrust - ThrParam) > 0.0001)
    //cout << "read: " << Thrust << " calc: " << ThrParam << endl; 
    // ----------------------------------------------------------------------------------
    
    pip_E_cms.clear(); pip_px_cms.clear(); pip_py_cms.clear(); pip_pz_cms.clear();
    pim_E_cms.clear(); pim_px_cms.clear(); pim_py_cms.clear(); pim_pz_cms.clear();
    pip_p.clear(); pim_p.clear(); pip_costheta.clear(); pim_costheta.clear();
    //ks_m_combine.clear(); ks_m_read.clear();
    if (isMC) {
        pip_isSignal.clear(); pim_isSignal.clear();
        pip_E_cms_gen.clear(); pip_px_cms_gen.clear(); pip_py_cms_gen.clear(); pip_pz_cms_gen.clear();
        pim_E_cms_gen.clear(); pim_px_cms_gen.clear(); pim_py_cms_gen.clear(); pim_pz_cms_gen.clear();
    }

    // using findKs to select good Ks from candidate in mdst_vee2 
    // see /home/belle/liyb/work/published/Xic2930/mode3/myutils.cc
    Mdst_vee2_Manager& VeeMgr = Mdst_vee2_Manager::get_manager();

    std::vector<Particle> kslist;
    
    for(std::vector<Mdst_vee2>::const_iterator iv=VeeMgr.begin(); iv!=VeeMgr.end(); iv++) {
        // Is it a K_short, or another V-particle ?
        if (iv->kind() != 1) continue;
        // Check for daughters
        if(!iv->chgd(0) || !iv->chgd(1)) continue;
        Particle Kshort(*iv);

        // use Fang Fang's Kshort finder to check K-short quality
        FindKs KSfinder;
        KSfinder.candidates(*iv, ip_position);
        int goodKsFlag = KSfinder.goodKs();

        if(goodKsFlag == 0) continue;     

        //saveKsInfo(Kshort,ksnb.fl());
        kslist.push_back(Kshort);
    }

      // 新增：保存Ks顶点坐标
    std::vector<double> ks_vx, ks_vy, ks_vz;

    // Accumulate cut flow statistics
    int n_total_ks = kslist.size();
    int n_after_pt = 0;
    int n_after_p = 0;
    int n_after_cosTheta = 0;
    int n_after_lepton = 0;
    int n_after_pid = 0;
    
    total_ks_candidates += n_total_ks;

    for (size_t i = 0; i < kslist.size(); ++i) {
        const Particle& Kshort = kslist[i];
        // 保存顶点坐标
        ks_vx.push_back(Kshort.momentum().decayVertex().x());
        ks_vy.push_back(Kshort.momentum().decayVertex().y());
        ks_vz.push_back(Kshort.momentum().decayVertex().z());

        // save info of Ks decay daughters pi+ pi-
        const Mdst_charged& dau1 = Kshort.child(0).mdstCharged();
        const Mdst_charged& dau2 = Kshort.child(1).mdstCharged();

        Hep3Vector p1(dau1.px(), dau1.py(), dau1.pz());
        Hep3Vector p2(dau2.px(), dau2.py(), dau2.pz());
        int q1 = dau1.charge();
        int q2 = dau2.charge();
        
        // Apply track quality cuts on Ks daughters (same as previous version, but no dr/dz cut)
        // Cut on pt and cosTheta
        double pt1 = p1.perp();
        double pt2 = p2.perp();
        double cosTheta1 = cos(p1.theta());
        double cosTheta2 = cos(p2.theta());
        
        if (pt1 < cuts::trkPt || pt2 < cuts::trkPt) 
            continue;
        n_after_pt++;
    
        //if (p1.mag() < cuts::trkP || p2.mag() < cuts::trkP)
        //    continue;
        n_after_p++;    
       
        if (cosTheta1 < cuts::trk_cosTheta_min || cosTheta1 > cuts::trk_cosTheta_max)
            continue;
        if (cosTheta2 < cuts::trk_cosTheta_min || cosTheta2 > cuts::trk_cosTheta_max)
            continue;
        n_after_cosTheta++;
        
        // Lepton identification: reject electrons and muons
        eid sel_e1(dau1);
        eid sel_e2(dau2);
        Muid_mdst muID1(dau1);
        Muid_mdst muID2(dau2);
        
        double e_id1 = sel_e1.prob(3, -1, 5);
        double e_id2 = sel_e2.prob(3, -1, 5);
        double mu_id1 = (muID1.Chi_2() > 0) ? muID1.Muon_likelihood() : 0;
        double mu_id2 = (muID2.Chi_2() > 0) ? muID2.Muon_likelihood() : 0;
        
        float e_cut = 0.85;
        float mu_cut = 0.9;
        
        // Reject if either daughter is identified as electron or muon
        bool isLepton1 = (mu_id1 > mu_cut) || (e_id1 > e_cut && mu_id1 < mu_cut);
        bool isLepton2 = (mu_id2 > mu_cut) || (e_id2 > e_cut && mu_id2 < mu_cut);
        
        if (isLepton1 || isLepton2)
            continue;
        n_after_lepton++;
        
        // PID check: ensure daughters are identified as pions
        atc_pid selKPi1(3, 1, 5, 3, 2);  // K/pi separation for dau1
        atc_pid selKPi2(3, 1, 5, 3, 2);  // K/pi separation for dau2
        atc_pid selPiP1(3, 1, 5, 2, 4);  // pi/proton separation for dau1
        atc_pid selPiP2(3, 1, 5, 2, 4);  // pi/proton separation for dau2
        
        double atcKPi1 = selKPi1.prob(dau1);
        double atcKPi2 = selKPi2.prob(dau2);
        double atcPiP1 = selPiP1.prob(dau1);
        double atcPiP2 = selPiP2.prob(dau2);
        
        // Require pion-like PID: atcKPi < 0.4 (not kaon) && atcPiP > 0.2 (not proton)
        bool isPion1 = (atcKPi1 < 0.4 && atcPiP1 > 0.2);
        bool isPion2 = (atcKPi2 < 0.4 && atcPiP2 > 0.2);
        
        if (!isPion1 || !isPion2)
            continue;
        n_after_pid++;
            
        /* this retrieve the 4-p without vertex fit correction
        double m_pi = xmass[2];
        double E1 = sqrt(p1.mag2() + m_pi*m_pi);
        double E2 = sqrt(p2.mag2() + m_pi*m_pi);
        HepLorentzVector p4_1(p1, E1);
        HepLorentzVector p4_2(p2, E2);
        */

        HepLorentzVector p4_1 = Kshort.child(0).momentum().p();
        HepLorentzVector p4_2 = Kshort.child(1).momentum().p(); 
        p4_1.boost(kinematics.CMBoost);
        p4_2.boost(kinematics.CMBoost);

        // 按电荷分配pi+ pi-
        if (q1 == 1) {
            pip_E_cms.push_back(p4_1.e());
            pip_px_cms.push_back(p4_1.px());
            pip_py_cms.push_back(p4_1.py());
            pip_pz_cms.push_back(p4_1.pz());
            pip_p.push_back(p1.mag());
            pip_costheta.push_back(cosTheta1);
        } else if (q1 == -1) {
            pim_E_cms.push_back(p4_1.e());
            pim_px_cms.push_back(p4_1.px());
            pim_py_cms.push_back(p4_1.py());
            pim_pz_cms.push_back(p4_1.pz());
            pim_p.push_back(p1.mag());
            pim_costheta.push_back(cosTheta1);
        }
        if (q2 == 1) {
            pip_E_cms.push_back(p4_2.e());
            pip_px_cms.push_back(p4_2.px());
            pip_py_cms.push_back(p4_2.py());
            pip_pz_cms.push_back(p4_2.pz());
            pip_p.push_back(p2.mag());
            pip_costheta.push_back(cosTheta2);
        } else if (q2 == -1) {
            pim_E_cms.push_back(p4_2.e());
            pim_px_cms.push_back(p4_2.px());
            pim_py_cms.push_back(p4_2.py());
            pim_pz_cms.push_back(p4_2.pz());
            pim_p.push_back(p2.mag());
            pim_costheta.push_back(cosTheta2);
        }
        
        // ----- test coordinate Ks reconstruction -----
        /* the difference is less than 10-3 GeV
        HepLorentzVector p4_combine = p4_1 + p4_2;
        HepLorentzVector p4_Ks = Kshort.momentum().p();
        p4_Ks.boost(kinematics.CMBoost);
        cout << "Ks mass from mdst_vee2: " << p4_Ks.m() << ", mass from combine daughters: " << p4_combine.m() << endl;
        ks_m_combine.push_back(p4_combine.m());
        ks_m_read.push_back(p4_Ks.m()); 
        */

        // ----------- save MC generate level info -----------
        if (isMC) {
            // dau1
            if (get_hepevt(dau1)) {
                const Gen_hepevt& gen = gen_level(get_hepevt(dau1));
                if (gen.mother() && abs(gen.mother().idhep()) == 310 && abs(gen.idhep()) == 211) {
                    HepLorentzVector p4_truth(gen.PX(), gen.PY(), gen.PZ(), gen.E());
                    p4_truth.boost(kinematics.CMBoost);
                    if (q1 == 1) {
                        pip_E_cms_gen.push_back(p4_truth.e());
                        pip_px_cms_gen.push_back(p4_truth.px());
                        pip_py_cms_gen.push_back(p4_truth.py());
                        pip_pz_cms_gen.push_back(p4_truth.pz());
                        pip_isSignal.push_back(true);
                    } else if (q1 == -1) {
                        pim_E_cms_gen.push_back(p4_truth.e());
                        pim_px_cms_gen.push_back(p4_truth.px());
                        pim_py_cms_gen.push_back(p4_truth.py());
                        pim_pz_cms_gen.push_back(p4_truth.pz());
                        pim_isSignal.push_back(true);
                    }
                } else {
                    if (q1 == 1) pip_isSignal.push_back(false);
                    else if (q1 == -1) pim_isSignal.push_back(false);
                }
            } else {
                if (q1 == 1) pip_isSignal.push_back(false);
                else if (q1 == -1) pim_isSignal.push_back(false);
            }
            // dau2
            if (get_hepevt(dau2)) {
                const Gen_hepevt& gen = gen_level(get_hepevt(dau2));
                if (gen.mother() && abs(gen.mother().idhep()) == 310 && abs(gen.idhep()) == 211) {
                    HepLorentzVector p4_truth(gen.PX(), gen.PY(), gen.PZ(), gen.E());
                    p4_truth.boost(kinematics.CMBoost);
                    if (q2 == 1) {
                        pip_E_cms_gen.push_back(p4_truth.e());
                        pip_px_cms_gen.push_back(p4_truth.px());
                        pip_py_cms_gen.push_back(p4_truth.py());
                        pip_pz_cms_gen.push_back(p4_truth.pz());
                        pip_isSignal.push_back(true);
                    } else if (q2 == -1) {
                        pim_E_cms_gen.push_back(p4_truth.e());
                        pim_px_cms_gen.push_back(p4_truth.px());
                        pim_py_cms_gen.push_back(p4_truth.py());
                        pim_pz_cms_gen.push_back(p4_truth.pz());
                        pim_isSignal.push_back(true);
                    }
                } else {
                    if (q2 == 1) pip_isSignal.push_back(false);
                    else if (q2 == -1) pim_isSignal.push_back(false);
                }
            } else {
                if (q2 == 1) pip_isSignal.push_back(false);
                else if (q2 == -1) pim_isSignal.push_back(false);
            }
        }
        // ------------------------------------------------------
    }

    // Accumulate to total statistics
    total_after_pt += n_after_pt;
    total_after_p += n_after_p; 
    total_after_cosTheta += n_after_cosTheta;
    total_after_lepton += n_after_lepton;
    total_after_pid += n_after_pid;
    
    // Print cumulative cut flow statistics every 10000 events
    if(!(countEvt%10000)) {
        cout << "\n=== Cumulative Ks Cut Flow Statistics (up to Event " << countEvt << ") ===" << endl;
        cout << "Total Ks candidates (after FindKs): " << total_ks_candidates << endl;
        cout << "After pt cut:                            " << total_after_pt 
             << " (rejected: " << total_ks_candidates - total_after_pt 
             << ", eff: " << (total_ks_candidates>0 ? 100.0*total_after_pt/total_ks_candidates : 0) << "%)" << endl;
        cout << "After p cut:                             " << total_after_p 
             << " (rejected: " << total_after_pt - total_after_p 
             << ", eff: " << (total_after_pt>0 ? 100.0*total_after_p/total_after_pt : 0) << "%)" << endl;
        cout << "After cosTheta cut:                      " << total_after_cosTheta 
             << " (rejected: " << total_after_pt - total_after_cosTheta 
             << ", eff: " << (total_after_pt>0 ? 100.0*total_after_cosTheta/total_after_pt : 0) << "%)" << endl;
        cout << "After lepton veto:                       " << total_after_lepton 
             << " (rejected: " << total_after_cosTheta - total_after_lepton 
             << ", eff: " << (total_after_cosTheta>0 ? 100.0*total_after_lepton/total_after_cosTheta : 0) << "%)" << endl;
        cout << "After pion PID cut:                      " << total_after_pid 
             << " (rejected: " << total_after_lepton - total_after_pid 
             << ", eff: " << (total_after_lepton>0 ? 100.0*total_after_pid/total_after_lepton : 0) << "%)" << endl;
        cout << "Overall efficiency:                      " 
             << (total_ks_candidates>0 ? 100.0*total_after_pid/total_ks_candidates : 0) << "%" << endl;
        cout << "================================================================" << endl;
    }
 
    // require at least one K_S candidate
    if (pip_px_cms.size() * pim_px_cms.size() == 0) return;

    m_info.evtNo = evtNo;
    m_info.runNo = runNo;
    m_info.expNo = expNo;
    m_info.Evis = Evis;
    m_info.Esum = Esum;
    m_info.Psum = Psum;
    m_info.Pz   = Pz;
    m_info.R2   = R2;
    m_info.HeavyJetMass = HeavyJetMass;
    m_info.Ntrk = Ntrk;
    m_info.Ncls = Ncls;
    m_info.sqrts = kinematics.cm.m();
    double thrust_vec[3] = { thr.mag(), thr.z()/thr.mag(), thr.phi() };
    memcpy(m_info.thrust, thrust_vec, sizeof(thrust_vec));

    if (isMC) {
        fillMCTruthForTopo();
    }

    tree->Fill();

    return;
}


void KsSpinAlignment_NullTest::readMC()
{
    Vp3 allParticlesCMS_truth;

    allParticlesCMS_truth.clear();
    pip_E_cms_truth.clear();
    pip_px_cms_truth.clear();
    pip_py_cms_truth.clear();
    pip_pz_cms_truth.clear();
    pim_E_cms_truth.clear();
    pim_px_cms_truth.clear();
    pim_py_cms_truth.clear();
    pim_pz_cms_truth.clear();
    pip_theta.clear();
    pim_theta.clear();
    pip_p_truth.clear();
    pim_p_truth.clear();
    pip_pt_truth.clear();
    pim_pt_truth.clear();

    int Nquark = 0;
    HepLorentzVector q_momentum(0.,0.,0.,0.);

    Gen_hepevt_Manager& gen_hep_Mgr = Gen_hepevt_Manager::get_manager();

    for(Gen_hepevt_Manager::iterator gen_it = gen_hep_Mgr.begin(); gen_it != gen_hep_Mgr.end(); ++gen_it) {
        int pid = gen_it->idhep();
        int status = gen_it->isthep();
        
        // for daughter particle retrieval
        Gen_hepevt_Manager& gen_hepevt_mgr = Gen_hepevt_Manager::get_manager();

        if (pid == 310)  // K_S (PDG code 310)
        {
            int ndaug = (gen_it->daLast() - gen_it->daFirst()) + 1;
            //if (ndaug != 2) continue; // K_S -> pi+ pi- only
            if (ndaug != 2){
                cout << "Warning: K_S decay daughters number is " << ndaug << endl;
                continue;
            }

            for (int i = 0; i < ndaug; ++i){
                Panther_ID ID(gen_it->daFirst() + i);
                Gen_hepevt& daughter = gen_hepevt_mgr(ID);
                
                int daughter_pid = daughter.idhep();
                //cout << "daughter: "<< i <<" pid: " << daughter_pid << endl;
                if (abs(daughter_pid) != 211) continue; // only pi+ and pi-
                
                HepLorentzVector daughter_v4 (daughter.PX(), daughter.PY(), daughter.PZ(), daughter.E());
                double theta_lab = daughter_v4.vect().theta()*180.0/acos(-1.0);
                double p = daughter_v4.vect().mag();
                double pt = daughter_v4.vect().perp();
                daughter_v4.boost(kinematics.CMBoost);

                if (daughter_pid == 211){  // pi+
                    pip_E_cms_truth.push_back(daughter_v4.e());
                    pip_px_cms_truth.push_back(daughter_v4.px());
                    pip_py_cms_truth.push_back(daughter_v4.py());
                    pip_pz_cms_truth.push_back(daughter_v4.pz());
                    pip_theta.push_back(theta_lab);
                    pip_p_truth.push_back(p);
                    pip_pt_truth.push_back(pt);
                }
                else if (daughter_pid == -211){  // pi-
                    pim_E_cms_truth.push_back(daughter_v4.e());
                    pim_px_cms_truth.push_back(daughter_v4.px());
                    pim_py_cms_truth.push_back(daughter_v4.py());
                    pim_pz_cms_truth.push_back(daughter_v4.pz());
                    pim_theta.push_back(theta_lab);
                    pim_p_truth.push_back(p);
                    pim_pt_truth.push_back(pt);
                }
            }
        }

        HepLorentzVector v4tmp (gen_it->PX(),gen_it->PY(),gen_it->PZ(),gen_it->E()); 
        v4tmp.boost(kinematics.CMBoost);

        if (status == 1) {  // stable final state particles
            if (abs(pid) == 12 || abs(pid) == 14 || abs(pid) == 16) // exclude neutrinos
                continue;

            //if (pid == 130) // K_L0
            //    continue;

            //if (gen_it->E() > 0.1) { // same cut as in reconstructed level
            allParticlesCMS_truth.push_back(v4tmp.vect());
            //}
        }

        if (abs(pid) == 1 || abs(pid) == 2 || abs(pid) == 3 || abs(pid) == 4) {
            Nquark++;
            if(Nquark == 1){ // primary quark would in first 
                q_momentum = v4tmp;
            }
            //int particle_index = std::distance(gen_hep_Mgr.begin(), gen_it);
        }
        
    }
    
    Hep3Vector t_truth = thrust(allParticlesCMS_truth.begin(), allParticlesCMS_truth.end(), retSelf);

    double thrust_vec[3] = { t_truth.mag(), t_truth.z()/t_truth.mag(), t_truth.phi() };
    for(int i=0; i<3; i++){
        thrust_truth.push_back(thrust_vec[i]);
    }

    // qqbar axis
    double qqbar_vec[2] = { 
        q_momentum.z() / q_momentum.vect().mag(), 
        q_momentum.phi() 
    };
    for(int i=0; i<2; i++){
        qqbar_axis.push_back(qqbar_vec[i]);
    }

    n_ks_truth = pip_E_cms_truth.size();
    if (n_ks_truth < 1) return; // ensure there is at least a K_S in truth level

    mctree->Fill();

    return;
}

void KsSpinAlignment_NullTest::fillMCTruthForTopo() {
    // Initialize
    m_nMCGen = 0;
    for (int i = 0; i < MAX_MC_PARTICLES; i++) {
        m_MCGenPDG[i] = std::numeric_limits<double>::quiet_NaN();
        m_MCGenMothIndex[i] = std::numeric_limits<double>::quiet_NaN();
    }

    Gen_hepevt_Manager& gen_mgr = Gen_hepevt_Manager::get_manager();

    // First pass: filter valid generator-level particles and build ID to new index map
    std::vector<const Gen_hepevt*> valid_particles;
    std::map<int, int> id_to_idx;
    int idx = 0;
    for (Gen_hepevt_Manager::iterator it = gen_mgr.begin(); it != gen_mgr.end(); ++it) {
        int pdg = it->idhep();
        int status = it->isthep();
        // Only keep particles with valid PDG and generator status
        if (pdg == 0) continue;
        if (abs(pdg) > 1e7) continue; // skip junk
        if (!(status == 1 || status == 2)) continue; // only final/intermediate
        valid_particles.push_back(&(*it));
        id_to_idx[it->get_ID()] = idx;
        idx++;
        if (idx >= MAX_MC_PARTICLES) break;
    }

    // Second pass: fill arrays
    for (size_t i = 0; i < valid_particles.size(); ++i) {
        const Gen_hepevt* it = valid_particles[i];
        m_MCGenPDG[i] = static_cast<double>(it->idhep());
        // Get mother index
        if (it->mother()) {
            int mother_id = it->mother().get_ID();
            if (id_to_idx.find(mother_id) != id_to_idx.end()) {
                m_MCGenMothIndex[i] = static_cast<double>(id_to_idx[mother_id]);
            } else {
                m_MCGenMothIndex[i] = -1.0;
            }
        } else {
            m_MCGenMothIndex[i] = -1.0;  // No mother
        }
    }
    m_nMCGen = valid_particles.size();
}

void KsSpinAlignment_NullTest::disp_stat(const char*){
    return;
} 


void KsSpinAlignment_NullTest::term(){
    // Print final cut flow statistics
    cout << "\n=== Final Ks Cut Flow Statistics ===" << endl;
    cout << "Total Ks candidates (after FindKs): " << total_ks_candidates << endl;
    cout << "After pt cut:                            " << total_after_pt 
         << " (rejected: " << total_ks_candidates - total_after_pt 
         << ", eff: " << (total_ks_candidates>0 ? 100.0*total_after_pt/total_ks_candidates : 0) << "%)" << endl;
    cout << "After cosTheta cut:                      " << total_after_cosTheta 
         << " (rejected: " << total_after_pt - total_after_cosTheta 
         << ", eff: " << (total_after_pt>0 ? 100.0*total_after_cosTheta/total_after_pt : 0) << "%)" << endl;
    cout << "After lepton veto:                       " << total_after_lepton 
         << " (rejected: " << total_after_cosTheta - total_after_lepton 
         << ", eff: " << (total_after_cosTheta>0 ? 100.0*total_after_lepton/total_after_cosTheta : 0) << "%)" << endl;
    cout << "After pion PID cut:                      " << total_after_pid 
         << " (rejected: " << total_after_lepton - total_after_pid 
         << ", eff: " << (total_after_lepton>0 ? 100.0*total_after_pid/total_after_lepton : 0) << "%)" << endl;
    cout << "Overall efficiency:                      " 
         << (total_ks_candidates>0 ? 100.0*total_after_pid/total_ks_candidates : 0) << "%" << endl;
    cout << "====================================" << endl;
    
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

