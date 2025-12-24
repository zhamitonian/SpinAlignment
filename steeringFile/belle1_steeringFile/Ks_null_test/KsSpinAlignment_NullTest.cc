#include "belle.h"
#include <cmath>
#include <algorithm>

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

version : v1.0.0
Date    : 2025.12.24
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

    if(!(countEvt%10000)) 
    cout << "evt " <<countEvt<< " expNo "<< expNo << ", runNo "<< runNo << ", evtNo "<< evtNo <<endl;
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
    

    // start our selection
    Vp3 hadrons_p3_cms;
    Vp4 hadrons_p4_cms;
    hadrons_p3_cms.clear();
    hadrons_p4_cms.clear();
    pip_E_cms.clear();
    pip_px_cms.clear();
    pip_py_cms.clear();
    pip_pz_cms.clear();
    pim_E_cms.clear();
    pim_px_cms.clear();
    pim_py_cms.clear();
    pim_pz_cms.clear();

    pip_p.clear();
    pim_p.clear();

    if (isMC){
        pip_isSignal.clear();
        pim_isSignal.clear();

        pip_E_cms_gen.clear();
        pip_px_cms_gen.clear();
        pip_py_cms_gen.clear();
        pip_pz_cms_gen.clear();
        pim_E_cms_gen.clear();
        pim_px_cms_gen.clear();
        pim_py_cms_gen.clear();
        pim_pz_cms_gen.clear();

    }

    for(Mdst_charged_Manager::iterator ch_it=charged_mgr.begin(); ch_it!=charged_mgr.end(); ch_it++){
        const Mdst_charged &chg = *ch_it;

        Hep3Vector vec_p3(chg.p(0), chg.p(1), chg.p(2));

        double dr, dz, refitPx, refitPy, refitPz;
        getDrDz(ch_it, 0, dr, dz, refitPx, refitPy, refitPz);

        if(fabs(dr) > cuts::vertexR || fabs(dz) > cuts::vertexZ) 
            continue;
        double pt = vec_p3.perp();
        if (pt < cuts::trkPt) 
            continue;
        //if (vec_p3.mag() < cuts::trkP)
        //    continue;
        double cosTheta = cos(vec_p3.theta());
        if (cosTheta < cuts::trk_cosTheta_min || cosTheta > cuts::trk_cosTheta_max)
            continue;

        // identification
        bool isLepton = false;
        int charge = chg.charge();
        int massHyp = 2;

        eid sel_e(*ch_it);  
        atc_pid selKPi(3,1,5,3,2);  //K/pi separation
        atc_pid selPiP(3,1,5,2,4); 
        atc_pid selKP(3,1,5,3,4);
        double mu_id = 0;
        Muid_mdst muID(*ch_it);
        if (muID.Chi_2() >0) 
            mu_id = muID.Muon_likelihood();

        double atcKPi=selKPi.prob(*ch_it);
        double atcKP=selKP.prob(*ch_it);
        double atcPiP=selPiP.prob(*ch_it);

        float e_cut=0.85;
        float mu_cut=0.9;
        double e_id=sel_e.prob(3,-1,5);

        if(mu_id > mu_cut) {
            massHyp = 1; // muon
            isLepton = true;
        }
        else if(e_id > e_cut && mu_id < mu_cut) {
            massHyp = 0; // electron
            isLepton = true;
        }

        if(!isLepton) {
            if(atcKPi > 0.6 && atcKP > 0.2) {
                massHyp = 3; // kaon
            }
            else if(atcKPi < 0.6 && atcPiP > 0.2) {
                massHyp = 2; // pion
            }
            else if(atcKP < 0.2 && atcPiP < 0.2) {
                massHyp = 4; // proton
            }
        }

        // ---------------------------------------------------------
        // match MC truth and retrieve mother particle info
        if (isMC == 1 && massHyp == 2) {  // Check pions (massHyp == 2)
            const Gen_hepevt &gen = gen_level(get_hepevt(chg));
            Gen_hepevt &mother = gen.mother();
            int mc_mother_pdg = mother.idhep();
            int mc_pdg = gen.idhep();

            HepLorentzVector p4_truth(gen.PX(), gen.PY(), gen.PZ(), gen.E());
            p4_truth.boost(kinematics.CMBoost);

            bool isPion = false;
            bool fromKs = false;

            if (abs(mc_mother_pdg) == 310)  // K_S (PDG = 310)
                fromKs = true;

            if (mc_pdg * charge == 211)  // pi± (PDG = ±211)
                isPion = true;

            if (charge == 1) {
                if (isPion && fromKs) {
                    pip_isSignal.push_back(true);
                    pip_E_cms_gen.push_back(p4_truth.e());
                    pip_px_cms_gen.push_back(p4_truth.px());
                    pip_py_cms_gen.push_back(p4_truth.py());
                    pip_pz_cms_gen.push_back(p4_truth.pz());
                } else {
                    pip_isSignal.push_back(false);
                }
            }
            else if(charge == -1)  {
                if (isPion && fromKs) {
                    pim_isSignal.push_back(true);
                    pim_E_cms_gen.push_back(p4_truth.e());
                    pim_px_cms_gen.push_back(p4_truth.px());
                    pim_py_cms_gen.push_back(p4_truth.py());
                    pim_pz_cms_gen.push_back(p4_truth.pz());
                } else {
                    pim_isSignal.push_back(false);
                }
            }
        }
        
        // ---------------------------------------------------------
        // store track kinematic variables
        double E = sqrt(vec_p3.mag2() + xmass[massHyp]*xmass[massHyp]); 
        HepLorentzVector p4_boosted(vec_p3, E);
        p4_boosted.boost(kinematics.CMBoost);

        if (massHyp == 2){  // Store pions (massHyp == 2)
            if(charge == 1) {
                pip_E_cms.push_back(p4_boosted.e());
                pip_px_cms.push_back(p4_boosted.px());
                pip_py_cms.push_back(p4_boosted.py());
                pip_pz_cms.push_back(p4_boosted.pz());
                pip_p.push_back(vec_p3.mag());
            }
            else if(charge == -1)  {
                pim_E_cms.push_back(p4_boosted.e());
                pim_px_cms.push_back(p4_boosted.px());
                pim_py_cms.push_back(p4_boosted.py());
                pim_pz_cms.push_back(p4_boosted.pz());
                pim_p.push_back(vec_p3.mag());
            }
        }

        hadrons_p3_cms.push_back(p4_boosted.vect());
        hadrons_p4_cms.push_back(p4_boosted);

    }

    // ensure there has at least a pi+ pi- pair
    if (pip_px_cms.size() * pim_px_cms.size() ==0) return;  
    
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
            //    allParticlesCMS_truth.push_back(v4tmp.vect());
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