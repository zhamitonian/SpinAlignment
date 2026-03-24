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
check the basf2 hadronic slection implementation
comparing to the belle1 hadronic selection
- change from the Ks scripte , so the class name , etc
does not change, nerver mind
-----------------------------------------

version : v1.0.0
Date    : 2026.01.27
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
    kinematics.yAngle = 0.013241569153711791;

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

    tree->Branch("trk_p", &trk_p);
    tree->Branch("trk_costheta", &trk_costheta);
    tree->Branch("trk_phi", &trk_phi);
    tree->Branch("gam_p", &gam_p);
    tree->Branch("gam_costheta", &gam_costheta);
    tree->Branch("gam_phi", &gam_phi);
    tree->Branch("gam2_p", &gam2_p);
    tree->Branch("gam2_costheta", &gam2_costheta);
    tree->Branch("gam2_phi", &gam2_phi);
   
    tree->Branch("trk_p_CMS", &trk_p_CMS);
    tree->Branch("trk_costheta_CMS", &trk_costheta_CMS);
    tree->Branch("trk_phi_CMS", &trk_phi_CMS);
    tree->Branch("gam_p_CMS", &gam_p_CMS);
    tree->Branch("gam_costheta_CMS", &gam_costheta_CMS);
    tree->Branch("gam_phi_CMS", &gam_phi_CMS);
    return;
}


void KsSpinAlignment_NullTest::event(BelleEvent* evptr, int* status){
    (void)evptr; (void)status;

    int  expNo=Belle_event_Manager::get_manager().begin()->ExpNo();
    int  evtNo=Belle_event_Manager::get_manager().begin()->EvtNo();
    int  runNo=Belle_event_Manager::get_manager().begin()->RunNo();
    const HepPoint3D ip_position = IpProfile::e_position();

    Evtcls_flag_Manager& FlagMgr = Evtcls_flag_Manager::get_manager();
    Evtcls_hadronic_flag_Manager& ClassMgr = Evtcls_hadronic_flag_Manager::get_manager();
    int EvtClsFlg = FlagMgr.begin()->flag(0);
    int EvtClsFlgB = ClassMgr.begin()->hadronic_flag(2);
    if (EvtClsFlg < 0 )
    {
        cout << "event " << evtNo << " in run " << runNo << " exp " << expNo << " is not hadron type!" << endl;
        cout << EvtClsFlg << endl;
        return; // only process hadron type events
    }

    if(!EvtClsFlgB){
        cout << "event " << evtNo << " in run " << runNo << " exp " << expNo << " is not hadronB skimed data!" << endl;
        cout << EvtClsFlgB << endl;
        return; // only process hadronB skimed data)
    }

    if(!IpProfile::usable()) {
        cout <<" ip not usable ..." << endl;
        return;
    }
    countEvt ++;

    trk_p.clear(); trk_costheta.clear(); trk_phi.clear();
    gam_p.clear(); gam_costheta.clear(); gam_phi.clear();
    gam2_p.clear(); gam2_costheta.clear(); gam2_phi.clear();

    trk_p_CMS.clear(); trk_costheta_CMS.clear(); trk_phi_CMS.clear();
    gam_p_CMS.clear(); gam_costheta_CMS.clear(); gam_phi_CMS.clear();

    Mdst_charged_Manager& charged_mgr = Mdst_charged_Manager::get_manager();
    Mdst_gamma_Manager& gamma_mgr = Mdst_gamma_Manager::get_manager();

    //---------------------------------------------------------
    Vp3 particles_p3_cms;
    Vp4 particles_p4_cms;
    particles_p3_cms.clear();
    particles_p4_cms.clear();

    double Evis = 0; double Pz = 0; double HeavyJetMass = 0; 
    double Ntrk = 0; double Ncls = 0;

    // ----------------------------------------------------------------------------------
    /*
    Evtcls_hadron_charged_Manager& hadron_chg_mgr = Evtcls_hadron_charged_Manager::get_manager();
    for(Evtcls_hadron_charged_Manager::iterator it = hadron_chg_mgr.begin();
        it != hadron_chg_mgr.end(); it++){

        if (!it->charged()) continue;  // Check if charged is valid

        const Mdst_charged& Charged = it->charged();
        Hep3Vector P(Charged.px(), Charged.py(), Charged.pz());
        
        double E = sqrt(P.mag2() + xmass[2]*xmass[2]);
        HepLorentzVector P_cms(P,E);        
        // lab
        trk_p.push_back(P_cms.vect().mag());
        trk_costheta.push_back(cos(P_cms.vect().theta()));
        trk_phi.push_back(P_cms.vect().phi());
        // cms
        P_cms.boost(kinematics.CMBoost);
        Hep3Vector vec = P_cms.vect();
        vec.rotateY(-kinematics.yAngle);
        P_cms = HepLorentzVector(vec, P_cms.e());

        particles_p3_cms.push_back(P_cms.vect());
        particles_p4_cms.push_back(P_cms);

        Evis += P_cms.e();
        Pz += P_cms.pz();
    }
   
    Ntrk = particles_p3_cms.size();
    */

    /*  in evtcls_hadron_neutral , it actually cut E > 0.105 , not 0.1
    // 对 gamma 也类似
    Evtcls_hadron_neutral_Manager& hadron_neu_mgr = Evtcls_hadron_neutral_Manager::get_manager();
    
    for(Evtcls_hadron_neutral_Manager::iterator it = hadron_neu_mgr.begin();
        it != hadron_neu_mgr.end(); it++){
        
        if (!it->gamma()) continue;  // Check if gamma is valid
        
        const Mdst_gamma& Gamma = it->gamma();
        Hep3Vector P(Gamma.px(), Gamma.py(), Gamma.pz());
        
        double E = P.mag();
        HepLorentzVector vec_p4(P, E);
        // lab
        gam_p.push_back(vec_p4.vect().mag());
        gam_costheta.push_back(cos(vec_p4.vect().theta()));
        gam_phi.push_back(vec_p4.vect().phi());
        // cms
        vec_p4.boost(kinematics.CMBoost);
        Hep3Vector vec = vec_p4.vect();
        vec.rotateY(-kinematics.yAngle);
        vec_p4 = HepLorentzVector(vec, vec_p4.e());
        
        particles_p3_cms.push_back(vec_p4.vect());
        particles_p4_cms.push_back(vec_p4);

        Evis += vec_p4.e();
    }

    */

    for(Mdst_charged_Manager::iterator it=charged_mgr.begin(); it!=charged_mgr.end(); it++){
        const Mdst_charged& Charged = *it;

        Hep3Vector P(Charged.px(), Charged.py(), Charged.pz());
        double dr, dz, refitPx, refitPy, refitPz;
        getDrDz(it, 0, dr, dz, refitPx, refitPy, refitPz);

        if(P.perp() > 0.1 && fabs(dr) < 2.0 && fabs(dz) < 4.0){

            double E = sqrt(P.mag2() + xmass[2]*xmass[2]);
            HepLorentzVector P_cms(P,E);  
            // lab
            trk_p.push_back(P_cms.vect().mag());
            trk_costheta.push_back(cos(P_cms.vect().theta()));
            trk_phi.push_back(P_cms.vect().phi());
            // cms
            P_cms.boost(kinematics.CMBoost);
            Hep3Vector vec = P_cms.vect();
            //vec.rotateY(-kinematics.yAngle);
            P_cms = HepLorentzVector(vec, P_cms.e());

            trk_p_CMS.push_back(P_cms.vect().mag());
            trk_costheta_CMS.push_back(cos(P_cms.vect().theta()));
            trk_phi_CMS.push_back(P_cms.vect().phi());

            particles_p3_cms.push_back(P_cms.vect());
            particles_p4_cms.push_back(P_cms);

            Evis += P_cms.e();
            Pz += P_cms.pz();
        }
    }

    Ntrk = particles_p3_cms.size();

    for(Mdst_gamma_Manager::iterator gam = gamma_mgr.begin(); gam != gamma_mgr.end(); gam++){
        const Mdst_gamma &gamma = *gam;
        Mdst_ecl &ecl = gamma.ecl();
        Hep3Vector P(gamma.p(0), gamma.p(1), gamma.p(2));

        double E = P.mag();

        if (ecl.quality() ==0 ){
            if ( P.mag() >= 0.1) {         // E cut
        HepLorentzVector vec_p4(P, E);
        // lab
        gam_p.push_back(vec_p4.vect().mag());
        gam_costheta.push_back(cos(vec_p4.vect().theta()));
        gam_phi.push_back(vec_p4.vect().phi());
        // cms
        vec_p4.boost(kinematics.CMBoost);
        Hep3Vector vec = vec_p4.vect();
        //vec.rotateY(-kinematics.yAngle);
        vec_p4 = HepLorentzVector(vec, vec_p4.e());

        gam_p_CMS.push_back(vec_p4.vect().mag());
        gam_costheta_CMS.push_back(cos(vec_p4.vect().theta()));
        gam_phi_CMS.push_back(vec_p4.vect().phi());
        
        particles_p3_cms.push_back(vec_p4.vect());
        particles_p4_cms.push_back(vec_p4);

        Evis += vec_p4.e();
            }
        }    
    }

    Ncls = particles_p3_cms.size() - Ntrk;

    //Hep3Vector thr = thrust(particles_p3_cms.begin(), particles_p3_cms.end(), retSelf);
    Hep3Vector thr = thrust(particles_p3_cms.begin(), particles_p3_cms.end(), SelfFunc(Hep3Vector()));
    Hep3Vector tn  = thr.unit();
    double ThrParam = thr.mag();
    HeavyJetMass = calculateHeavyJetMassEnergy(particles_p4_cms, tn).first;

    // ----------------------------------------------------------------------------------
    // Additional hadron event selection (for HadronB skimed data)
    Evtcls_hadron_info_Manager& hadronInfo_mgr = Evtcls_hadron_info_Manager::get_manager();
    //double Evis = hadronInfo_mgr.begin()->Evis();
    double Evis_old = hadronInfo_mgr.begin()->Evis();
    double Esum = hadronInfo_mgr.begin()->Esum();
    double Psum = hadronInfo_mgr.begin()->Psum();
    double Pz_old   = hadronInfo_mgr.begin()->Pz();
    double R2   = hadronInfo_mgr.begin()->R2();
    double HeavyJetMass_old = hadronInfo_mgr.begin()->HeavyJetMass();
    int Ntrk_old = hadronInfo_mgr.begin()->Ntrk();
    int Ncls_old = hadronInfo_mgr.begin()->Ncls();
    double Thrust = hadronInfo_mgr.begin()->Thrust();
    if(Evis < cuts::Evis)
        return;

    if(Ncls !=  Ncls_old)
    {
      cout << "Ncls not match! calc: " << Ncls << " evtcls: " << Ncls_old << endl;  
    }
    if (Ntrk != Ntrk_old)
    {
      cout << "Ntrk not match! calc: " << Ntrk << " evtcls: " << Ntrk_old << endl;  
    }

    // ----------------------------------------------------------------------------------
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

