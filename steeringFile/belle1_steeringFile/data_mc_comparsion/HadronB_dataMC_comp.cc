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

#include "./HadronB_dataMC_comp.h"
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
#include EVTVTX_H
#include HEPEVT_H

using namespace std;

/*
-----------------------------------------
hadronB data vs MC comparison

version : v1.0.1
Date    : 2025.12.04
Author  : Zhen Wang
*/

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

HadronB_dataMC_comp::HadronB_dataMC_comp(){
    output_filename=new char[256];
    return;
}

void HadronB_dataMC_comp::init(int *status){
    cout << "in HadronB_dataMC_comp::init()" << endl;
    hist_def();
    *status=0;
    return;
}

void HadronB_dataMC_comp::begin_run(BelleEvent* evptr, int* status){
    (void)evptr; (void)status;

    IpProfile::begin_run();
    BeamEnergy::begin_run();

    
    kinematics.ler_e = BeamEnergy::E_LER();
    kinematics.her_e = BeamEnergy::E_HER();
    kinematics.x_angle = BeamEnergy::Cross_angle();

    if(kinematics.ler_e < 3.0 || kinematics.her_e < 7.0 || kinematics.ler_e > 5.0 || kinematics.her_e > 9.0){
        cout<<"someting is wrong ler_e is "<< kinematics.ler_e << " her_e is "<< kinematics.her_e <<endl;
        return;
    }

    kinematics.cm = HepLorentzVector(-kinematics.her_e * sin(kinematics.x_angle), 0., -kinematics.her_e * cos(kinematics.x_angle) + kinematics.ler_e, kinematics.her_e + kinematics.ler_e);
    kinematics.CMBoost = kinematics.cm.boostVector();
    kinematics.sqrts = kinematics.cm.m();

    return;
}


void HadronB_dataMC_comp::end_run(BelleEvent* evptr, int* status){
    (void)evptr; (void)status;
    return; 
}


void HadronB_dataMC_comp::hist_def(){
    output_file = new TFile(output_filename, "RECREATE");
    
    // Create histograms for event variables
    h_Evis = new TH1F("h_Evis", "Visible Energy;E_{vis} (GeV);[MeV]", 100, 1, 13);
    h_HeavyJetMass_over_Evis = new TH1F("h_HeavyJetMass_over_Evis", "Heavy Jet Mass / E_{vis};M_{jet}/E_{vis};[]", 60, 0.1, 0.7);
    h_Esum = new TH1F("h_Esum", "Sum Energy;E_{sum} (GeV);[MeV]", 100, 0, 10);
    h_Psum = new TH1F("h_Psum", "Sum Momentum;P_{sum} (GeV/c);[MeV]", 100, 0, 12);
    h_Pz = new TH1F("h_Pz", "Z Momentum;P_{z} (GeV/c);[MeV]", 100, -5, 5);
    h_R2 = new TH1F("h_R2", "R2;R2;[]", 100, 0, 1);
    h_HeavyJetMass = new TH1F("h_HeavyJetMass", "Heavy Jet Mass;M_{jet} (GeV/c^{2});[MeV]", 78, 0.5, 6);
    h_Ntrk = new TH1F("h_Ntrk", "Number of Tracks;N_{trk};[]", 20, 0, 20);
    h_Ncls = new TH1F("h_Ncls", "Number of Clusters;N_{cls};[]", 25, 0, 25);
    h_thrust = new TH1F("h_thrust", "Thrust;Thrust;[]", 100, 0.5, 1.0);
    h_thrust_costheta = new TH1F("h_thrust_costheta", "Thrust cos#theta;cos#theta;[]", 100, -1, 1);
    h_thrust_phi = new TH1F("h_thrust_phi", "Thrust #phi;#phi (rad);[rad]", 100, -3.15, 3.15);
    h_HeavyJetEnergy = new TH1F("h_HeavyJetEnergy", "Heavy Jet Energy;E_{jet} (GeV);[MeV]", 75, 1.5, 7.5);
    h_Sphericity = new TH1F("h_Sphericity", "Sphericity;S;[]", 100, 0, 1);
    h_Aplanarity = new TH1F("h_Aplanarity", "Aplanarity;A;[]", 100, 0, 0.5);
    
    // Create histograms for all charged tracks (lab frame)
    h_trk_dr = new TH1F("h_trk_dr", "Track dr;dr (cm);[cm]", 100, -1, 1);
    h_trk_dz = new TH1F("h_trk_dz", "Track dz;dz (cm);[cm]", 100, -4, 4);
    h_trk_p = new TH1F("h_trk_p", "Track Momentum;p (GeV/c);[MeV]", 100, 0, 5);
    h_trk_costheta = new TH1F("h_trk_costheta", "Track cos#theta;cos#theta;[]", 100, -1, 1);
    h_trk_phi = new TH1F("h_trk_phi", "Track #phi;#phi (rad);[rad]", 100, -3.15, 3.15);
    
    // Create histograms for photons
    h_photon_p = new TH1F("h_photon_p", "Photon Momentum;p (GeV/c);[MeV]", 50, 0, 2.5);
    h_photon_costheta = new TH1F("h_photon_costheta", "Photon cos#theta;cos#theta;[]", 40, -1, 1);
    h_photon_phi = new TH1F("h_photon_phi", "Photon #phi;#phi (rad);[rad]", 40, -3.15, 3.15);
    
    // Create histograms for primary vertex
    h_PrimeR = new TH1F("h_PrimeR", "Primary Vertex R;R (cm);[cm]", 40, 0, 0.2);
    h_PrimeZ = new TH1F("h_PrimeZ", "Primary Vertex Z;Z (cm);[cm]", 100, -2, 2);
}


void HadronB_dataMC_comp::event(BelleEvent* evptr, int* status){
    (void)evptr; (void)status;

    int  expNo=Belle_event_Manager::get_manager().begin()->ExpNo();
    int  evtNo=Belle_event_Manager::get_manager().begin()->EvtNo();
    int  runNo=Belle_event_Manager::get_manager().begin()->RunNo();

    Evtcls_hadronic_flag_Manager& ClassMgr = Evtcls_hadronic_flag_Manager::get_manager();

    int EvtClsFlgB = ClassMgr.begin()->hadronic_flag(2);

    if(!EvtClsFlgB){
        return; // only process hadronB skimed data)
    }

    if(!(countEvt%10000)) 
    cout << "evt " <<countEvt<< " expNo "<< expNo << ", runNo "<< runNo << ", evtNo "<< evtNo <<endl;
    if(!IpProfile::usable()) {
        cout <<" ip not usable ..." << endl;
        return;
    }
    countEvt ++;

    // Get primary vertex position (same method as EventClassify.cc)
    HepPoint3D PrimaryVtx;
    int VtxQuality = -99;
    
    Evtvtx_primary_vertex_Manager& VtxMgr = Evtvtx_primary_vertex_Manager::get_manager();
    std::vector<Evtvtx_primary_vertex>::iterator vtx_it = VtxMgr.begin();
    
    if(vtx_it != VtxMgr.end() && *vtx_it) {
        VtxQuality = (*vtx_it).quality();
        if(VtxQuality >= 2) {
            PrimaryVtx = HepPoint3D((*vtx_it).PV_x(), (*vtx_it).PV_y(), (*vtx_it).PV_z());
        } else {
            PrimaryVtx = IpProfile::position(1);  // fallback to IP profile
        }
    } else {
        PrimaryVtx = IpProfile::position(1);  // fallback to IP profile
    }
    
    PrimeR = sqrt(PrimaryVtx.x() * PrimaryVtx.x() + PrimaryVtx.y() * PrimaryVtx.y());
    PrimeZ = PrimaryVtx.z();

    //---------------------------------------------------------
    Vp3 particles_p3_cms;
    Vp4 particles_p4_cms;
    particles_p3_cms.clear();
    particles_p4_cms.clear();
    
    trk_dr.clear();
    trk_dz.clear();
    trk_p.clear();
    trk_costheta.clear();
    trk_phi.clear();
    trk_charge.clear();
    
    photon_p.clear();
    photon_costheta.clear();
    photon_phi.clear();

Evtcls_hadron_charged_Manager& hadron_chg_mgr = Evtcls_hadron_charged_Manager::get_manager();
for(Evtcls_hadron_charged_Manager::iterator it = hadron_chg_mgr.begin();
    it != hadron_chg_mgr.end(); it++){

    if (!it->charged()) continue;  // Check if charged is valid

    const Mdst_charged& Charged = it->charged();
    Hep3Vector P(Charged.px(), Charged.py(), Charged.pz());
    
    // Store track variables in lab frame (no cuts applied)
    double p = P.mag();
    double cosTheta = cos(P.theta());
    double phi = P.phi();
    int charge = Charged.charge();
    
    // Get dr, dz directly from helix parameters (matching EventClassify.cc implementation)
    Mdst_trk& Trk = Charged.trk();
    Mdst_trk_fit& Trkfit = Trk.mhyp(2); // pion hypothesis
    double dr = Trkfit.helix(0);
    double dz = Trkfit.helix(3);
    
    trk_dr.push_back(dr);
    trk_dz.push_back(dz);
    trk_p.push_back(p);
    trk_costheta.push_back(cosTheta);
    trk_phi.push_back(phi);
    trk_charge.push_back(charge);
    
    // Fill track histograms
    h_trk_dr->Fill(dr);
    h_trk_dz->Fill(dz);
    h_trk_p->Fill(p);
    h_trk_costheta->Fill(cosTheta);
    h_trk_phi->Fill(phi);
    
    // For CMS momentum calculation
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
    
    // Store photon variables in lab frame
    double p_photon = P.mag();
    double cosTheta_photon = cos(P.theta());
    double phi_photon = P.phi();
    photon_p.push_back(p_photon);
    photon_costheta.push_back(cosTheta_photon);
    photon_phi.push_back(phi_photon);
    
    // Fill photon histograms
    h_photon_p->Fill(p_photon);
    h_photon_costheta->Fill(cosTheta_photon);
    h_photon_phi->Fill(phi_photon);
    
    double E = P.mag();
    HepLorentzVector vec_p4(P, E);
    vec_p4.boost(kinematics.CMBoost);
    
    particles_p3_cms.push_back(vec_p4.vect());
    particles_p4_cms.push_back(vec_p4);
}

    Hep3Vector thr = thrust(particles_p3_cms.begin(), particles_p3_cms.end(), SelfFunc(Hep3Vector()));
    Hep3Vector tn  = thr.unit();
    double ThrParam = thr.mag();
    
    // Calculate Heavy Jet Mass and Energy
    std::pair<double, double> jetResult = calculateHeavyJetMassEnergy(particles_p4_cms, thr);
    double heavyJetMass = jetResult.first;
    double heavyJetEnergy = jetResult.second;
    
    // Calculate Sphericity and Aplanarity
    std::pair<double, double> shapeResult = calculateSphericityAplanarity(particles_p4_cms);
    double sphericity = shapeResult.first;
    double aplanarity = shapeResult.second;

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
    
    // Validation: Compare calculated values with read values
    // Typical differences: Thrust < 0.001, HeavyJetMass < 0.02 (acceptable for physics analysis)
    
    // Fill histograms before Evis cut
    h_Evis->Fill(Evis);
    if(Evis < cuts::Evis)
        return;
    // ----------------------------------------------------------------------------------
    
    // Fill histograms for event variables
    h_Esum->Fill(Esum);
    h_Psum->Fill(Psum);
    h_Pz->Fill(Pz);
    h_R2->Fill(R2);
    h_HeavyJetMass->Fill(HeavyJetMass);
    h_Ntrk->Fill(Ntrk);
    h_Ncls->Fill(Ncls);
    h_thrust->Fill(thr.mag());
    h_thrust_costheta->Fill(thr.z()/thr.mag());
    h_thrust_phi->Fill(thr.phi());
    h_HeavyJetEnergy->Fill(heavyJetEnergy);
    h_Sphericity->Fill(sphericity);
    h_Aplanarity->Fill(aplanarity);
    h_PrimeR->Fill(PrimeR);
    h_PrimeZ->Fill(PrimeZ);
    if(Evis > 0) {
        h_HeavyJetMass_over_Evis->Fill(HeavyJetMass / Evis);
    }

    return;
}


void HadronB_dataMC_comp::disp_stat(const char*){
    return;
} 


void HadronB_dataMC_comp::term(){
    output_file->Write();
    output_file->Close();
    return;
}


std::pair<double, double> HadronB_dataMC_comp::calculateHeavyJetMassEnergy(const Vp4& particles, const Hep3Vector& thrustAxis) {
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


std::pair<double, double> HadronB_dataMC_comp::calculateSphericityAplanarity(const std::vector<HepLorentzVector>& particles) {
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

void HadronB_dataMC_comp::other(int* , BelleEvent*, int* ){
    return;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif


// version log 
// v1.0.0 : 
//  - First version for hadronB data vs MC comparison
// v1.0.1 :
//  - modify binning and range of histograms