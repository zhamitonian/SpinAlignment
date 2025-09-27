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

#include "./SpinAlignment.h"
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

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

Hep3Vector& retSelf(Hep3Vector& vec)
{
	return vec;
};


SpinAlignment::SpinAlignment(){
    output_filename=new char[256];

    r1.SetSeed(100);
    r5.SetSeed(500);
    r6.SetSeed(600);
    r7.SetSeed(700);
    r8.SetSeed(800);
    r9.SetSeed(900);

    return;
}


void SpinAlignment::init(int *status){
    cout << "in SpinAlignment::init()" << endl;
    hist_def();
    *status=0;
    return;
}


void SpinAlignment::begin_run(BelleEvent* evptr, int* status){
    // loading calibration constants and lookup tables needed for particle ID ?
    eid::init_data();
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
    kinematics.firstElectronCM = HepLorentzVector(kinematics.her_e * sin(kinematics.x_angle), 0., kinematics.her_e * cos(kinematics.x_angle), kinematics.her_e);
    kinematics.secondElectronCM = HepLorentzVector(0., 0., -kinematics.ler_e, kinematics.ler_e);

    return;
}


void SpinAlignment::end_run(BelleEvent* evptr, int* status){
    (void)evptr; (void)status;
    return; 
}


void SpinAlignment::hist_def(){
    output_file = new TFile(output_filename, "RECREATE");
    tree = new TTree("event", "event");

    // macro to add branches
    #define ADDBRANCH__(name, var, type)   tree->Branch(#name, &m_info.var, #type)
    #define ADDBRANCH(x, type)            ADDBRANCH__(x, x, x/type)
    #define ADDBARRAY__(name, var, n, type) ADDBRANCH__(name, var, name[n]/type)
    #define ADDBARRAY(x, n, type)          ADDBARRAY__(x, x, n, type)

    ADDBRANCH(evtNo, I);
    ADDBRANCH(runNo, I);
    ADDBRANCH(expNo, I);
    ADDBRANCH(Q, D);
    ADDBRANCH(Ecms, D);
    ADDBRANCH(Energy_cms, D);
    ADDBRANCH(Evis_cms, D);
    ADDBRANCH(ECLEnergyWO, D);
    ADDBRANCH(ECLEnergy, D);
    ADDBRANCH(BalancePz_cms, D);
    ADDBRANCH(HeavyJetMass, D);
    ADDBRANCH(HeavyJetEnergy, D);
    ADDBRANCH(sphericity, D);
    ADDBRANCH(aplanarity, D);
    ADDBRANCH(nGood, I);
    ADDBRANCH(nPip, I);
    ADDBRANCH(nPim, I);
    ADDBRANCH(nCluster, I);
    ADDBRANCH(nPhoton, I);
    ADDBRANCH(z, D);
    ADDBRANCH(pt, D);
    ADDBARRAY(foxWolfram, 5, D);
    ADDBARRAY(cms_vecP, 4, D);
    ADDBARRAY(thrust, 3, D);

    // vector<double> type branches
    tree->Branch("photon_p", &photon_p);
    tree->Branch("photon_theta", &photon_theta);
    tree->Branch("photon_phi", &photon_phi);
    tree->Branch("pip_p", &pip_p);
    tree->Branch("pim_p", &pim_p);
    tree->Branch("pip_theta", &pip_theta);
    tree->Branch("pim_theta", &pim_theta);
    tree->Branch("pip_phi", &pip_phi);
    tree->Branch("pim_phi", &pim_phi);

    return;
}


void SpinAlignment::event(BelleEvent* evptr, int* status){
    (void)evptr; (void)status;

    int  expNo=Belle_event_Manager::get_manager().begin()->ExpNo();
    int  evtNo=Belle_event_Manager::get_manager().begin()->EvtNo();
    int  runNo=Belle_event_Manager::get_manager().begin()->RunNo();

    if(!(countEvt%100)) 
    cout << "evt " <<countEvt<< " expNo "<< expNo << ", runNo "<< runNo << ", evtNo "<< evtNo <<endl;

    if(!IpProfile::usable()) {
        cout <<" ip not usable ..." << endl;
        return;
    }
    const HepPoint3D ip_position = IpProfile::e_position();
    const HepSymMatrix ip_osition_err = IpProfile::e_position_err();
    countEvt ++;

    Vp3 allParticles;
    Vp3 allParticles_Boosted;
    Vp4 Vp4s_Boosted;
    Vp4 Vp4s;
    allParticles.clear();
    allParticles_Boosted.clear();
    photon_p.clear();
    photon_theta.clear();
    photon_phi.clear();
    pip_p.clear();
    pim_p.clear();
    pip_theta.clear();
    pim_theta.clear();
    pip_phi.clear();
    pim_phi.clear();

    double Energy_cms = 0;
    double Evis_cms = 0;
    double ECLEnergyWO= 0;
    double ECLEnergy = 0;
    double BalancePz_cms = 0;
    double heavyJetMass = 0;
    double heavyJetEnergy = 0;
    double sphericity = 0;
    double aplanarity = 0;
    int nPhoton = 0;
    int nCluster= 0;
    int nPip = 0;
    int nPim = 0;

    Mdst_charged_Manager& charged_mgr = Mdst_charged_Manager::get_manager();
    //Mdst_trk_Manager& trk_mgr = Mdst_trk_Manager::get_manager();  not used
    Mdst_gamma_Manager& gamma_mgr = Mdst_gamma_Manager::get_manager();
    Mdst_ecl_aux_Manager& eclaux_mgr = Mdst_ecl_aux_Manager::get_manager();

    char particleName[10];
    strcpy(particleName, "unknown");

    for(Mdst_charged_Manager::iterator ch_it=charged_mgr.begin(); ch_it!=charged_mgr.end(); ch_it++){
        const Mdst_charged &chg = *ch_it;
        //const Mdst_ecl &ecl= chg.ecl();

        int trkID = (*ch_it).get_ID();

        double dr, dz, refitPx, refitPy, refitPz;
        getDrDz(ch_it, 0, dr, dz, refitPx, refitPy, refitPz);

        if(fabs(dr) > cuts::vertexR || fabs(dz) > cuts::vertexZ) 
            continue;

        //Particle p(chg, particleName);
        Hep3Vector vec_p3(chg.p(0), chg.p(1), chg.p(2));
        double E = sqrt(vec_p3.mag2() + xmass[2]*xmass[2]); 

        double pt = vec_p3.perp();
        if (pt < cuts::minPt) 
            continue;

        HepLorentzVector vec_p4(vec_p3, E);
        HepLorentzVector vec_p4_boosted(vec_p3, E);
        vec_p4_boosted.boost(kinematics.CMBoost);

        int charge = chg.charge();
        if(charge == 1) {
            strcpy(particleName, "pi+");
            nPip ++;
            pip_p.push_back(vec_p3.mag());
            pip_theta.push_back(vec_p3.theta());
            pip_phi.push_back(vec_p3.phi());
        }
        else if(charge == -1)  {
            strcpy(particleName, "pi-");
            nPim ++;
            pim_p.push_back(vec_p3.mag());
            pim_theta.push_back(vec_p3.theta());
            pim_phi.push_back(vec_p3.phi());
            //pim_vecP.push_back({vec_p4.px(), vec_p4.py(), vec_p4.pz(), vec_p4.e()});
        }

        allParticles.push_back(vec_p3);
        allParticles_Boosted.push_back(vec_p4_boosted.vect());
        Vp4s_Boosted.push_back(vec_p4_boosted);
        Vp4s.push_back(vec_p4);

        Energy_cms +=  vec_p4_boosted.e();
        Evis_cms += vec_p4_boosted.e();
        BalancePz_cms += vec_p4_boosted.pz();
        //ECLEnergy += ecl.e();
        //ECLEnergyWO += ecl.e();

    }

        int nGood = allParticles_Boosted.size(); 
        if (nGood < 3) return;

    for(Mdst_gamma_Manager::iterator gam = gamma_mgr.begin(); gam != gamma_mgr.end(); gam++){
        const Mdst_gamma &gamma = *gam;
        Mdst_ecl &ecl = gamma.ecl();

        Hep3Vector vec_p3(gamma.p(0), gamma.p(1), gamma.p(2));

        //Mdst_ecl_aux &aux = eclaux_mgr(Panther_ID(gamma.ecl().get_ID()));

        //double e9oe25 = aux.e9oe25();
        double gammaE = vec_p3.mag();
        double ecl_E = ecl.energy();
        photon_p.push_back(gammaE);
        photon_theta.push_back(vec_p3.theta());
        photon_phi.push_back(vec_p3.phi());
        //cout << gammaE - ecl_E << endl;
        
        if(gammaE < cuts::minPhotonE) 
            continue;

        double gam_theta = vec_p3.theta() * 180/ PI; 
        bool thetaInCDCAcceptance = false;
        if(gam_theta > 17.0  && gam_theta < 150.0 ){
            thetaInCDCAcceptance = true;
        } 
        else{
            cout << "a photon out of CDC acceptance with theta "<< gam_theta << endl;
        }

        HepLorentzVector vec_p4(vec_p3, gammaE);
        HepLorentzVector vec_p4_boosted(vec_p3, gammaE);
        vec_p4_boosted.boost(kinematics.CMBoost);

        Energy_cms += vec_p4_boosted.e();
        ECLEnergyWO += vec_p4.e();
        if (thetaInCDCAcceptance){
            Evis_cms += vec_p4_boosted.e();
            ECLEnergy += vec_p4.e();
            BalancePz_cms += vec_p4_boosted.pz();
            nPhoton ++;
        }
        nCluster ++;
        allParticles.push_back(vec_p3);
        allParticles_Boosted.push_back(vec_p4_boosted.vect());
        Vp4s_Boosted.push_back(vec_p4_boosted);
        Vp4s.push_back(vec_p4);

    }
    if (nCluster - nPhoton > 1)
    cout << Evis_cms << " " << Energy_cms << " " << nPhoton << " "<< nCluster << endl;

    Hep3Vector t_cms = thrust(allParticles_Boosted.begin(), allParticles_Boosted.end(), retSelf);
    Hep3Vector t = thrust(allParticles.begin(), allParticles.end(), retSelf);
    FoxWolfram fw = foxwolfram(allParticles_Boosted.begin(), allParticles_Boosted.end(), retSelf);
    FoxWolfram fw1 = foxwolfram(allParticles.begin(), allParticles.end(), retSelf);

    HepLorentzVector total_p4 = HepLorentzVector(0,0,0,0);
    for (size_t i = 0; i < Vp4s.size(); i++) {
        total_p4 += Vp4s[i];
    }
    double Ecms = total_p4.m();

    std::pair<double, double> jetMassEnergy = calculateHeavyJetMassEnergy(Vp4s_Boosted, t_cms.unit());
    heavyJetMass = jetMassEnergy.first;
    heavyJetEnergy = jetMassEnergy.second;
    std::pair<double, double> sph_apl = calculateSphericityAplanarity(Vp4s_Boosted);
    sphericity = sph_apl.first;
    aplanarity = sph_apl.second;
    

    m_info.evtNo = evtNo;
    m_info.runNo = runNo;
    m_info.expNo = expNo;
    m_info.Ecms = Ecms;
    m_info.nGood = nGood;
    m_info.nPip = nPip;
    m_info.nPim = nPim;
    m_info.HeavyJetMass = heavyJetMass;
    m_info.HeavyJetEnergy = heavyJetEnergy;
    m_info.sphericity = sphericity;
    m_info.aplanarity = aplanarity;
    m_info.nCluster= nCluster;
    m_info.nPhoton = nPhoton;
    m_info.Evis_cms = Evis_cms;
    m_info.BalancePz_cms = BalancePz_cms;
    m_info.Energy_cms = Energy_cms;
    m_info.ECLEnergy = ECLEnergy;
    m_info.ECLEnergyWO = ECLEnergyWO;

    double foxwolfram[5] = { fw.R(0), fw.R(1), fw.R(2), fw.R(3), fw.R(4) };
    memcpy(m_info.foxWolfram, foxwolfram, sizeof(foxwolfram));
    double thrust_vec[3] = { t_cms.mag(), t_cms.z()/t_cms.mag(), t_cms.phi() };
    memcpy(m_info.thrust, thrust_vec, sizeof(thrust_vec));
    double cm_vec[4] = { kinematics.cm.px(), kinematics.cm.py(), kinematics.cm.pz(), kinematics.cm.e() };
    memcpy(m_info.cms_vecP, cm_vec, sizeof(cm_vec));

    tree->Fill();

    return;
}


void SpinAlignment::disp_stat(const char*){
    return;
} 


void SpinAlignment::term(){
    output_file->Write();
    output_file->Close();
    return;
}



std::pair<double, double> SpinAlignment::calculateHeavyJetMassEnergy(const Vp4& particles, const Hep3Vector& thrustAxis) {
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


std::pair<double, double> SpinAlignment::calculateSphericityAplanarity(const std::vector<HepLorentzVector>& particles) {
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

 void SpinAlignment::getDrDz(Mdst_charged_Manager::iterator chr_it, int masshyp, double& dr, double& dz, double& refitPx, double& refitPy, double& refitPz){
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


void SpinAlignment::other(int* , BelleEvent*, int* ){
    return;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif