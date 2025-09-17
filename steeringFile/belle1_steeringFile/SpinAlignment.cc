#include "belle.h"
#include <cmath>

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
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<Particle> VParticle;
typedef std::vector<Hep3Vector> Vp3;


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
    ler_e = 0;
    her_e = 0;
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

    ler_e=BeamEnergy::E_LER();
    her_e=BeamEnergy::E_HER();

    if(ler_e <3.0 || her_e <7.0 || ler_e > 5.0 || her_e > 9.0){
        cout<<"someting is wrong ler_e is "<< ler_e << " her_e is "<<her_e<<endl;
        return;
    }

    double x_angle = BeamEnergy::Cross_angle();
    kinematics::cm=HepLorentzVector(-her_e*sin(x_angle), 0., -her_e*cos(x_angle) + ler_e, her_e + ler_e);
    kinematics::CMBoost=kinematics::cm.boostVector();
    kinematics::firstElectronCM=HepLorentzVector(her_e*sin(x_angle), 0., her_e*cos(x_angle), her_e);
    kinematics::secondElectronCM=HepLorentzVector(0., 0. , -ler_e, ler_e);

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
    ADDBRANCH(Q, F);
    ADDBRANCH(visEnergy, F);
    ADDBRANCH(Evis, F);
    ADDBRANCH(nGood, I);
    ADDBRANCH(nPip, I);
    ADDBRANCH(nPim, I);
    ADDBRANCH(nCluster, I);
    //ADDBRANCH(ler_e, F);
    //ADDBRANCH(her_e, F);
    ADDBRANCH(thrust, F);
    ADDBRANCH(thrust_theta, F);
    ADDBRANCH(thrust_phi, F);
    ADDBRANCH(z, F);
    ADDBRANCH(pt, F);

    ADDBARRAY(cms_vecP, 4, F);
    
    tree->Branch("nPhoton", &nPhoton, "nPhoton/I");
    tree->Branch("photon_p", photon_p, "photon_p[nPhoton]/F");
    tree->Branch("photon_theta", photon_phi, "photon_theta[nPhoton]/F");
    tree->Branch("photon_phi", photon_phi, "photon_phi[nPhoton]/F");
    tree->Branch("nPip", &nPip, "nPip/I");
    tree->Branch("nPim", &nPim, "nPim/I");
    tree->Branch("pip_vecP", pip_vecP, "pip_vecP[nPip][4]/F");
    tree->Branch("pim_vecP", pim_vecP, "pim_vecP[nPim][4]/F");
    tree->Branch("pip_theta", pip_theta, "pip_theta[nPip]/F");
    tree->Branch("pim_theta", pim_theta, "pim_theta[nPim]/F");
    tree->Branch("pip_phi", pip_phi, "pip_phi[nPip]/F");
    tree->Branch("pim_phi", pim_phi, "pim_phi[nPim]/F");

    return;
}


void SpinAlignment::event(BelleEvent* evptr, int* status){
    (void)evptr; (void)status;

    int  expNo=Belle_event_Manager::get_manager().begin()->ExpNo();
    int  evtNo=Belle_event_Manager::get_manager().begin()->EvtNo();
    int  runNo=Belle_event_Manager::get_manager().begin()->RunNo();

    if(!(countEvt%10000)) 
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
    allParticles.clear();
    allParticles_Boosted.clear();

    double Energy_cms = 0;
    double Evis_cms = 0;
    double BalancePz_cms = 0;

    Mdst_charged_Manager& charged_mgr = Mdst_charged_Manager::get_manager();
    //Mdst_trk_Manager& trk_mgr = Mdst_trk_Manager::get_manager();  not used
    Mdst_gamma_Manager& gamma_mgr = Mdst_gamma_Manager::get_manager();
    Mdst_ecl_aux_Manager& eclaux_mgr = Mdst_ecl_aux_Manager::get_manager();

    char particleName[10];
    strcpy(particleName, "unknown");

    for(Mdst_charged_Manager::iterator ch_it=charged_mgr.begin(); ch_it!=charged_mgr.end(); ch_it++){
        const Mdst_charged &chg = *ch_it;

        int charge = chg.charge();
        if(charge > 0) {
            strcpy(particleName, "pi+");
            nPip ++;
        }
        else {
            strcpy(particleName, "pi-");
            nPim ++;
        }
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
        vec_p4_boosted.boost(kinematics::CMBoost);

        allParticles.push_back(vec_p3);
        allParticles_Boosted.push_back(vec_p4_boosted.vect());

        Energy_cms += vec_p4_boosted.e();
        Evis_cms += vec_p4_boosted.e();
        BalancePz_cms += vec_p4_boosted.pz();
    }

        int nGood = allParticles_Boosted.size(); 
        if (nGood < 3) return;

    for(Mdst_gamma_Manager::iterator gam = gamma_mgr.begin(); gam != gamma_mgr.end(); gam++){
        const Mdst_gamma &gamma = *gam;

        Hep3Vector vec_p3(gamma.p(0), gamma.p(1), gamma.p(2));

        Mdst_ecl_aux &aux = eclaux_mgr(Panther_ID(gamma.ecl().get_ID()));

        double e9oe25 = aux.e9oe25();
        double gammaE = vec_p3.mag();

        if(gammaE < cuts::minPhotonE) 
            continue;

        double gam_theta = vec_p3.theta();
        bool thetaInCDCAcceptance = false;
        if(gam_theta > 17 * PI / 180 && gam_theta < 150 * PI / 180) 
            thetaInCDCAcceptance = true;

        HepLorentzVector vec_p4(vec_p3, gammaE);
        HepLorentzVector vec_p4_boosted(vec_p3, gammaE);
        vec_p4_boosted.boost(kinematics::CMBoost);

        Energy_cms += vec_p4_boosted.e();
        if (thetaInCDCAcceptance){
            Evis_cms += vec_p4_boosted.e();
            BalancePz_cms += vec_p4_boosted.pz();
            nPhoton ++;
        }
        nCluster ++;
    }

    Thrust t_cms = thrustall(allParticles_Boosted.begin(), allParticles_Boosted.end(), retSelf);
    Thrust t = thrustall(allParticles.begin(), allParticles.end(), retSelf);

    // don't know if this snytax is correct
    kinematics::thrust = t.thru;
    //kinematics::thrust_theta = t.theta;
    //kinematics::thrust_phi = t.phi;

    m_info.evtNo = evtNo;
    m_info.runNo = runNo;
    m_info.expNo = expNo;

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