#include "belle.h"
#include <cmath>
//#include "TH1F.h"

#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "particle/Particle.h"
#include "particle/utility.h"
#include "particle/combination.h"
//#include "LambdaAna/LambdaAna.h"
#include "./LambdaAna.h"
#include "./AnaConsts.h"
#include "./FindLambdaC.h"
#include "./FindD0Dp.h"
#include <eid/eid.h>
#include "TTree.h"

#include <panther/panther.h>
//#include <mdst/Evtcls_hadron_info.h>
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
//#include "TMath.h"

#include MDST_H
#include BELLETDF_H
#include HEPEVT_H

using namespace std;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<float> Vfloat;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<Particle> VParticle;

	Hep3Vector& retSelf(Hep3Vector& vec)
	 {
   	return vec;
 	 };

  //Vint weightedLambda; 
//bool eventstore=false;
  //HepLorentzVector m_Lambda_mix;
  //HepLorentzVector m_proton_mix;
  //HepLorentzVector m_proton_mix_tmp;
  //HepLorentzVector m_pion_mix;
  //int hemiflag_mix_tmp;
  //int hemiflag_mix;
  //int m_prcharge_mix;
  //int m_prcharge_mix_tmp;
  //Hep3Vector thrust_mix;
  //Hep3Vector thrust_mix_tmp;
  HepLorentzVector m_Lambda_mix;
  Vp4 m_proton_mix;
  Vp4 m_proton_mix_lab;
  HepLorentzVector m_proton_mix_tmp;
  HepLorentzVector m_proton_mix_lab_tmp;
  HepLorentzVector m_pion_mix;
  HepLorentzVector m_pion_mix_lab;
  int  hemiflag_mix_tmp;
  Vint hemiflag_mix;
  Vint m_prcharge_mix;
  int m_prcharge_mix_tmp;
  Vp3 thrust_mix;
  Hep3Vector thrust_mix_tmp;
  long int Msize=20000;

LambdaAna::LambdaAna() {
  puts( "***** [LambdaAna] Module loaded successfully." );
  //strcpy(rFileName,"notInitialized.root");
  m_debug=false;
  m_mc=0;// data or MC sample
  m_rm_mctree=1;
  m_rm_tree=0;
  m_rm_mixtree=0;
  m_PID = true;
  //m_mc=new char[20];
  m_output_filename=new char[256];
  m_output_filename_genTree=new char[256];
  weightedProton.clear();
  weightedPion.clear();
  weightedProton_zptcom.clear();
  weightedPion_zptcom.clear();
  weightedProton_zptcom_b.clear();
  weightedPion_zptcom_b.clear();
  weightedProton_zptcom_c.clear();
  weightedPion_zptcom_c.clear();
  weightedFactor_lambdaId.clear();
  weightedFactor_factor.clear();
  //comweighted_Lambda.clear();
  //comweighted_Hadron.clear();
  //comweightedProton.clear();
  //comweightedPion.clear();
  //comweighted_Hadron_href.clear();
  //comweightedProton_href.clear();
  //comweightedPion_href.clear();
  comweighted_Hadron_zptcom_href.clear();
  comweightedProton_zptcom_href.clear();
  comweightedPion_zptcom_href.clear();
  comweighted_Hadron_zptcom_href_b.clear();
  comweightedProton_zptcom_href_b.clear();
  comweightedPion_zptcom_href_b.clear();
  r1.SetSeed(100);
  //r2.SetSeed(200);
  //r3.SetSeed(300);
  //r4.SetSeed(400);
  r5.SetSeed(500);
  r6.SetSeed(600);
  r7.SetSeed(700);
  r8.SetSeed(800);
  r9.SetSeed(900);
  // my_parameter = 0.0;
  // strcpy( my_filename, "test" );
  return;
}

void LambdaAna::init(int *) {
  printf( "\n[LambdaAna] parameters\n" );
  eler=0;
  eher=0;
  //m_rm_mctree=true;
  m_thrustMagCut=0.8;
  m_zhcut=0.2;
  m_zhcut_up=0.9;
  m_inputPolarization=-0.05;
  m_LamDecayPar = 0.642;
  m_antiLamDecayPar = -0.71;
  m_PI = 3.1415926;
   //m_file=new TFile(rFileName,"recreate");
  //printf( "my_parameter = %f\n", my_parameter );
  printf( "IsMCSample = %I\n", m_mc);
  printf( "output_filename = %s\n", m_output_filename);
  // Ptype dummy("VPHO");
  //m_proton_mix.setPx(0);
  //m_proton_mix.setPy(0);
  //m_proton_mix.setPz(0);
  //m_proton_mix.setE(0);
  m_pion_mix.setPx(0);
  m_pion_mix.setPy(0);
  m_pion_mix.setPz(0);
  m_pion_mix.setE(0);
  m_pion_mix_lab.setPx(0);
  m_pion_mix_lab.setPy(0);
  m_pion_mix_lab.setPz(0);
  m_pion_mix_lab.setE(0);
  m_Lambda_mix.setPx(0);
  m_Lambda_mix.setPy(0);
  m_Lambda_mix.setPz(0);
  m_Lambda_mix.setE(0);
  return;
}

void LambdaAna::begin_run(BelleEvent* evptr, int *status) {
  eid::init_data();
  (void)evptr; (void)status;

   IpProfile::begin_run();
   BeamEnergy::begin_run();

   eler=BeamEnergy::E_LER();
   eher=BeamEnergy::E_HER();
   
    if(eler <3.0 || eher <7.0 || eler > 5.0 || eher > 9.0)
      {
        //validRun=false;
        cout<<"someting is wrong eler is "<< eler << " eher is "<<eher<<endl;
        return;
      }
    //else
     // {
        //validRun=true;
      //}


    //double theta(0.022);
    double theta=BeamEnergy::Cross_angle();
    kinematics::cm=HepLorentzVector(-eher*sin(theta),0.,-eher*cos(theta)+eler,eher+eler);
    kinematics::CMBoost=kinematics::cm.boostVector();
    kinematics::firstElectronCM=HepLorentzVector(eher*sin(theta),0.,eher*cos(theta),eher);
    kinematics::secondElectronCM=HepLorentzVector(0.,0.,-eler,eler);

    //cout<<"first Electron Lab: "<< kinematics::firstElectronCM.px() << ",  "<< kinematics::firstElectronCM.py()<<", "<<kinematics::firstElectronCM.pz()<<endl;
    //cout<<"boost vector "<<kinematics::CMBoost.x() << ", "<< kinematics::CMBoost.y() <<", "<<kinematics::CMBoost.z()<<endl;

    HepLorentzVector CMBoost2=kinematics::cm;
    CMBoost2.boost(-kinematics::CMBoost); 
    kinematics::Q=CMBoost2.t();
    kinematics::firstElectronCM.boost(kinematics::CMBoost);
    kinematics::secondElectronCM.boost(kinematics::CMBoost);

    //cout<<"first ElectronCM: "<< kinematics::firstElectronCM.px() << ",  "<< kinematics::firstElectronCM.py()<<", "<<kinematics::firstElectronCM.pz()<<endl;
    //cout<<"second ElectronCM: "<< kinematics::secondElectronCM.px() << ",  "<< kinematics::secondElectronCM.py()<<", "<<kinematics::secondElectronCM.pz()<<endl;

  return;
}

void LambdaAna::end_run(BelleEvent* evptr, int *status) {
  (void)evptr; (void)status;
  return;
}

void LambdaAna::hist_def() {
// extern BelleTupleManager *BASF_Histogram;
// BelleTupleManager *tm = BASF_Histogram;
 //m_tup_evt = tm->ntuple( "evtTup", "NLambda ngood runNo evtNo expNo Flag thrust thrx thry thrz thrx_mc thry_mc thrz_mc thrdiff eler eher",1);
//m_tup = tm->ntuple( "lTup", "ngood Q Lmass Ltheta Lpx Lpy Lpz Le ppx ppy ppz pe pipx pipy pipz pie Lpxlab Lpylab Lpzlab rawmass hemiflag kind Isgood chisq Lmag runNo evtNo expNo thrust thrx thry thrz thrx_mc thry_mc thrz_mc eler eher tangle thetap coslp costhetap costp cosTest z pt prmag pimag prtheta pitheta nkpoh nkmoh npipoh npimoh",2);
  //m_tup_mc_evt = tm->ntuple( "mcEvtTup", "NLam thr_mc",3);
 // m_tup_mc = tm->ntuple( "mcTup", "thr_mc thep_mc z_mc Lmother Lm_mc Lm_ch costheta cosl_mc lp_new lp_check tp_test tp_new pcharge hflag Lmag_mc pimag_mc prmag_mc pithe_mc lthe_mc prthe_mc  pithe_cm prthe_cm lthe_cm tan_mc drop",4);
  //hmass = tm->histogram("mass",200,0.9,1.3,100);
  //hthr = tm->histogram("thrust",100,0.49,1.01,101);

  m_file=new TFile(m_output_filename,"recreate");
  if(m_rm_mctree!=1)m_mctree= new TTree("mcLambda","My Transversity MC Tree");
  if(m_rm_tree!=1)m_tree= new TTree("Lambda","My Transversity Tree");
  //m_tree= new TTree("Lambda","My Transversity Data Tree");
  if(m_rm_mixtree!=1)m_mixtree= new TTree("mix","My Transversity Mix Tree");
  //event leval tree
  m_evttree= new TTree("evttree","eveTree");
  m_evttree->Branch("nLambda_b",&m_nLambda_b,"m_nLambda_b/I");
  m_evttree->Branch("nLambda",&m_nLambda,"m_nLambda/I");
  m_evttree->Branch("nLc",&m_nLc,"m_nLc/I");
  m_evttree->Branch("nLc_mode0",&m_nLc_mode0,"m_nLc_mode0/I");
  m_evttree->Branch("nSigma0",&m_nSigma0,"m_nSigma0/I");

  #define ADDBRANCH__(name,var,type)   m_tree->Branch(#name,&m_info.var,#type)
  #define ADDBARRAY__(name,var,n,type) ADDBRANCH__(name,var,name[n]/type)
  #define ADDBRANCH(x,type)            ADDBRANCH__(x,x,x/type)
  #define ADDBARRAY(x,n,type)          ADDBARRAY__(x,x,n,type)

  #define ADDMCBRANCH__(name,var,type)   m_mctree->Branch(#name,&m_mc_info.var,#type)
  #define ADDMCBARRAY__(name,var,n,type) ADDMCBRANCH__(name,var,name[n]/type)
  #define ADDMCBRANCH(x,type)            ADDMCBRANCH__(x,x,x/type)
  #define ADDMCBARRAY(x,n,type)          ADDMCBARRAY__(x,x,n,type)

  //and create all branches
  if(m_rm_tree!=1) {
  ADDBRANCH(evtNo, I);
  ADDBRANCH(runNo, I);
  ADDBRANCH(expNo, I);
  ADDBRANCH(Q, F);
  ADDBRANCH(visEnergy, F);
  ADDBRANCH(Evis, F);
  ADDBRANCH(ngood, I);
  ADDBRANCH(eler, F);
  ADDBRANCH(eher, F);
  ADDBRANCH(thrust, F);
  ADDBRANCH(tangle, F);
  ADDBRANCH(thrdiff, F);
  ADDBRANCH(throff, F);
  ADDBRANCH(thrshift, F);
  ADDBRANCH(numquark, I);
  ADDBRANCH(refdirdiff, F);
  //ADDBRANCH(pt, F);
  ADDBRANCH(Lmass, F);
  ADDBRANCH(rawmass, F);
  ADDBRANCH(masscheck, F);
  ADDBRANCH(z, F);
  ADDBRANCH(pt, F);
  ADDBRANCH(z_mc, F);
  ADDBRANCH(pt_mc, F);
  ADDBRANCH(pt_mcThr, F);
  ADDBRANCH(hemiflag, I);
  ADDBRANCH(kind, I);
  ADDBRANCH(Isgood, I);
  ADDBRANCH(chisq, F);
  ADDBRANCH(pimag, F);
  ADDBRANCH(prmag, F);
  ADDBRANCH(lamag, F);
  ADDBRANCH(pitheta, F);
  ADDBRANCH(prtheta, F);
  ADDBRANCH(latheta, F);
  ADDBRANCH(pitheta_cms, F);
  ADDBRANCH(prtheta_cms, F);
  ADDBRANCH(latheta_cms, F);
  ADDBRANCH(lacostheta_cms, F);
  ADDBRANCH(y, F);
  ADDBRANCH(thetap, F);
  ADDBRANCH(cosl, F);
  ADDBRANCH(costheta, F);
  ADDBRANCH(sintheta, F);
  ADDBRANCH(costheta_mcThrust, F);
  ADDBRANCH(weight, F);
  ADDBRANCH(drop, I);
  //ADDBRANCH(drop_pt, I);
  ADDBRANCH(drop_zpt, I);
  ADDBRANCH(drop_zpt_b, I);//enlarge
  ADDBRANCH(drop_zpt_c, I);//enlarge
  ADDBRANCH(matchedTrkIndex, I);
  ADDBRANCH(costheta_test, F);
  ADDBRANCH(costheta_fake, F);
  ADDBRANCH(costheta_fake_b, F);
  ADDBRANCH(costheta_fake_b_mcT, F);
  //ADDBRANCH(costheta_fake_c, F);
  ADDBRANCH(costheta_d, F);
  //ADDBRANCH(costheta_d_test, F);
  ADDBRANCH(costheta_fake_d, F);
  //ADDBRANCH(thetap_fake, F);

  ADDBRANCH(nkpoh, I);
  ADDBRANCH(nkmoh, I);
  ADDBRANCH(npipoh, I);
  ADDBRANCH(npimoh, I);

  ADDBRANCH(nMatch, I);
  ADDBRANCH(mother, I);
  ADDBRANCH(isthep, I);
  ADDBRANCH(mothermother, I);
  ADDBRANCH(cosl_true, F);
  ADDBRANCH(cost_d_true, F);
  ADDBRANCH(cost_true, F);
  ADDBRANCH(refdif_true, F);
  ADDBRANCH(cost_truemom_recThr, F);
  ADDBRANCH(mom_diff, F);
  ADDBRANCH(angle_diff, F);

  ADDBRANCH(nMatch_b, I);
  ADDBRANCH(mother_b, I);
  ADDBRANCH(isthep_b, I);
  //ADDBRANCH(mothermother_b, I);

  ADDBRANCH(Lp4_lab_phi, F);
  ADDBRANCH(Lp4_lab_theta, F);
  ADDBRANCH(Prop4_lab_phi, F);
  ADDBRANCH(Prop4_lab_theta, F);
  ADDBRANCH(Pionp4_lab_phi, F);
  ADDBRANCH(Pionp4_lab_theta, F);
  ADDBRANCH(Lp4_raw_phi, F);
  ADDBRANCH(Lp4_raw_theta, F);
  ADDBRANCH(Prop4_raw_phi, F);
  ADDBRANCH(Prop4_raw_theta, F);
  ADDBRANCH(Pionp4_raw_phi, F);
  ADDBRANCH(Pionp4_raw_theta, F);

  ADDBRANCH(Lp4_match_phi, F);//in lab. sys.
  ADDBRANCH(Lp4_match_theta, F);
  ADDBRANCH(Prop4_match_phi, F);
  ADDBRANCH(Prop4_match_theta, F);
  ADDBRANCH(Pionp4_match_phi, F);
  ADDBRANCH(Pionp4_match_theta, F);

  ADDBRANCH(thrust_phi, F);
  ADDBRANCH(thrust_theta, F);
  ADDBRANCH(thrust_mc_phi, F);
  ADDBRANCH(thrust_mc_theta, F);
  ADDBRANCH(qq_axis_phi, F);
  ADDBRANCH(qq_axis_theta, F);

  ADDBRANCH(thrust_phi_lab, F);
  ADDBRANCH(thrust_theta_lab, F);
  ADDBRANCH(thrust_mc_phi_lab, F);
  ADDBRANCH(thrust_mc_theta_lab, F);

  ADDBARRAY(thrust_dir, 3, F);
  ADDBARRAY(thrust_dir_mc, 3, F);
  ADDBARRAY(qq_axis_mc, 3, F);
  ADDBARRAY(Lp4, 4, F);
  ADDBARRAY(Prop4, 4, F);
  ADDBARRAY(Piop4, 4, F);
  ADDBARRAY(Lp4_lab, 4, F);
  ADDBARRAY(Prop4_lab, 4, F);
  ADDBARRAY(Piop4_lab, 4, F);

  m_tree->Branch("m_index",&m_index,"m_index/I");
  m_tree->Branch("cost_qq",&m_cost_qq,"m_cost_qq/F");
  m_tree->Branch("prob_proton",&m_prob_proton,"m_prob_proton/F");
  m_tree->Branch("prob_pk",&m_prob_pk,"m_prob_pk/F");
  m_tree->Branch("fl",&m_fl,"m_fl/F");

  m_tree->Branch("nMatchedQuark",&m_nMatchedQuark,"m_nMatchedQuark/I");
  m_tree->Branch("matchedFlavor",&m_matchedFlavor,"m_matchedFlavor/I");

  m_tree->Branch("nphoton",&nphoton,"nphoton/I");
  m_tree->Branch("Ephoton",Ephoton,"Ephoton[nphoton]/F");
  m_tree->Branch("theta_photon",theta_photon,"theta_photon[nphoton]/F");
  m_tree->Branch("phi_photon",phi_photon,"phi_photon[nphoton]/F");
  //m_tree->Branch("helicity",helicity,"helicity[nphoton]/F");
  m_tree->Branch("inv_lamgam",inv_lamgam,"inv_lamgam[nphoton]/F");

  m_tree->Branch("m_nlambdac",&m_nlambdac,"m_nlambdac/I");
  m_tree->Branch("m_lc_mass",&m_lc_mass,"m_lc_mass[m_nlambdac]/F");
  m_tree->Branch("m_lc_charge",&m_lc_charge,"m_lc_charge[m_nlambdac]/I");
  m_tree->Branch("m_lc_mode",&m_lc_mode,"m_lc_mode[m_nlambdac]/I");
  m_tree->Branch("m_nD",&m_nD,"m_nD/I");
  m_tree->Branch("m_D_mass",&m_D_mass,"m_D_mass[m_nD]/F");
  m_tree->Branch("m_D_charge",&m_D_charge,"m_D_charge[m_nD]/I");
  m_tree->Branch("m_D_mode",&m_D_mode,"m_D_mode[m_nD]/I");

  m_tree->Branch("comkp_cost",comkp_cost,"comkp_cost[nkpoh]/F");
  m_tree->Branch("comkm_cost",comkm_cost,"comkm_cost[nkmoh]/F");
  m_tree->Branch("compip_cost",compip_cost,"compip_cost[npipoh]/F");
  m_tree->Branch("compim_cost",compim_cost,"compim_cost[npimoh]/F");

  m_tree->Branch("comkp_pt_href",comkp_pt_href,"comkp_pt_href[nkpoh]/F");
  m_tree->Branch("comkm_pt_href",comkm_pt_href,"comkm_pt_href[nkmoh]/F");
  m_tree->Branch("compip_pt_href",compip_pt_href,"compip_pt_href[npipoh]/F");
  m_tree->Branch("compim_pt_href",compim_pt_href,"compim_pt_href[npimoh]/F");

  m_tree->Branch("comkp_cost_href",comkp_cost_href,"comkp_cost_href[nkpoh]/F");
  m_tree->Branch("comkm_cost_href",comkm_cost_href,"comkm_cost_href[nkmoh]/F");
  m_tree->Branch("compip_cost_href",compip_cost_href,"compip_cost_href[npipoh]/F");
  m_tree->Branch("compim_cost_href",compim_cost_href,"compim_cost_href[npimoh]/F");

  m_tree->Branch("comkp_theta_GJ",comkp_theta_GJ,"comkp_theta_GJ[nkpoh]/F");
  m_tree->Branch("comkm_theta_GJ",comkm_theta_GJ,"comkm_theta_GJ[nkmoh]/F");
  m_tree->Branch("compip_theta_GJ",compip_theta_GJ,"compip_theta_GJ[npipoh]/F");
  m_tree->Branch("compim_theta_GJ",compim_theta_GJ,"compim_theta_GJ[npimoh]/F");

  m_tree->Branch("comkp_cost_href_true",comkp_cost_href_true,"comkp_cost_href_true[nkpoh]/F");
  m_tree->Branch("comkm_cost_href_true",comkm_cost_href_true,"comkm_cost_href_true[nkmoh]/F");
  m_tree->Branch("compip_cost_href_true",compip_cost_href_true,"compip_cost_href_true[npipoh]/F");
  m_tree->Branch("compim_cost_href_true",compim_cost_href_true,"compim_cost_href_true[npimoh]/F");

  m_tree->Branch("comkp_pt_href_true",comkp_pt_href_true,"comkp_pt_href_true[nkpoh]/F");
  m_tree->Branch("comkm_pt_href_true",comkm_pt_href_true,"comkm_pt_href_true[nkmoh]/F");
  m_tree->Branch("compip_pt_href_true",compip_pt_href_true,"compip_pt_href_true[npipoh]/F");
  m_tree->Branch("compim_pt_href_true",compim_pt_href_true,"compim_pt_href_true[npimoh]/F");

  m_tree->Branch("comkp_theta_GJ_true",comkp_theta_GJ_true,"comkp_theta_GJ_true[nkpoh]/F");
  m_tree->Branch("comkm_theta_GJ_true",comkm_theta_GJ_true,"comkm_theta_GJ_true[nkmoh]/F");
  m_tree->Branch("compip_theta_GJ_true",compip_theta_GJ_true,"compip_theta_GJ_true[npipoh]/F");
  m_tree->Branch("compim_theta_GJ_true",compim_theta_GJ_true,"compim_theta_GJ_true[npimoh]/F");

  m_tree->Branch("kp_nmatch",kp_nmatch,"kp_nmatch[nkpoh]/I");
  //m_tree->Branch("kp_drop",kp_drop,"kp_drop[nkpoh]/I");
  m_tree->Branch("km_nmatch",km_nmatch,"km_nmatch[nkmoh]/I");
  //m_tree->Branch("km_drop",km_drop,"km_drop[nkmoh]/I");
  m_tree->Branch("pip_nmatch",pip_nmatch,"pip_nmatch[npipoh]/I");
  //m_tree->Branch("pip_drop",pip_drop,"pip_drop[npipoh]/I");
  m_tree->Branch("pim_nmatch",pim_nmatch,"pim_nmatch[npimoh]/I");
  //m_tree->Branch("pim_drop",pim_drop,"pim_drop[npimoh]/I");
  m_tree->Branch("kpz_true",kpz_true,"kpz_true[nkpoh]/F");
  m_tree->Branch("kmz_true",kmz_true,"kmz_true[nkmoh]/F");
  m_tree->Branch("pipz_true",pipz_true,"pipz_true[npipoh]/F");
  m_tree->Branch("pimz_true",pimz_true,"pimz_true[npimoh]/F");
  m_tree->Branch("kp_matchedType",m_kp_matchedType,"m_kp_matchedType[nkpoh]/I");
  m_tree->Branch("km_matchedType",m_km_matchedType,"m_km_matchedType[nkmoh]/I");
  m_tree->Branch("pip_matchedType",m_pip_matchedType,"m_pip_matchedType[npipoh]/I");
  m_tree->Branch("pim_matchedType",m_pim_matchedType,"m_pim_matchedType[npimoh]/I");

 // m_tree->Branch("kp_drop_href",kp_drop_href,"kp_drop_href[nkpoh]/I");
 // m_tree->Branch("km_drop_href",km_drop_href,"km_drop_href[nkmoh]/I");
 // m_tree->Branch("pip_drop_href",pip_drop_href,"pip_drop_href[npipoh]/I");
 // m_tree->Branch("pim_drop_href",pim_drop_href,"pim_drop_href[npimoh]/I");

  m_tree->Branch("kp_drop_zpt_href",kp_drop_zpt_href,"kp_drop_zpt_href[nkpoh]/I");
  m_tree->Branch("km_drop_zpt_href",km_drop_zpt_href,"km_drop_zpt_href[nkmoh]/I");
  m_tree->Branch("pip_drop_zpt_href",pip_drop_zpt_href,"pip_drop_zpt_href[npipoh]/I");
  m_tree->Branch("pim_drop_zpt_href",pim_drop_zpt_href,"pim_drop_zpt_href[npimoh]/I");

  m_tree->Branch("kp_drop_zpt_href_b",kp_drop_zpt_href_b,"kp_drop_zpt_href_b[nkpoh]/I");//enlarge
  m_tree->Branch("km_drop_zpt_href_b",km_drop_zpt_href_b,"km_drop_zpt_href_b[nkmoh]/I");
  m_tree->Branch("pip_drop_zpt_href_b",pip_drop_zpt_href_b,"pip_drop_zpt_href_b[npipoh]/I");
  m_tree->Branch("pim_drop_zpt_href_b",pim_drop_zpt_href_b,"pim_drop_zpt_href_b[npimoh]/I");

  m_tree->Branch("kp_mom_diff",kp_mom_diff,"kp_mom_diff[nkpoh]/F");
  m_tree->Branch("km_mom_diff",km_mom_diff,"km_mom_diff[nkmoh]/F");
  m_tree->Branch("pip_mom_diff",pip_mom_diff,"pip_mom_diff[npipoh]/F");
  m_tree->Branch("pim_mom_diff",pim_mom_diff,"pim_mom_diff[npimoh]/F");

  m_tree->Branch("kp_angle_diff",kp_angle_diff,"kp_angle_diff[nkpoh]/F");
  m_tree->Branch("km_angle_diff",km_angle_diff,"km_angle_diff[nkmoh]/F");
  m_tree->Branch("pip_angle_diff",pip_angle_diff,"pip_angle_diff[npipoh]/F");
  m_tree->Branch("pim_angle_diff",pim_angle_diff,"pim_angle_diff[npimoh]/F");

  m_tree->Branch("kp_oa_lam_cms",kp_oa_lam_cms,"kp_oa_lam_cms[nkpoh]/F");//open angle with the lambda in the e.e. cms frame
  m_tree->Branch("km_oa_lam_cms",km_oa_lam_cms,"km_oa_lam_cms[nkmoh]/F");
  m_tree->Branch("pip_oa_lam_cms",pip_oa_lam_cms,"pip_oa_lam_cms[npipoh]/F");
  m_tree->Branch("pim_oa_lam_cms",pim_oa_lam_cms,"pim_oa_lam_cms[npimoh]/F");

  m_tree->Branch("kp_lam_mass",kp_lam_mass,"kp_lam_mass[nkpoh]/F");//possible resonsance  lambda + hadron 
  m_tree->Branch("km_lam_mass",km_lam_mass,"km_lam_mass[nkmoh]/F");
  m_tree->Branch("pip_lam_mass",pip_lam_mass,"pip_lam_mass[npipoh]/F");
  m_tree->Branch("pim_lam_mass",pim_lam_mass,"pim_lam_mass[npimoh]/F");

  m_tree->Branch("kpz",kpz,"kpz[nkpoh]/F");
  m_tree->Branch("kmz",kmz,"kmz[nkmoh]/F");
  m_tree->Branch("pipz",pipz,"pipz[npipoh]/F");
  m_tree->Branch("pimz",pimz,"pimz[npimoh]/F");

  m_tree->Branch("kp_theta",m_kp_theta,"m_kp_theta[nkpoh]/F");
  m_tree->Branch("km_theta",m_km_theta,"m_km_theta[nkmoh]/F");
  m_tree->Branch("pip_theta",m_pip_theta,"m_pip_theta[npipoh]/F");
  m_tree->Branch("pim_theta",m_pim_theta,"m_pim_theta[npimoh]/F");

  m_tree->Branch("kp_prob",m_kp_prob,"m_kp_prob[nkpoh]/F");
  m_tree->Branch("km_prob",m_km_prob,"m_km_prob[nkmoh]/F");

  m_tree->Branch("kp_tangle",kp_oat,"kp_oat[nkpoh]/F");//open angle between the light hardron and thrust axis
  m_tree->Branch("km_tangle",km_oat,"km_oat[nkmoh]/F");
  m_tree->Branch("pip_tangle",pip_oat,"pip_oat[npipoh]/F");
  m_tree->Branch("pim_tangle",pim_oat,"pim_oat[npimoh]/F");

  m_tree->Branch("kp_p4",kp_p4,"kp_p4[nkpoh][4]/F");//in the e+e- cms frame
  m_tree->Branch("km_p4",km_p4,"km_p4[nkmoh][4]/F");
  m_tree->Branch("pip_p4",pip_p4,"pip_p4[npipoh][4]/F");
  m_tree->Branch("pim_p4",pim_p4,"pim_p4[npimoh][4]/F");
 }


if(m_rm_mixtree!=1){
  m_mixtree->Branch("Lmass",&m_Lmass_mix,"m_Lmass_mix/F");
  m_mixtree->Branch("z",&m_z_mix,"m_z_mix/F");
  m_mixtree->Branch("pt",&m_pt_mix,"m_pt_mix/F");
  m_mixtree->Branch("cost",&m_cost_mix,"m_cost_mix/F");
  m_mixtree->Branch("cost_b",&m_cost_mix_b,"m_cost_mix_b/F");
  m_mixtree->Branch("kind",&m_kind_mix,"m_kind_mix/I");
  m_mixtree->Branch("prmag",&m_prmag_mix,"m_prmag_mix/F");
  m_mixtree->Branch("pimag",&m_pimag_mix,"m_pimag_mix/F");
  //m_mixtree->Branch("angle",&m_oa_mix,"m_oa_mix/F");
  m_mixtree->Branch("lamag",&m_lamag_mix,"m_lamag_mix/F");//in the lab
  m_mixtree->Branch("lamag_cms",&m_lamag_mix_cms,"m_lamag_mix_cms/F");//in the lab
 }

  //mc info tree, before loose any efficiencies
if(m_rm_mctree!=1){
  ADDMCBRANCH(thr_mc, F);
  ADDMCBRANCH(thep_mc, F);
  ADDMCBRANCH(z_mc, F);
  ADDMCBRANCH(pt_mc, F);
  ADDMCBRANCH(Lm_mc, F);
  ADDMCBRANCH(Lm_ch, F);
  ADDMCBRANCH(costheta, F);
  ADDMCBRANCH(sintheta, F);
  //ADDMCBRANCH(costheta_fake, F);
  ADDMCBRANCH(costheta_fake_b, F);
  //ADDMCBRANCH(costheta_fake_c, F);
  ADDMCBRANCH(costheta_d, F);
  ADDMCBRANCH(costheta_fake_d, F);
  //ADDMCBRANCH(thetap_fake, F);
  ADDMCBRANCH(cosl_mc, F);
  //ADDMCBRANCH(lp_new, F);
  //ADDMCBRANCH(lp_check, F);
  //ADDMCBRANCH(tp_test, F);
  //ADDMCBRANCH(tp_new, F);
  ADDMCBRANCH(Lmag_mc, F);
  ADDMCBRANCH(pimag_mc, F);
  ADDMCBRANCH(prmag_mc, F);
  ADDMCBRANCH(pithe_mc, F);
  ADDMCBRANCH(prthe_mc, F);
  ADDMCBRANCH(lthe_mc, F);
  ADDMCBRANCH(pithe_cm, F);
  ADDMCBRANCH(prthe_cm, F);
  ADDMCBRANCH(lthe_cm, F);
  ADDMCBRANCH(tan_mc, F);

  ADDMCBRANCH(hflag, I);
  ADDMCBRANCH(pcharge, I);
  ADDMCBRANCH(Lmother, I);
  ADDMCBRANCH(Listhep, I);
  ADDMCBRANCH(Lmothermother, I);
  ADDMCBRANCH(Nleaf, I);
  ADDMCBRANCH(drop, I);//alone Lambda or anti-Lambda
  //ADDMCBRANCH(drop_pt, I);//alone Lambda or anti-Lambda
  ADDMCBRANCH(drop_zpt, I);//alone Lambda or anti-Lambda
  ADDMCBRANCH(drop_zpt_b, I);//alone Lambda or anti-Lambda
  ADDMCBRANCH(drop_zpt_c, I);//alone Lambda or anti-Lambda
  ADDMCBRANCH(nkp_mc, I);
  ADDMCBRANCH(nkm_mc, I);
  ADDMCBRANCH(npip_mc, I);
  ADDMCBRANCH(npim_mc, I);
  ADDMCBRANCH(nquark, I);
  m_mctree->Branch("weight_mc",&m_weight_mc,"m_weight_mc/F");
  m_mctree->Branch("quarkhemi",quarkhemi,"quarkhemi[nquark]/I");
  m_mctree->Branch("quarkId",quarkId,"quarkId[nquark]/I");
  m_mctree->Branch("quarkAng",m_quarkAng,"m_quarkAng[nquark]/F");

  m_mctree->Branch("kpz_mc",kpz_mc,"kpz_mc[nkp_mc]/F");
  m_mctree->Branch("kmz_mc",kmz_mc,"kmz_mc[nkm_mc]/F");
  m_mctree->Branch("pipz_mc",pipz_mc,"pipz_mc[npip_mc]/F");
  m_mctree->Branch("pimz_mc",pimz_mc,"pimz_mc[npim_mc]/F");

  //m_mctree->Branch("kp_cost_cms_mc",kp_cost_cms_mc,"kp_cost_cms_mc[nkp_mc]/F");
  //m_mctree->Branch("km_cost_cms_mc",km_cost_cms_mc,"km_cost_cms_mc[nkm_mc]/F");
  //m_mctree->Branch("pip_cost_cms_mc",pip_cost_cms_mc,"pip_cost_cms_mc[npip_mc]/F");
  //m_mctree->Branch("pim_cost_cms_mc",pim_cost_cms_mc,"pim_cost_cms_mc[npim_mc]/F");

  m_mctree->Branch("comkp_cost_mc",comkp_cost_mc,"comkp_cost_mc[nkp_mc]/F");
  m_mctree->Branch("comkm_cost_mc",comkm_cost_mc,"comkm_cost_mc[nkm_mc]/F");
  m_mctree->Branch("compip_cost_mc",compip_cost_mc,"compip_cost_mc[npip_mc]/F");
  m_mctree->Branch("compim_cost_mc",compim_cost_mc,"compim_cost_mc[npim_mc]/F");

  m_mctree->Branch("comkp_cost_href_mc",comkp_cost_href_mc,"comkp_cost_href_mc[nkp_mc]/F");
  m_mctree->Branch("comkm_cost_href_mc",comkm_cost_href_mc,"comkm_cost_href_mc[nkm_mc]/F");
  m_mctree->Branch("compip_cost_href_mc",compip_cost_href_mc,"compip_cost_href_mc[npip_mc]/F");
  m_mctree->Branch("compim_cost_href_mc",compim_cost_href_mc,"compim_cost_href_mc[npim_mc]/F");

  m_mctree->Branch("comkp_pt_href_mc",comkp_pt_href_mc,"comkp_pt_href_mc[nkp_mc]/F");
  m_mctree->Branch("comkm_pt_href_mc",comkm_pt_href_mc,"comkm_pt_href_mc[nkm_mc]/F");
  m_mctree->Branch("compip_pt_href_mc",compip_pt_href_mc,"compip_pt_href_mc[npip_mc]/F");
  m_mctree->Branch("compim_pt_href_mc",compim_pt_href_mc,"compim_pt_href_mc[npim_mc]/F");

  m_mctree->Branch("comkp_theta_GJ_mc",comkp_theta_GJ_mc,"comkp_theta_GJ_mc[nkp_mc]/F");
  m_mctree->Branch("comkm_theta_GJ_mc",comkm_theta_GJ_mc,"comkm_theta_GJ_mc[nkm_mc]/F");
  m_mctree->Branch("compip_theta_GJ_mc",compip_theta_GJ_mc,"compip_theta_GJ_mc[npip_mc]/F");
  m_mctree->Branch("compim_theta_GJ_mc",compim_theta_GJ_mc,"compim_theta_GJ_mc[npim_mc]/F");

  //m_mctree->Branch("comkp_drop",comkp_drop,"comkp_drop[nkp_mc]/I");
  //m_mctree->Branch("comkm_drop",comkm_drop,"comkm_drop[nkm_mc]/I");
  //m_mctree->Branch("compip_drop",compip_drop,"compip_drop[npip_mc]/I");
  //m_mctree->Branch("compim_drop",compim_drop,"compim_drop[npim_mc]/I");

  //m_mctree->Branch("comkp_drop_href",comkp_drop_href,"comkp_drop_href[nkp_mc]/I");
  //m_mctree->Branch("comkm_drop_href",comkm_drop_href,"comkm_drop_href[nkm_mc]/I");
  //m_mctree->Branch("compip_drop_href",compip_drop_href,"compip_drop_href[npip_mc]/I");
  //m_mctree->Branch("compim_drop_href",compim_drop_href,"compim_drop_href[npim_mc]/I");

  m_mctree->Branch("comkp_zpt_drop_href",comkp_zpt_drop_href,"comkp_zpt_drop_href[nkp_mc]/I");
  m_mctree->Branch("comkm_zpt_drop_href",comkm_zpt_drop_href,"comkm_zpt_drop_href[nkm_mc]/I");
  m_mctree->Branch("compip_zpt_drop_href",compip_zpt_drop_href,"compip_zpt_drop_href[npip_mc]/I");
  m_mctree->Branch("compim_zpt_drop_href",compim_zpt_drop_href,"compim_zpt_drop_href[npim_mc]/I");

  m_mctree->Branch("comkp_zpt_drop_href_b",comkp_zpt_drop_href_b,"comkp_zpt_drop_href_b[nkp_mc]/I");
  m_mctree->Branch("comkm_zpt_drop_href_b",comkm_zpt_drop_href_b,"comkm_zpt_drop_href_b[nkm_mc]/I");
  m_mctree->Branch("compip_zpt_drop_href_b",compip_zpt_drop_href_b,"compip_zpt_drop_href_b[npip_mc]/I");
  m_mctree->Branch("compim_zpt_drop_href_b",compim_zpt_drop_href_b,"compim_zpt_drop_href_b[npim_mc]/I");
  m_mctree->Branch("angle_qqth",&angle_qqthr,"angle_qqthr/F");
  m_mctree->Branch("costheta_mc_qq",&m_costheta_mc_qq,"m_costheta_mc_qq/F");
 }

  /*m_mctree->Branch("thr_mc",thr_mc,"thr_mc/F");
  */

  /*int* LpionCounter=new int;
  float* memLoc=new float[1200];
  treeData->Branch(LpionCounter,LpionCounter,"LpionCounter/I");
  treeData->Branch(LpionCos,memLoc,"LpionCos[LpionCounter]/F");
*/

  return;
}

void LambdaAna::event(BelleEvent* evptr, int* status) {
  (void)evptr; (void)status;

   int  expNo=Belle_event_Manager::get_manager().begin()->ExpNo();
   int  evtNo=Belle_event_Manager::get_manager().begin()->EvtNo();
   int  runNo=Belle_event_Manager::get_manager().begin()->RunNo();
   ppion.clear();
   pkaon.clear();
   //cout<<"************************************************"<<evtNo<<endl;
   //cout<<"evt No: "<<evtNo<<endl;

    if(!(countEvt%10000)) 
    cout << "evt " <<countEvt<< " expNo "<< expNo << ", runNo "<< runNo << ", evtNo "<< evtNo <<endl;
    
   //weightedLambda.clear(); comweighted_Lambda.clear(); comweighted_Hadron.clear();
   weightedProton.clear(); weightedPion.clear();
   weightedProton_zptcom.clear(); weightedPion_zptcom.clear(); weightedProton_zptcom_b.clear(); weightedPion_zptcom_b.clear(); weightedProton_zptcom_c.clear(); weightedPion_zptcom_c.clear();
    weightedFactor_lambdaId.clear(); weightedFactor_factor.clear();
   //comweightedProton.clear(); comweightedPion.clear(); comweighted_Hadron.clear();
   //comweightedProton_href.clear(); comweightedPion_href.clear(); comweighted_Hadron_href.clear();
   comweightedProton_zptcom_href.clear(); comweightedPion_zptcom_href.clear(); comweighted_Hadron_zptcom_href.clear();
   comweightedProton_zptcom_href_b.clear(); comweightedPion_zptcom_href_b.clear(); comweighted_Hadron_zptcom_href_b.clear();
   Hep3Vector truthThrustAxisCMS(0,0,0);
   Hep3Vector qqAxisCMS(0,0,0); int numquark;
   if(m_mc==1) { truthThrustAxisCMS = readMC(qqAxisCMS, numquark);/* cout<<" reading MC "<<endl; */}///before event selection, keep all events before detector efficience
   //if(comweightedProton.size() != comweighted_Hadron.size()) cout<<"somethig is wrong, size of dropped pair is inconsistent "<<endl;
   if(comweightedProton_zptcom_href.size() != comweighted_Hadron_zptcom_href.size()) cout<<"somethig is wrong, size of dropped pair is inconsistent "<<endl;
   //cout<<" qqAxisCMS : "<<qqAxisCMS.x()<<", "<< qqAxisCMS.y() <<" , "<<qqAxisCMS.z()<<endl;

   ///read MC
  /* Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();
   for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
    {
     if(gen_it->get_ID()!=bMesonId);
    }
  */

    if(!IpProfile::usable()) {cout <<" ip not usable ..." << endl;return;}
    
    //const HepPoint3D ip_position = IpProfile::position();
    //const HepSymMatrix ip_position_err = IpProfile::position_err();
    const HepPoint3D ip_position = IpProfile::e_position();
    const HepSymMatrix ip_osition_err = IpProfile::e_position_err();
    countEvt ++;

    vector<Hep3Vector> allParticlesBoosted; allParticlesBoosted.clear();
    vector<Hep3Vector> allParticlesNonBoosted; allParticlesNonBoosted.clear();
    vector<float> nonBoostedE;
    float visEnergy_lab=0;
    float EnergyCMS = 0;
    float EvisCMS = 0;
    float BalancePzCMS = 0;

    Mdst_charged_Manager& Mdst_charged_Mgr=Mdst_charged_Manager::get_manager();
    Mdst_vee2_Manager& Mdst_vee2_Mgr = Mdst_vee2_Manager::get_manager();
    Mdst_trk_Manager& mdst_trk_Mgr=Mdst_trk_Manager::get_manager();
    //if(Mdst_vee2_Mgr.size()<1) { cout<<"size of vee2 "<<Mdst_vee2_Mgr.size()<<endl; return;  }//no Lambda candidate
    FindLambda findLtmp; int nLcan=0;
    Vint trkIndex_used; trkIndex_used.clear();

    atc_pid selKPi(3,1,5,3,2);  //K/pi separation
    atc_pid selPiP(3,1,5,2,4); 
    atc_pid selKP(3,1,5,3,4);
    Vint ipion; ipion.clear();
    Vint ikaon; ikaon.clear();
    //Vint ipionp; ipionp.clear();
    //Vint ipionm; ipionm.clear();
    Vp4 vp4pion; vp4pion.clear();
    Vp4 vp4kaon; vp4kaon.clear();
    Vp4 vp4pion_lab; vp4pion_lab.clear();
    Vp4 vp4kaon_lab; vp4kaon_lab.clear();
    Vfloat vkaonprob; vkaonprob.clear();
    pkaon.clear();
    ppion.clear();
    Vint pionchr; pionchr.clear();
    Vint kaonchr; kaonchr.clear();
    Vint pionTrkID; pionTrkID.clear();
    Vint kaonTrkID; kaonTrkID.clear();
    //Vint ikaonp; ikaonp.clear();
    //Vint ikaonm; ikaonm.clear();
    char ptypName[200];
    strcpy(ptypName,"unknown");
    int cnt_pion=0;
    int ii=0;

     for(Mdst_charged_Manager::iterator chr_it=Mdst_charged_Mgr.begin();chr_it!=Mdst_charged_Mgr.end();chr_it++,ii++)
      {
        double m_mass=m_pi;
        int massHyp=2;
        bool isLepton=false;
        bool isPionKaon=false;
        bool isPion=false;
        double m_theta=0;
        double m_phi=0;
        double m_qt=0;
        double m_z=0;
        bool positivelyIdentified=false;
        int charge=(*chr_it).charge();
        const Mdst_charged &chg = *chr_it;
         //defaults...
        if(charge>0)
          {
            strcpy(ptypName,"PI+");
          }
        else
          {
            strcpy(ptypName,"PI-");
          }
        int TrkID=(*chr_it).get_ID();
        eid sel_e(*chr_it);
        double mu_id=0;

        Muid_mdst muID(*chr_it);
        if(muID.Chi_2()>0) mu_id=muID.Muon_likelihood();

        double atcKPi=selKPi.prob(*chr_it);
        double atcKP=selKP.prob(*chr_it);
        double atcPiP=selPiP.prob(*chr_it);

        float e_cut=0.8;
        float mu_cut=0.9;
        double e_id=sel_e.prob(3,-1,5);

        double dr, dz, refitPx, refitPy, refitPz;
        getDrDz(chr_it, massHyp,dr,dz, refitPx, refitPy, refitPz);
         if ( fabs(dr) > 1  || fabs(dz) > 3 )//slides from kibayashi 
            continue;
        
 Particle p(chg,string(ptypName));
  //const Gen_hepevt &h_p=p.genHepevt();// seems not used any more
 const Gen_hepevt &h_p = gen_level(get_hepevt(p.mdstCharged()));

  //if(h_p) cout<<" find h_p"<<endl; 
  if(isPion) cnt_pion++;
  Hep3Vector h3Vect((*chr_it).p(0),(*chr_it).p(1),(*chr_it).p(2));
  double pt = h3Vect.perp();        
  if (pt < 0.1) continue; // 100 MeV cut on pt

  double E=sqrt(m_mass*m_mass+h3Vect.mag2());
  HepLorentzVector LabVec(h3Vect,E);
  HepLorentzVector boostedVec(h3Vect,E);
  boostedVec.boost(kinematics::CMBoost);
  if(positivelyIdentified&&massHyp==2) 
 {ipion.push_back(ii); vp4pion.push_back(boostedVec); vp4pion_lab.push_back(LabVec); pionchr.push_back(charge); pionTrkID.push_back(TrkID); ppion.push_back(p);}
  if(positivelyIdentified&&massHyp==3) {ikaon.push_back(ii); vkaonprob.push_back(atcKP); vp4kaon.push_back(boostedVec);vp4kaon_lab.push_back(LabVec); kaonchr.push_back(charge); kaonTrkID.push_back(TrkID); pkaon.push_back(p);}
  allParticlesBoosted.push_back(boostedVec.vect());
  allParticlesNonBoosted.push_back(h3Vect);
  nonBoostedE.push_back(E);

  visEnergy_lab+=E; //in the lab frame
  EnergyCMS += boostedVec.e();
  EvisCMS += boostedVec.e();
  BalancePzCMS += boostedVec.pz();

}//end of loop  charged particles


  int ngood = allParticlesBoosted.size();
  if (ngood < 3) return;

    Mdst_gamma_Manager& gamma_mgr=Mdst_gamma_Manager::get_manager();
    Mdst_ecl_aux_Manager& eclaux_mgr = Mdst_ecl_aux_Manager::get_manager();
    int gammaCount=0; Vint vphoton; vphoton.clear();
    int NosPhoton = 0, NosCluster = 0;
    std::vector<Mdst_gamma>::const_iterator i_begin = gamma_mgr.begin();
    //int igam=0;
   for(int igam =0; igam<gamma_mgr.size();igam++)
      {
        std::vector<Mdst_gamma>::const_iterator i = i_begin + igam;
        Hep3Vector h(i->px(),i->py(),i->pz());
        const Mdst_gamma& gam=*i;
        int id=(int)gam.get_ID();
        double px=gam.px();
        double py=gam.py();
        double pz=gam.pz();
	Mdst_ecl_aux &aux =eclaux_mgr(Panther_ID(gam.ecl().get_ID()));
        //if(gam.ecl().quality()!=0) continue;
        double e9oe25 =aux.e9oe25();
        double gammaE=sqrt(px*px+py*py+pz*pz);

       if(gammaE< 0.1) continue;
        HepLorentzVector boostedVec(px,py,pz,gammaE);
        Hep3Vector photVec(px,py,pz);
   
       float photTheta= photVec.theta()*180/m_pi;
       bool thetaInCDCAcceptance = false;
       if(photTheta >17 && photTheta < 150)
          thetaInCDCAcceptance = true;

        vphoton.push_back(igam);

        boostedVec.boost(kinematics::CMBoost);
        //v_gammaE.push_back(gammaE);

        allParticlesNonBoosted.push_back(photVec);
        //photon
        //allPB_particleClass.push_back(-1);
        allParticlesBoosted.push_back(boostedVec.vect());
        nonBoostedE.push_back(gammaE);
        // if(boostedVec.vect().mag()>0.1)
        //fjParticles.push_back(PseudoJet(boostedVec.px(),boostedVec.py(),boostedVec.pz(),boostedVec.e()));
        //allPB_E.push_back(boostedVec.e());
        visEnergy_lab+=gammaE;//in the labe frame
        EnergyCMS += boostedVec.e();
        if(thetaInCDCAcceptance)
        {
        EvisCMS += boostedVec.e();
        BalancePzCMS += boostedVec.pz();
        NosPhoton++;
        }
        gammaCount++;
        NosCluster++;
}

       Thrust t=thrustall(allParticlesBoosted.begin(),allParticlesBoosted.end(),retSelf);
       Thrust labThrust=thrustall(allParticlesNonBoosted.begin(),allParticlesNonBoosted.end(),retSelf); 

    /*
    if(cos(labThrust.axis.theta()) > cuts::maxLabThrustCosTheta || cos(labThrust.axis.theta())<cuts::minLabThrustCosTheta)
      {
       // exitEvent();
        return;
      }
    */
      kinematics::thrustDirCM=t.axis;
      kinematics::thrustZReverted=true;
      kinematics::thrustDirLab=labThrust.axis;
      kinematics::thrustMag=t.thru;

    /*
     //if(t.thru<m_thrustMagCut) return;
     if(t.thru>m_thrustMagCut) return;// we do a test
    */
      int randomFlag=-1;
     
     //cout<<" rand() "<< rand()<<endl;
      if(rand() % 100 <50) 
     { kinematics::thrustZReverted=false;
      kinematics::thrustDirCM=-1*t.axis; 
      kinematics::thrustDirLab=-1*labThrust.axis;
      randomFlag=1;
     }
     

      //hthr->accumulate(kinematics::thrustMag); 

    //  if(kinematics::thrustMag<cuts::minThrust)  return;


    Hep3Vector tDiff=t.axis-thrust(allParticlesBoosted.begin(),allParticlesBoosted.end(),retSelf);
    Evtcls_hadron_info_Manager& hadronInfo_mgr = Evtcls_hadron_info_Manager::get_manager();
    float visEnergyOnFile=hadronInfo_mgr.begin()->Evis();
    kinematics::E_miss=kinematics::Q-visEnergyOnFile;


   /* for(int i=0;i<allParticlesBoosted.size();i++)
      {
        float m_z=2*allPB_E[i]/kinematics::Q;
        float hrustProj=kinematics::thrustDirCM.dot((*it)->p().vect())/(kinematics::thrustDirCM.mag()*(*it)->p().vect().mag());
        if(hrustProj>0) sphere[i]=0;
        else  sphere[i]=1;
      }
*/
   //atc_pid selKPi(3,1,5,3,2);  //K/pi separation
   //atc_pid selPiP(3,1,5,2,4); 
   //vector<int> v_pion; v_pion.clear();
  //int cnt_pion=0;
/*  for (Mdst_charged_Manager::iterator chr_it=Mdst_charged_Mgr.begin();chr_it!=Mdst_charged_Mgr.end();chr_it++)
  {
   double charge=(*chr_it).charge();
   int id = (int) chr_it->get_ID();
   double atcKPi=selKPi.prob(*chr_it);
   double atcPiP=selPiP.prob(*chr_it);
   if(!(atcKPi<0.4&&atcPiP>0.6)) continue;
   cnt_pion++;

   if(charge>0)  { strcpy(ptypName,"PI+");}
   else {strcpy(ptypName,"PI-");
   }
  }
*/

  // if(cnt_pion<1) return;

    FindLambda findL;
    int NLambda=0;

     /*if goodLambda() = 1 : best S/sqrt(S+N) 
       with PID(atc_pid(3,1,5,4,2)>0.6) cut
       if goodLambda() = 2 : best S/sqrt(S+N) 
       without PID cut         
     */
   double diffangle = kinematics::thrustDirCM.angle(kinematics::thrustDirCMS_MC);
   double offangle = kinematics::thrustDirCM.angle(qqAxisCMS);
   double shiftangle = kinematics::thrustDirCMS_MC.angle(qqAxisCMS);
   for(Mdst_vee2_Manager::iterator vee_it=Mdst_vee2_Mgr.begin();vee_it!=Mdst_vee2_Mgr.end();vee_it++)
      {
    findL.candidates((*vee_it),ip_position);
    int goodlamda = findL.goodLambda();
    //if(goodlamda!=1&&goodlamda!=2) continue; 
    //if(goodlamda!=2) continue; 
        cout << "test" << endl;

    int m_kind = findL.kind(); 
    float m_dr  = findL.dr();   
    float m_dz  = findL.dz();     
    float m_dphi = findL.dphi();   
    float m_zdist = findL.zdist(); 
    m_fl = findL.fl();     
    //double m_pmag = findL.pmag(); 
    double m_chisq= findL.chisq(); 
    double mom = (*vee_it).p(4);
    //double mass = (*vee_it).p(4).mass();
    HepLorentzVector p4( (*vee_it).px(), (*vee_it).py(),(*vee_it).pz(),(*vee_it).energy());
    HepLorentzVector p4CMS = p4;
    p4CMS.boost(kinematics::CMBoost);
    //need check the order in Mdst_vee2 
    int kind = (*vee_it).kind();//Lambda or Lambda-bar; proton pi- or anti-proton pi+
    int pcharge=0; 
    double protonE,pionE;  HepLorentzVector pp4, pip4, pp4_raw, pip4_raw ;
    HepLorentzVector pp4CMS, pip4CMS;
    Hep3Vector ch0 ((*vee_it).chgd(0).px(),(*vee_it).chgd(0).py(),(*vee_it).chgd(0).pz());
    Hep3Vector ch1 ((*vee_it).chgd(1).px(),(*vee_it).chgd(1).py(),(*vee_it).chgd(1).pz());
    Particle L(*vee_it);
    Particle &child0 = L.relation().child(0);
    Particle &child1 = L.relation().child(1);
    Particle Pproton;
    Particle Ppion;
    HepLorentzVector ch0_fit = child0.p();
    HepLorentzVector ch1_fit = child1.p();
    int LpionTrkID=-1;
    int LprotonTrkID=-1;
    //cout<<" ch0 and ch1 "<< ch0 << " , "<< ch1<<endl;
    m_prob_proton=-10;
    m_prob_pk=-10;
    atc_pid selProton(3,1,5,4,2);
    atc_pid selPK(3,1,5,4,3);
    if(kind==2)
    {
    LpionTrkID=(*vee_it).chgd(1).get_ID();
    LprotonTrkID=(*vee_it).chgd(0).get_ID();
    protonE = sqrt(pow((*vee_it).chgd(0).px(),2)+pow((*vee_it).chgd(0).py(),2)+pow((*vee_it).chgd(0).pz(),2)+xmass[4]*xmass[4]);
    pionE = sqrt(pow((*vee_it).chgd(1).px(),2)+pow((*vee_it).chgd(1).py(),2)+pow((*vee_it).chgd(1).pz(),2)+xmass[2]*xmass[2]);
    pp4_raw.setVect(ch0),
    pp4_raw.setE(protonE);
    pip4_raw.setVect(ch1),
    pip4_raw.setE(pionE);
    pp4=ch0_fit;
    pip4=ch1_fit;
    pcharge=1;
    Pproton = L.relation().child(0);
    Ppion = L.relation().child(1);
    m_prob_proton = selProton.prob(&((*vee_it).chgd(0)));
    m_prob_pk= selPK.prob(&((*vee_it).chgd(0)));
    }
    else { ///kind ==3
    LpionTrkID=(*vee_it).chgd(0).get_ID();
    LprotonTrkID=(*vee_it).chgd(1).get_ID();
    protonE = sqrt(pow((*vee_it).chgd(1).px(),2)+pow((*vee_it).chgd(1).py(),2)+pow((*vee_it).chgd(1).pz(),2)+xmass[4]*xmass[4]);
    pionE = sqrt(pow((*vee_it).chgd(0).px(),2)+pow((*vee_it).chgd(0).py(),2)+pow((*vee_it).chgd(0).pz(),2)+xmass[2]*xmass[2]);
    pp4_raw.setVect(ch1),
    pp4_raw.setE(protonE);
    pip4_raw.setVect(ch0),
    pip4_raw.setE(pionE);
    pp4=ch1_fit;
    pip4=ch0_fit;
    pcharge=-1;
    Pproton = L.relation().child(1);
    Ppion = L.relation().child(0);
    m_prob_proton = selProton.prob(&((*vee_it).chgd(1)));
    m_prob_pk= selPK.prob(&((*vee_it).chgd(1)));
     }
    if(m_PID&&m_prob_proton<0.6) continue;

    pp4CMS=pp4; pip4CMS=pip4;
    pp4CMS.boost(kinematics::CMBoost);
    pip4CMS.boost(kinematics::CMBoost);

    HepLorentzVector p4Lab = pp4 + pip4;
    //HepLorentzVector p4CMScheck = p4Lab;
    //p4CMScheck.boost(kinematics::CMBoost);
    double mass=p4.m();
    double mass_raw=(pp4_raw+pip4_raw).m();
    HepLorentzVector Lp4_raw=pp4_raw+pip4_raw;
    double mass_check=(pp4+pip4).m();
    if(m_debug) cout<<"mass "<< mass << ", raw mass "<< mass_raw<<" constrained mass "<<mass_check<<endl; 
    if(m_debug)  cout<<"Lp4 check "<< p4Lab << ", from vee " << p4 <<endl; 
    if(m_debug) {cout<<" m_mass "<< mass <<endl;}
    if(m_debug) {cout<<"thrust axis :"<<kinematics::thrustDirCM.x()<< " , "<< kinematics::thrustDirCM.y()<< " , "<<kinematics::thrustDirCM.z()<<endl;}
    
    float m_z=2*p4CMS.e()/kinematics::Q;
    float hrustProj=kinematics::thrustDirCM.dot(p4CMS.vect())/(kinematics::thrustDirCM.mag()*p4CMS.vect().mag());
    double pt= p4CMS.vect().perp(kinematics::thrustDirCM);
    double pt_mcThr= p4CMS.vect().perp(kinematics::thrustDirCMS_MC);
    //int ptsign= p4CMS.vect().(kinematics::thrustDirCM);

    double tangle = p4CMS.vect().angle(kinematics::thrustDirCM);
    int hemiflag=-10;
    if(hrustProj>0) hemiflag=1;
    else { hemiflag=-1;
    tangle = p4CMS.vect().angle(-kinematics::thrustDirCM);
    }

    int nkaonp_otherHemi=0, nkaonm_otherHemi=0;
    int npionp_otherHemi=0, npionm_otherHemi=0;
    Vfloat vkpz, vkmz, vpipz, vpimz; vkpz.clear(); vkmz.clear();vpipz.clear();vpimz.clear();
    Vfloat vkp_oat, vkm_oat, vpip_oat, vpim_oat; //open angle from reversed thrust axis
    Vfloat vkp_mom_diff, vkm_mom_diff, vpip_mom_diff, vpim_mom_diff;
    Vfloat vkp_angle_diff, vkm_angle_diff, vpip_angle_diff, vpim_angle_diff;
    vkp_mom_diff.clear(); vkm_mom_diff.clear(); vpip_mom_diff.clear(); vpim_mom_diff.clear();
    vkp_angle_diff.clear(); vkm_angle_diff.clear(); vpip_angle_diff.clear(); vpim_angle_diff.clear();
    Vint vkp_nmatch, vkm_nmatch, vpip_nmatch, vpim_nmatch;
    Vint vkp_trackIndex, vkm_trackIndex, vpip_trackIndex, vpim_trackIndex;
    vkp_nmatch.clear(); vkm_nmatch.clear(); vpip_nmatch.clear(); vpim_nmatch.clear();
    vkp_trackIndex.clear(); vkm_trackIndex.clear(); vpip_trackIndex.clear(); vpim_trackIndex.clear();
    Vp4 vkp_oh_p4; vkp_oh_p4.clear();//kp in the opposite hemi in the cms
    Vp4 vkm_oh_p4; vkm_oh_p4.clear();
    Vp4 vpip_oh_p4; vpip_oh_p4.clear();
    Vp4 vpim_oh_p4; vpim_oh_p4.clear();
    Vp4 vkp_oh_p4_matched; vkp_oh_p4_matched.clear();//kp in the opposite hemi in the cms
    Vp4 vkm_oh_p4_matched; vkm_oh_p4_matched.clear();
    Vp4 vpip_oh_p4_matched; vpip_oh_p4_matched.clear();
    Vp4 vpim_oh_p4_matched; vpim_oh_p4_matched.clear();
    Vint kp_matchedType; kp_matchedType.clear();
    Vint km_matchedType; km_matchedType.clear();
    Vint pip_matchedType; pip_matchedType.clear();
    Vint pim_matchedType; pim_matchedType.clear();
    Vfloat vatcKp ; vatcKp.clear();
    Vfloat vatcKm ; vatcKm.clear();
    Vint Vpionindex; Vpionindex.clear();
    Vint Vkaonindex; Vkaonindex.clear();
    GetNhadron(hemiflag, vp4kaon, vp4pion,vp4kaon_lab,vp4pion_lab, kaonchr, pionchr, LprotonTrkID, LpionTrkID, pionTrkID,kaonTrkID, vkpz, vkmz, vpipz, vpimz, vkp_oat, vkm_oat, vpip_oat, vpim_oat, vkp_nmatch, vkm_nmatch, vpip_nmatch, vpim_nmatch, vkp_trackIndex, vkm_trackIndex, vpip_trackIndex, vpim_trackIndex, vkp_mom_diff,vkm_mom_diff,vpip_mom_diff,vpim_mom_diff,vkp_angle_diff, vkm_angle_diff, vpip_angle_diff, vpim_angle_diff,vkp_oh_p4,vkm_oh_p4,vpip_oh_p4,vpim_oh_p4, vkp_oh_p4_matched,vkm_oh_p4_matched,vpip_oh_p4_matched,vpim_oh_p4_matched, kp_matchedType, km_matchedType, pip_matchedType, pim_matchedType, Vkaonindex, Vpionindex);
    if(vkp_nmatch.size() != vkp_trackIndex.size() || vkp_trackIndex.size()!=vkp_mom_diff.size()  ) cout<<"something is wrong, size of nMatch and size of trackIndex is inconsistent "<<endl;
    if(vkp_oh_p4.size()!=vkp_nmatch.size() || vkp_oh_p4.size()!=vkpz.size() ) cout<<"something is wrong. number of kp in the opposite hemi is not consistent! "<<vkp_oh_p4.size() <<" , "<<vkpz.size()<<endl;
    nkaonp_otherHemi = vkpz.size();
    nkaonm_otherHemi = vkmz.size();
    npionp_otherHemi = vpipz.size();
    npionm_otherHemi = vpimz.size();
    for(int kk = 0; kk< Vkaonindex.size();kk++)
    {
    if(kaonchr[Vkaonindex[kk]]>0) vatcKp.push_back(vkaonprob[Vkaonindex[kk]]);
    if(kaonchr[Vkaonindex[kk]]<0) vatcKm.push_back(vkaonprob[Vkaonindex[kk]]);
    }

    if(vatcKp.size()!=nkaonp_otherHemi) cout<<"ERROR ! "<<endl;
    if(vatcKm.size()!=nkaonm_otherHemi) cout<<"ERROR ! "<<endl;

    ///calculate the costheta_p 
    double lptheta =-1;
    double thetap_fake=-1;
    double costheta_fake=-1;
    double costheta_test=-1;
    double costheta_fake_b=-1;
    double costheta_fake_b_mcT=-1;
    double costheta_fake_c=-1;
    double costheta_d=-1;
    double costheta_fake_d=-1;
    Hep3Vector refdir_mcThrust(0,0,0);
    double thetap_mcThrust = CalTheta(kinematics::thrustDirCMS_MC, kinematics::thrustDirLab_MC, p4, pp4, pip4, lptheta,costheta_test, thetap_fake,costheta_fake, costheta_fake_b_mcT, costheta_fake_c, costheta_d, costheta_fake_d, refdir_mcThrust);
    //Hep3Vector refdir_tmp(0,0,0);
    //double theta_truemom_recThr= CalTheta(kinematics::thrustDirCM, kinematics::thrustDirLab, Lp4Matched_lab, pp4Matched_lab, pip4Matched_lab, lptheta, thetap_fake,costheta_fake,costheta_fake_b, costheta_fake_c, costheta_d,costheta_fake_d,refdir_tmp);
    Hep3Vector refdir(0,0,0);
    ///overwrite thetap_fake,costheta_fake
    double thetap = CalTheta(kinematics::thrustDirCM, kinematics::thrustDirLab, p4, pp4, pip4, lptheta, costheta_test, thetap_fake,costheta_fake,costheta_fake_b, costheta_fake_c, costheta_d,costheta_fake_d,refdir);
    double costheta = cos(thetap); 
    double sintheta = sin(thetap); 
    double costheta_mcThrust = cos(thetap_mcThrust); 
    //double cost_truemom_recThr = cos(theta_truemom_recThr);
    double cosl = cos(lptheta);
    double protonmag = pp4.vect().mag();
    double pionmag = pip4.vect().mag();
    double protonTheta = pp4.vect().cosTheta();
    double pionTheta = pip4.vect().cosTheta();
    double protonTheta_cms = pp4CMS.vect().cosTheta();
    double pionTheta_cms = pip4CMS.vect().cosTheta();
    double LTheta_lab = p4.vect().cosTheta();
    double LcosTheta_cms = p4CMS.vect().cosTheta();
    double LTheta_cms = p4CMS.vect().theta();
    double my_y = (p4CMS.vect().dot(kinematics::firstElectronCM.vect()))/(fabs(p4CMS.vect().dot(kinematics::thrustDirCM))); 
   
   //for mixed events, the logic in this code is, only the last Lambda event is used to mix with that in other event
   m_proton_mix_tmp = pp4CMS;
   m_proton_mix_lab_tmp = pp4;
   m_prcharge_mix_tmp = pcharge;
   hemiflag_mix_tmp = hemiflag;
   thrust_mix_tmp = kinematics::thrustDirCM;
   m_pion_mix = pip4CMS;
   m_pion_mix_lab = pip4;

    NLambda++;
    //hmass->accumulate(mass);

  
  m_info.evtNo=evtNo;
  m_info.runNo=runNo;
  m_info.expNo=expNo;
  m_info.ngood=ngood;
  m_info.Q=kinematics::Q;
  m_info.eler=eler;
  m_info.eher=eher;
  m_info.visEnergy=visEnergy_lab;
  m_info.Evis=visEnergyOnFile;
  m_info.thrust=kinematics::thrustMag;
  m_info.thrust_dir[0]=kinematics::thrustDirCM.x();
  m_info.thrust_dir[1]=kinematics::thrustDirCM.y();
  m_info.thrust_dir[2]=kinematics::thrustDirCM.z();
  m_info.thrust_phi=kinematics::thrustDirCM.phi();
  m_info.thrust_theta=kinematics::thrustDirCM.theta();
  m_info.thrust_phi_lab=kinematics::thrustDirLab.phi();
  m_info.thrust_theta_lab=kinematics::thrustDirLab.theta();
  m_info.thrust_dir_mc[0]=truthThrustAxisCMS.x();
  m_info.thrust_dir_mc[1]=truthThrustAxisCMS.y();
  m_info.thrust_dir_mc[2]=truthThrustAxisCMS.z();
  m_info.qq_axis_mc[0]=qqAxisCMS.x();
  m_info.qq_axis_mc[1]=qqAxisCMS.y();
  m_info.qq_axis_mc[2]=qqAxisCMS.z();
  m_info.thrust_mc_phi=truthThrustAxisCMS.phi();
  m_info.thrust_mc_theta=truthThrustAxisCMS.theta();
  m_info.qq_axis_phi=qqAxisCMS.phi();
  m_info.qq_axis_theta=qqAxisCMS.theta();
  m_info.thrust_mc_phi_lab=kinematics::thrustDirLab_MC.phi();
  m_info.thrust_mc_theta_lab=kinematics::thrustDirLab_MC.theta();
  m_info.tangle=tangle;
  m_info.thrdiff=diffangle;
  m_info.numquark=numquark;
  m_info.throff=offangle;
  m_info.thrshift=shiftangle;

  m_info.Lmass=mass;
  m_info.rawmass=mass_raw;
  m_info.masscheck=mass_check;
  m_info.z=m_z;
  m_info.pt=pt;
  m_info.pt_mcThr=pt_mcThr;
  m_info.hemiflag=hemiflag;
  m_info.kind=m_kind;
  m_info.Isgood=goodlamda;
  m_info.chisq=m_chisq;
  m_info.Lp4[0]=p4CMS.px();
  m_info.Lp4[1]=p4CMS.py();
  m_info.Lp4[2]=p4CMS.pz();
  m_info.Lp4[3]=p4CMS.e();
  m_info.Prop4[0]=pp4CMS.px();
  m_info.Prop4[1]=pp4CMS.py();
  m_info.Prop4[2]=pp4CMS.pz();
  m_info.Prop4[3]=pp4CMS.e();
  m_info.Piop4[0]=pip4CMS.px();
  m_info.Piop4[1]=pip4CMS.py();
  m_info.Piop4[2]=pip4CMS.pz();
  m_info.Piop4[3]=pip4CMS.e();
  m_info.Lp4_lab[0]=p4.px();
  m_info.Lp4_lab[1]=p4.py();
  m_info.Lp4_lab[2]=p4.pz();
  m_info.Lp4_lab[3]=p4.e();
  m_info.Lp4_lab_phi=p4.vect().phi();
  m_info.Lp4_lab_theta=p4.vect().theta();
  m_info.Prop4_lab_phi=pp4.vect().phi();
  m_info.Prop4_lab_theta=pp4.vect().theta();
  m_info.Pionp4_lab_phi=pip4.vect().phi();
  m_info.Pionp4_lab_theta=pip4.vect().theta();
  m_info.Lp4_raw_phi=Lp4_raw.vect().phi();
  m_info.Lp4_raw_theta=Lp4_raw.vect().theta();
  m_info.Prop4_raw_phi=pp4_raw.vect().phi();
  m_info.Prop4_raw_theta=pp4_raw.vect().theta();
  m_info.Pionp4_raw_phi=pip4_raw.vect().phi();
  m_info.Pionp4_raw_theta=pip4_raw.vect().theta();
  m_info.Prop4_lab[0]=pp4.px();
  m_info.Prop4_lab[1]=pp4.py();
  m_info.Prop4_lab[2]=pp4.pz();
  m_info.Prop4_lab[3]=pp4.e();
  m_info.Piop4_lab[0]=pip4.px();
  m_info.Piop4_lab[1]=pip4.py();
  m_info.Piop4_lab[2]=pip4.pz();
  m_info.Piop4_lab[3]=pip4.e();
  m_info.prmag=protonmag;//lab frame
  m_info.pimag=pionmag;//lab frame
  m_info.lamag=p4.vect().mag();//lab frame
  m_info.prtheta=protonTheta;//lab frame
  m_info.pitheta=pionTheta;//lab frame
  m_info.latheta=LTheta_lab;//lab frame
  m_info.prtheta_cms=protonTheta_cms;//CMS frame
  m_info.pitheta_cms=pionTheta_cms;//CMS frame
  m_info.latheta_cms=LTheta_cms;//CMS frame
  m_info.y=my_y;//CMS frame
  m_info.lacostheta_cms=LcosTheta_cms;//CMS frame

  //m_info.prOAThr_Lsys=;//proton is in Lambda frame
  //m_info.prOAbeam_Lsys=;//proton is in Lambda frame

  m_info.refdirdiff=refdir.angle(refdir_mcThrust); 
  m_info.thetap=thetap; 
  m_info.costheta=costheta;
  m_info.sintheta=sintheta;
  m_info.costheta_mcThrust=costheta_mcThrust;
  //m_info.cost_truemom_recThr=cost_truemom_recThr;
  //m_info.thetap_fake=thetap_fake; 
  m_info.costheta_test=costheta_test;
  m_info.costheta_fake=costheta_fake;
  m_info.costheta_fake_b=costheta_fake_b;
  m_info.costheta_fake_b_mcT=costheta_fake_b_mcT;
  //m_info.costheta_fake_c=costheta_fake_c;
  m_info.costheta_d=costheta_d;
  m_info.costheta_fake_d=costheta_fake_d;
  m_info.cosl=cosl;
  m_info.nkpoh=nkaonp_otherHemi;
  m_info.nkmoh=nkaonm_otherHemi;
  m_info.npipoh=npionp_otherHemi;
  m_info.npimoh=npionm_otherHemi;
  
  for(int tt=0;tt<nkaonp_otherHemi;tt++) 
  { kpz[tt]= vkpz[tt]; kp_oat[tt]= vkp_oat[tt]; kp_nmatch[tt]=vkp_nmatch[tt]; m_kp_matchedType[tt] = kp_matchedType[tt]; kp_mom_diff[tt]=vkp_mom_diff[tt]; kp_angle_diff[tt]=vkp_angle_diff[tt]; comkp_cost[tt]=costheta; kp_p4[tt][0]=vkp_oh_p4[tt].px(); kp_p4[tt][1]=vkp_oh_p4[tt].py(); kp_p4[tt][2]=vkp_oh_p4[tt].pz(); kp_p4[tt][3]=vkp_oh_p4[tt].e(); 
kp_oa_lam_cms[tt]=vkp_oh_p4[tt].vect().angle(p4CMS.vect());
kp_lam_mass[tt] = (vkp_oh_p4[tt] + p4CMS).m();
float m_kpz_true=2*vkp_oh_p4_matched[tt].e()/kinematics::Q;
kpz_true[tt]=m_kpz_true;
m_kp_prob[tt] = vatcKp[tt];
m_kp_theta[tt]=vkp_oh_p4[tt].vect().theta();
}
  for(int tt=0;tt<nkaonm_otherHemi;tt++) 
  { kmz[tt]= vkmz[tt]; km_oat[tt]= vkm_oat[tt]; km_nmatch[tt]=vkm_nmatch[tt]; m_km_matchedType[tt] = km_matchedType[tt]; km_mom_diff[tt]=vkm_mom_diff[tt]; km_angle_diff[tt]=vkm_angle_diff[tt]; comkm_cost[tt]=costheta; km_p4[tt][0]=vkm_oh_p4[tt].px(); km_p4[tt][1]=vkm_oh_p4[tt].py(); km_p4[tt][2]=vkm_oh_p4[tt].pz(); km_p4[tt][3]=vkm_oh_p4[tt].e();
km_oa_lam_cms[tt]=vkm_oh_p4[tt].vect().angle(p4CMS.vect());
km_lam_mass[tt] = (vkm_oh_p4[tt]+p4CMS).m();
float m_kmz_true=2*vkm_oh_p4_matched[tt].e()/kinematics::Q;
kmz_true[tt]=m_kmz_true;
m_km_prob[tt] = vatcKm[tt];
m_km_theta[tt]=vkm_oh_p4[tt].vect().theta();
}
  for(int tt=0;tt<npionp_otherHemi;tt++) 
  { pipz[tt]= vpipz[tt]; pip_oat[tt]= vpip_oat[tt]; pip_nmatch[tt]=vpip_nmatch[tt]; m_pip_matchedType[tt] = pip_matchedType[tt]; pip_mom_diff[tt]=vpip_mom_diff[tt]; pip_angle_diff[tt]=vpip_angle_diff[tt]; compip_cost[tt]=costheta; pip_p4[tt][0]=vpip_oh_p4[tt].px(); pip_p4[tt][1]=vpip_oh_p4[tt].py();pip_p4[tt][2]=vpip_oh_p4[tt].pz(); pip_p4[tt][3]=vpip_oh_p4[tt].e(); pip_oa_lam_cms[tt]=vpip_oh_p4[tt].vect().angle(p4CMS.vect());
pip_lam_mass[tt] = (vpip_oh_p4[tt]+p4CMS).m();
float m_pipz_true=2*vpip_oh_p4_matched[tt].e()/kinematics::Q;
pipz_true[tt]=m_pipz_true;
m_pip_theta[tt]=vpip_oh_p4[tt].vect().theta();
}
  for(int tt=0;tt<npionm_otherHemi;tt++) 
  { pimz[tt]= vpimz[tt]; pim_oat[tt]= vpim_oat[tt]; pim_nmatch[tt]=vpim_nmatch[tt]; m_pim_matchedType[tt] = pim_matchedType[tt]; pim_mom_diff[tt]=vpim_mom_diff[tt]; pim_angle_diff[tt]=vpim_angle_diff[tt]; compim_cost[tt]=costheta; pim_p4[tt][0]=vpim_oh_p4[tt].px();pim_p4[tt][1]=vpim_oh_p4[tt].py(); pim_p4[tt][2]=vpim_oh_p4[tt].pz();pim_p4[tt][3]=vpim_oh_p4[tt].e();  pim_oa_lam_cms[tt]=vpim_oh_p4[tt].vect().angle(p4CMS.vect()); 
pim_lam_mass[tt] = (vpim_oh_p4[tt]+p4CMS).m();
float m_pimz_true=2*vpim_oh_p4_matched[tt].e()/kinematics::Q;
pimz_true[tt]=m_pimz_true;
m_pim_theta[tt]=vpim_oh_p4[tt].vect().theta();
}

  //take pion/kaon in the opposei hemi as the reference, that pion-Lambda plane or kaon-Lambda plane
   for(int tt=0;tt<nkaonp_otherHemi;tt++)
   {
   ///vkp_oh_p4 alread in cms frame
   //float  pt_href=0;
   comkp_cost_href[tt] = Cal_HadronFrame(vkp_oh_p4[tt], p4, pp4, comkp_pt_href[tt]);
   comkp_theta_GJ[tt] = Cal_HadronGJFrame(vkp_oh_p4[tt], p4);
   }
  
    for(int tt=0;tt<nkaonm_otherHemi;tt++)
   {
   comkm_cost_href[tt] = Cal_HadronFrame(vkm_oh_p4[tt], p4, pp4,comkm_pt_href[tt]);
   comkm_theta_GJ[tt] = Cal_HadronGJFrame(vkm_oh_p4[tt], p4);
   }

   for(int tt=0;tt<npionp_otherHemi;tt++)
   {
   compip_cost_href[tt] = Cal_HadronFrame(vpip_oh_p4[tt], p4, pp4, compip_pt_href[tt]);
   compip_theta_GJ[tt] = Cal_HadronGJFrame(vpip_oh_p4[tt], p4 );
   }

   for(int tt=0;tt<npionm_otherHemi;tt++)
   {
   compim_cost_href[tt] = Cal_HadronFrame(vpim_oh_p4[tt], p4, pp4,compim_pt_href[tt]);
   compim_theta_GJ[tt] = Cal_HadronGJFrame(vpim_oh_p4[tt], p4);
   }


  ///match to Lambdas in MC
  int nMatch=-1;
  int mother=0;
  int mothermother=0;
  int nMatch_b=-1;
  int mother_b=0;
  int mothermother_b=0;
  int isthep_b=0;
  int isthep=0;
  double cosl_true=-10; double cost_true=-10; double cost_d_true=-10; double cost_true_qq=-10; 
  double mom_diff=100; double angle_diff=-100;
  HepLorentzVector pp4Matched_lab(0,0,0,0); HepLorentzVector pp4Matched_lab_b(0,0,0,0);
  HepLorentzVector  pip4Matched_lab(0,0,0,0); HepLorentzVector  pip4Matched_lab_b(0,0,0,0);
  HepLorentzVector Lp4Matched_lab(0,0,0,0);  HepLorentzVector Lp4Matched_lab_b(0,0,0,0);
  Hep3Vector refdir_true(0,0,0);
  
  int Lmatched_trkIndex=-1;//panther ID
  int Prmatched_trkIndex=-1;//panther ID
  int Pimatched_trkIndex=-1;//panther ID
  float m_weight=1.0;
  int m_ifDrop= 0;
  int m_ifDrop_zptcom=0;
  int m_ifDrop_zptcom_b=0;
  int m_ifDrop_zptcom_c=0;
  int nMatchedQuark=0;
  int MatchedFlavor=0;
  int matched_trkIndex=-1;
  if(m_mc==1)nMatch_b = MatchToMCTruth(pcharge, mother_b, mothermother_b, isthep_b, pp4, pip4, pp4Matched_lab_b, pip4Matched_lab_b, Lp4Matched_lab_b, matched_trkIndex, mom_diff, angle_diff);
  if(m_mc==1)nMatch = MatchToMCTruth_new(pcharge, mother,mothermother, isthep, p4 ,Pproton, Ppion, L, cosl_true, cost_true, cost_true_qq, cost_d_true, Lmatched_trkIndex, Prmatched_trkIndex, Pimatched_trkIndex, mom_diff, angle_diff,Lp4Matched_lab, pp4Matched_lab, pip4Matched_lab, refdir_true, nMatchedQuark, MatchedFlavor);

   Hep3Vector refdir_tmp(0,0,0);
   double theta_truemom_recThr= CalTheta(kinematics::thrustDirCM, kinematics::thrustDirLab, Lp4Matched_lab, pp4Matched_lab, pip4Matched_lab, lptheta, costheta_test, thetap_fake,costheta_fake,costheta_fake_b, costheta_fake_c, costheta_d,costheta_fake_d,refdir_tmp);
   double cost_truemom_recThr = cos(theta_truemom_recThr);
  for(int jj=0;jj<weightedFactor_lambdaId.size();jj++)
  {
   if(Lmatched_trkIndex==weightedFactor_lambdaId[jj]) {m_weight = weightedFactor_factor[jj]; break;}
  }
 
  for(int jj=0;jj<weightedProton.size();jj++)
  {
  if (Prmatched_trkIndex==weightedProton[jj] || Pimatched_trkIndex==weightedPion[jj] ||Pimatched_trkIndex==weightedProton[jj] || Prmatched_trkIndex==weightedPion[jj] ) { m_ifDrop=1; break; }
  }

   for(int jj=0;jj<weightedProton_zptcom.size();jj++)
  {
  if (Prmatched_trkIndex==weightedProton_zptcom[jj] || Pimatched_trkIndex==weightedPion_zptcom[jj] ||Pimatched_trkIndex==weightedProton_zptcom[jj] || Prmatched_trkIndex==weightedPion_zptcom[jj] ) { m_ifDrop_zptcom=1; break; }
  }

 for(int jj=0;jj<weightedProton_zptcom_b.size();jj++)
  {
  if (Prmatched_trkIndex==weightedProton_zptcom_b[jj] || Pimatched_trkIndex==weightedPion_zptcom_b[jj] ||Pimatched_trkIndex==weightedProton_zptcom_b[jj] || Prmatched_trkIndex==weightedPion_zptcom_b[jj] ) { m_ifDrop_zptcom_b=1; break; }
  }

 for(int jj=0;jj<weightedProton_zptcom_c.size();jj++)
  {
  if (Prmatched_trkIndex==weightedProton_zptcom_c[jj] || Pimatched_trkIndex==weightedPion_zptcom_c[jj] ||Pimatched_trkIndex==weightedProton_zptcom_c[jj] || Prmatched_trkIndex==weightedPion_zptcom_c[jj] ) { m_ifDrop_zptcom_c=1; break; }
  }

  for(int tt=0;tt<nkaonp_otherHemi;tt++)
  {
  int hindex= vkp_trackIndex[tt];
  int m_dropPair=0;
  /*for(int dd=0;dd<comweightedProton.size(); dd++) 
  {
  if ((Prmatched_trkIndex==comweightedProton[dd]||Pimatched_trkIndex==comweightedPion[dd]||Pimatched_trkIndex==comweightedProton[dd]||Prmatched_trkIndex==comweightedPion[dd])&&hindex==comweighted_Hadron[dd]) { m_dropPair=1; break; }
  }
  kp_drop[tt] = m_dropPair;
  int m_dropPair_href=0;
 for(int dd=0;dd<comweightedProton_href.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton_href[dd]||Pimatched_trkIndex==comweightedPion_href[dd]||Pimatched_trkIndex==comweightedProton_href[dd]||Prmatched_trkIndex==comweightedPion_href[dd])&&hindex==comweighted_Hadron_href[dd]) { m_dropPair_href=1; break; }
  }
  kp_drop_href[tt] = m_dropPair_href;
  */

  int m_dropPair_href=0;
 for(int dd=0;dd<comweightedProton_zptcom_href.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton_zptcom_href[dd]||Pimatched_trkIndex==comweightedPion_zptcom_href[dd]||Pimatched_trkIndex==comweightedProton_zptcom_href[dd]||Prmatched_trkIndex==comweightedPion_zptcom_href[dd])&&hindex==comweighted_Hadron_zptcom_href[dd]) { m_dropPair_href=1; break; }
  }
  kp_drop_zpt_href[tt] = m_dropPair_href;

 int m_dropPair_href_b=0;
 for(int dd=0;dd<comweightedProton_zptcom_href_b.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton_zptcom_href_b[dd]||Pimatched_trkIndex==comweightedPion_zptcom_href_b[dd]||Pimatched_trkIndex==comweightedProton_zptcom_href_b[dd]||Prmatched_trkIndex==comweightedPion_zptcom_href_b[dd])&&hindex==comweighted_Hadron_zptcom_href_b[dd]) { m_dropPair_href_b=1; break; }
  }
  kp_drop_zpt_href_b[tt] = m_dropPair_href_b;
 
  }

 for(int tt=0;tt<nkaonm_otherHemi;tt++)
  {
  int hindex= vkm_trackIndex[tt];
  int m_dropPair=0;
  /*for(int dd=0;dd<comweightedProton.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton[dd]||Pimatched_trkIndex==comweightedPion[dd]||Pimatched_trkIndex==comweightedProton[dd]||Prmatched_trkIndex==comweightedPion[dd]) && hindex == comweighted_Hadron[dd]) { m_dropPair=1; break; }
  }
  km_drop[tt] = m_dropPair;
  int m_dropPair_href=0;
  for(int dd=0;dd<comweightedProton_href.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton_href[dd]||Pimatched_trkIndex==comweightedPion_href[dd]||Pimatched_trkIndex==comweightedProton_href[dd]||Prmatched_trkIndex==comweightedPion_href[dd]) && hindex == comweighted_Hadron_href[dd]) { m_dropPair_href=1; break; }
  }
  km_drop_href[tt] = m_dropPair_href;
  */

  int m_dropPair_href=0;
  for(int dd=0;dd<comweightedProton_zptcom_href.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton_zptcom_href[dd]||Pimatched_trkIndex==comweightedPion_zptcom_href[dd]||Pimatched_trkIndex==comweightedProton_zptcom_href[dd]||Prmatched_trkIndex==comweightedPion_zptcom_href[dd]) && hindex == comweighted_Hadron_zptcom_href[dd]) { m_dropPair_href=1; break; }
  }
  km_drop_zpt_href[tt] = m_dropPair_href;

  int m_dropPair_href_b=0;
  for(int dd=0;dd<comweightedProton_zptcom_href_b.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton_zptcom_href_b[dd]||Pimatched_trkIndex==comweightedPion_zptcom_href_b[dd]||Pimatched_trkIndex==comweightedProton_zptcom_href_b[dd]||Prmatched_trkIndex==comweightedPion_zptcom_href_b[dd]) && hindex == comweighted_Hadron_zptcom_href_b[dd]) { m_dropPair_href_b=1; break; }
  }
  km_drop_zpt_href_b[tt] = m_dropPair_href_b;
 
  }

  for(int tt=0;tt<npionp_otherHemi;tt++)
  {
  //cout<<"pip trkIndex "<<hindex <<endl;
  int hindex= vpip_trackIndex[tt];
  int m_dropPair=0;
  /*for(int dd=0;dd<comweightedProton.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton[dd]||Pimatched_trkIndex==comweightedPion[dd]||Pimatched_trkIndex==comweightedProton[dd]||Prmatched_trkIndex==comweightedPion[dd]) && hindex == comweighted_Hadron[dd]) { m_dropPair=1; break; }
  }
  pip_drop[tt] = m_dropPair;
   
  int m_dropPair_href=0;
  for(int dd=0;dd<comweightedProton_href.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton_href[dd]||Pimatched_trkIndex==comweightedPion_href[dd]||Pimatched_trkIndex==comweightedProton_href[dd]||Prmatched_trkIndex==comweightedPion_href[dd]) && hindex == comweighted_Hadron_href[dd]) { m_dropPair_href=1; break; }
  }
  pip_drop_href[tt] = m_dropPair_href;
  */

  int m_dropPair_href=0;
  for(int dd=0;dd<comweightedProton_zptcom_href.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton_zptcom_href[dd]||Pimatched_trkIndex==comweightedPion_zptcom_href[dd]||Pimatched_trkIndex==comweightedProton_zptcom_href[dd]||Prmatched_trkIndex==comweightedPion_zptcom_href[dd]) && hindex == comweighted_Hadron_zptcom_href[dd]) { m_dropPair_href=1; break; }
  }
  pip_drop_zpt_href[tt] = m_dropPair_href;

  int m_dropPair_href_b=0;
  for(int dd=0;dd<comweightedProton_zptcom_href_b.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton_zptcom_href_b[dd]||Pimatched_trkIndex==comweightedPion_zptcom_href_b[dd]||Pimatched_trkIndex==comweightedProton_zptcom_href_b[dd]||Prmatched_trkIndex==comweightedPion_zptcom_href_b[dd]) && hindex == comweighted_Hadron_zptcom_href_b[dd]) { m_dropPair_href_b=1; break; }
  }
  pip_drop_zpt_href_b[tt] = m_dropPair_href_b;

 }

 for(int tt=0;tt<npionm_otherHemi;tt++)
  {
  int hindex= vpim_trackIndex[tt];
  /*int m_dropPair=0;
  for(int dd=0;dd<comweightedProton.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton[dd]||Pimatched_trkIndex==comweightedPion[dd]||Pimatched_trkIndex==comweightedProton[dd]||Prmatched_trkIndex==comweightedPion[dd]) && hindex == comweighted_Hadron[dd]) { m_dropPair=1; break; }
  }
  pim_drop[tt] = m_dropPair;
  int m_dropPair_href=0;
  for(int dd=0;dd<comweightedProton_href.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton_href[dd]||Pimatched_trkIndex==comweightedPion_href[dd]||Pimatched_trkIndex==comweightedProton_href[dd]||Prmatched_trkIndex==comweightedPion_href[dd]) && hindex == comweighted_Hadron_href[dd]) { m_dropPair_href=1; break; }
  }
  pim_drop_href[tt] = m_dropPair_href;
  */

  int  m_dropPair_href=0;
  for(int dd=0;dd<comweightedProton_zptcom_href.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton_zptcom_href[dd]||Pimatched_trkIndex==comweightedPion_zptcom_href[dd]||Pimatched_trkIndex==comweightedProton_zptcom_href[dd]||Prmatched_trkIndex==comweightedPion_zptcom_href[dd]) && hindex == comweighted_Hadron_zptcom_href[dd]) { m_dropPair_href=1; break; }
  }
  pim_drop_zpt_href[tt] = m_dropPair_href;

 int m_dropPair_href_b=0;
  for(int dd=0;dd<comweightedProton_zptcom_href_b.size(); dd++)
  {
  if ((Prmatched_trkIndex==comweightedProton_zptcom_href_b[dd]||Pimatched_trkIndex==comweightedPion_zptcom_href_b[dd]||Pimatched_trkIndex==comweightedProton_zptcom_href_b[dd]||Prmatched_trkIndex==comweightedPion_zptcom_href_b[dd]) && hindex == comweighted_Hadron_zptcom_href_b[dd]) { m_dropPair_href_b=1; break; }
  }
  pim_drop_zpt_href_b[tt] = m_dropPair_href_b;

  }

 //take pion/kaon in the opposei hemi as the reference, that pion-Lambda plane or kaon-Lambda plane
   float tmp_pt;
   for(int tt=0;tt<nkaonp_otherHemi;tt++)
   {
   ///vkp_oh_p4 alread in cms frame
   comkp_cost_href_true[tt] = Cal_HadronFrame(vkp_oh_p4_matched[tt], Lp4Matched_lab, pp4Matched_lab, comkp_pt_href_true[tt]);
   comkp_theta_GJ_true[tt] = Cal_HadronGJFrame(vkp_oh_p4_matched[tt], Lp4Matched_lab);
   }

    for(int tt=0;tt<nkaonm_otherHemi;tt++)
   {
   comkm_cost_href_true[tt] = Cal_HadronFrame(vkm_oh_p4_matched[tt], Lp4Matched_lab, pp4Matched_lab, comkm_pt_href_true[tt]);
   comkm_theta_GJ_true[tt] = Cal_HadronGJFrame(vkm_oh_p4_matched[tt], Lp4Matched_lab);
   }

   for(int tt=0;tt<npionp_otherHemi;tt++)
   {
   compip_cost_href_true[tt] = Cal_HadronFrame(vpip_oh_p4_matched[tt], Lp4Matched_lab, pp4Matched_lab, compip_pt_href_true[tt]);
   compip_theta_GJ_true[tt] = Cal_HadronGJFrame(vpip_oh_p4_matched[tt], Lp4Matched_lab);
   }

   for(int tt=0;tt<npionm_otherHemi;tt++)
   {
   compim_cost_href_true[tt] = Cal_HadronFrame(vpim_oh_p4_matched[tt], Lp4Matched_lab, pp4Matched_lab, compim_pt_href_true[tt]);
   compim_theta_GJ_true[tt] = Cal_HadronGJFrame(vpim_oh_p4_matched[tt], Lp4Matched_lab);
   }


 for(int igam =0; igam<vphoton.size();igam++)
      {
        std::vector<Mdst_gamma>::const_iterator i = i_begin + vphoton[igam];
        Hep3Vector h(i->px(),i->py(),i->pz());
        const Mdst_gamma& gam=*i;
        int id=(int)gam.get_ID();
        double px=gam.px();
        double py=gam.py();
        double pz=gam.pz();
        //Mdst_ecl_aux &aux =eclaux_mgr(Panther_ID(gam.ecl().get_ID()));
        //if(gam.ecl().quality()!=0) continue;
        //double e9oe25 =aux.e9oe25();
        double gammaE=sqrt(px*px+py*py+pz*pz);
        HepLorentzVector phoVec(px,py,pz,gammaE);
        Hep3Vector photVec3(px,py,pz);
        Ephoton[igam] = gammaE;
        theta_photon[igam] = photVec3.theta();//in lab.
        phi_photon[igam] = photVec3.phi();//in lab.
	inv_lamgam[igam] = (phoVec + p4).m();
  }

  nphoton=vphoton.size();

  //find lambdac --> lambda pi, lambda pi pi0 lambda 3pi
  FindLambdaC findLC;
  
  int nLambdaC = findLC.find( p4, ipion );  
  m_nlambdac = nLambdaC;

  if(nLambdaC>99) cout<<"warning number of lambdac "<<nLambdaC<<", ipion "<<ipion.size()<<endl;
  for(int ll = 0; ll < nLambdaC;++ll)
  {
   Particle plc = findLC.getLambdac(ll); 
   m_lc_mass[ll] = plc.p().m();
   //m_lc_charge[ll] = plc.charge();
   m_lc_mode[ll] = findLC.getMode(ll);
   m_lc_charge[ll] = findLC.getCharge(ll);
   //cout<<"charge "<< plc.charge()<<", "<<findLC.getCharge(ll)<<endl;
   if(ll>98) break;
  }


  //find D from 9 golden mode, exclude tracks from Lambda
  FindD0Dp findD;
  int nD = findD.find( ipion, ikaon,  ip_position);
  int cntD=0;
  for(int dd = 0; dd < nD; ++dd)
  {
  Particle Dmeson = findD.getDmeson(dd);
  HepLorentzVector p4D = Dmeson.p();
  p4D.boost(kinematics::CMBoost);
  int Dhemiflag=0;
  if(kinematics::thrustDirCM.dot(p4D.vect())>0) Dhemiflag=1;
  else Dhemiflag=-1;
  //cout<<"L hemi "<<hemiflag<<" D hemi "<<Dhemiflag<<endl;
  if(Dhemiflag*hemiflag>0) continue;
  m_D_mass[cntD] = Dmeson.p().m();
  m_D_mode[cntD] = findD.getMode(dd);
  m_D_charge[cntD] = findD.getCharge(dd);
  cntD++;
  if(cntD>98) break;//to protect
  }
  m_nD= cntD;
 
 //cout<<"nD "<<nD<<endl; 
  
  m_nMatchedQuark = nMatchedQuark;
  m_matchedFlavor= MatchedFlavor;
  m_info.cost_truemom_recThr=cost_truemom_recThr;
  m_info.nMatch=nMatch;
  m_info.nMatch_b=nMatch_b;
  m_info.cosl_true=cosl_true;
  m_info.cost_true=cost_true;
  m_info.refdif_true=refdir.angle(refdir_true);
  m_cost_qq=cost_true_qq;
  m_info.cost_d_true=cost_d_true;
  m_info.mom_diff=mom_diff;//momentum difference from MCtruch after Lambda match
  m_info.angle_diff=angle_diff;//angle between MCtruch after Lambda match
  m_info.Lp4_match_phi=Lp4Matched_lab.vect().phi();
  m_info.Lp4_match_theta=Lp4Matched_lab.vect().theta();
  m_info.Prop4_match_phi=pp4Matched_lab.vect().phi();
  m_info.Prop4_match_theta=pp4Matched_lab.vect().theta();
  m_info.Pionp4_match_phi=pip4Matched_lab.vect().phi();
  m_info.Pionp4_match_theta=pip4Matched_lab.vect().theta();
  Lp4Matched_lab.boost(kinematics::CMBoost);
  m_info.z_mc = 2*Lp4Matched_lab.e()/kinematics::Q;
  m_info.pt_mc = Lp4Matched_lab.perp(kinematics::thrustDirCMS_MC);
  m_info.mother=mother;
  m_info.mothermother=mothermother;
  m_info.mother_b=mother_b;
  m_info.isthep_b=isthep_b;
  //m_info.mothermother_b=mothermother_b;
  m_info.weight=m_weight;
  m_info.drop=m_ifDrop;
  m_info.drop_zpt=m_ifDrop_zptcom;
  m_info.drop_zpt_b=m_ifDrop_zptcom_b;
  m_info.drop_zpt_c=m_ifDrop_zptcom_c;
  m_info.matchedTrkIndex=Lmatched_trkIndex;
  m_index = NLambda; //index of the lambda candiates

  if(m_rm_tree!=1)m_tree->Fill();

}//end of loop vee

if(NLambda<1) return;
if(NLambda>20) std::cout<<"warning, Num of Lambda candidate: "<<NLambda<<std::endl;
if(m_debug) std::cout<<" Num of Lambda candidate: "<<NLambda<<std::endl;
//double diffangle = kinematics::thrustDirCM.angle( truthThrustAxisCMS );

//for mixed event;
//if(eventstore) 
//{
int recordEE=-1;
for(int ee=0;ee<m_proton_mix.size();ee++)
{

if(m_prcharge_mix_tmp*m_prcharge_mix[ee]>0)
{
double thrangle = thrust_mix_tmp.angle(thrust_mix[ee]);
if((thrangle<20./180.*pi&&hemiflag_mix_tmp*hemiflag_mix[ee]>0)|| (thrangle>160./180.*pi&&hemiflag_mix_tmp*hemiflag_mix[ee]<0))
{
m_Lambda_mix = m_proton_mix[ee] + m_pion_mix;///alreday in c.m.s
if(m_Lambda_mix.m()>1.095&&m_Lambda_mix.m()<1.135) 
{

//m_prmag_mix = m_proton_mix[ee].vect().mag();
//m_pimag_mix = m_pion_mix.vect().mag();
m_prmag_mix = m_proton_mix_lab[ee].vect().mag();
m_pimag_mix = m_pion_mix_lab.vect().mag();
//m_oa_mix  = m_pion_mix.vect().angle(m_proton_mix[ee].vect());//open angle of proton and pion in c.m.s
m_lamag_mix = (m_pion_mix_lab+m_proton_mix_lab[ee]).vect().mag();
m_lamag_mix_cms = m_Lambda_mix.vect().mag();
float hrustProj=thrust_mix[ee].dot(m_Lambda_mix.vect());
Hep3Vector refdir_mix= thrust_mix[ee].cross(m_Lambda_mix.vect());
Hep3Vector refdir_mix_b = refdir_mix;
if(hrustProj<0) refdir_mix = -1 *refdir_mix;
Hep3Vector boostvec = m_Lambda_mix.boostVector();
m_proton_mix[ee].boost(-1*boostvec);
int kind_mix =0;
if(m_prcharge_mix[ee]>0) kind_mix=2;
if(m_prcharge_mix[ee]<0) kind_mix=3;

m_cost_mix = cos(m_proton_mix[ee].vect().angle(refdir_mix));
m_cost_mix_b = cos(m_proton_mix[ee].vect().angle(refdir_mix_b));
m_Lmass_mix= m_Lambda_mix.m();
m_z_mix = 2*m_Lambda_mix.e()/kinematics::Q; 
m_pt_mix = m_Lambda_mix.vect().perp(thrust_mix[ee]); 
m_kind_mix=kind_mix;
if(m_rm_mixtree!=1)m_mixtree->Fill();
recordEE=ee;
}
}
}
if(recordEE>-1) break;
if(recordEE>Msize) cout<<"somethig is wrong, recordEE "<<recordEE<<endl;
}//end of loop previous protons
//}

//push more protons 
if(m_proton_mix.size()<Msize)
{
m_proton_mix.push_back(m_proton_mix_tmp);
m_proton_mix_lab.push_back(m_proton_mix_lab_tmp);
m_prcharge_mix.push_back(m_prcharge_mix_tmp);
thrust_mix.push_back(thrust_mix_tmp);
hemiflag_mix.push_back(hemiflag_mix_tmp);
}
//replace one of the previous protons 
else if(recordEE>-1) {
m_proton_mix[recordEE]= m_proton_mix_tmp;
m_proton_mix_lab[recordEE]= m_proton_mix_lab_tmp;
m_prcharge_mix[recordEE]= m_prcharge_mix_tmp;
thrust_mix[recordEE] = thrust_mix_tmp;
hemiflag_mix[recordEE]= hemiflag_mix_tmp;
}
//eventstore=true;

  //m_tup_evt->column("NLambda",NLambda);
  //m_tup_evt->column("ngood",ngood);
  //m_tup_evt->column("evtNo",evtNo);
  //m_tup_evt->column("runNo",runNo);
  //m_tup_evt->column("expNo",expNo);
  //m_tup_evt->column("eler",eler);
  //m_tup_evt->column("eher",eher);
  //m_tup_evt->column("Flag",randomFlag);
  //m_tup_evt->column("thrx",kinematics::thrustDirCM.x());
  //m_tup_evt->column("thry",kinematics::thrustDirCM.y());
  //m_tup_evt->column("thrz",kinematics::thrustDirCM.z());
  //m_tup_evt->column("thrx_mc",truthThrustAxisCMS.x());
  //m_tup_evt->column("thry_mc",truthThrustAxisCMS.y());
  //m_tup_evt->column("thrz_mc",truthThrustAxisCMS.z());
  //m_tup_evt->column("thrdiff",diffangle);
  //m_tup_evt->column("thrust",kinematics::thrustMag);

  //m_tup_evt->dumpData();


//treeData->
//treeData->Fill();

  return;
}

void LambdaAna::term() {
// if(m_rm_mctree)m_mctree->Delete();
 m_file->Write();
 m_file->Close();
 //m_tree->Write();
 //m_mctree->Write();
 //m_file->Write();
  return;
}

void LambdaAna::disp_stat(const char*)
{
}

 void LambdaAna::getDrDz(Mdst_charged_Manager::iterator chr_it, int masshyp, double& dr, double& dz, double& refitPx, double& refitPy, double& refitPz)
  {
    Mdst_trk& mdsttrk = chr_it->trk();
    Mdst_trk_fit &mdsttrkfit=mdsttrk.mhyp(masshyp);
    HepPoint3D pivot(mdsttrkfit.pivot_x(),mdsttrkfit.pivot_y(),mdsttrkfit.pivot_z());

    HepVector a( 5, 0 );
    a[0] = mdsttrkfit.helix( 0 ); // helix parameters defined at the pivot
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
    //    cout <<"helix momentum: "<< helix.momentum()<<endl;
    refitPx=helix.momentum().x();
    refitPy=helix.momentum().y();
    refitPz=helix.momentum().z();
    //HepLorentzVector boostedVec(helix.momentum(),sqrt(helix.momentum().mag2()+m_mass*m_mass));
    dr  = helix.dr();
    dz  = helix.dz();
  }

 //sector selection for photons from Ami
  int LambdaAna::GetECLSector(double theta){

    const double fwd1 = 12.4;
    const double fwd2 = 31.4;
    const double barrel1 = 32.2;
    const double barrel2 = 128.7;
    const double bwd1 = 130.7;
    const double bwd2 = 155.1;

    //theta = theta*TMath::RadToDeg();
     theta = theta/3.14*180;

    if(theta<=fwd1) return 1;
    if(theta>fwd1 && theta<fwd2) return 2;
    else if(theta>=fwd2 && theta<=barrel1) return 3; //gap
    else if(theta>barrel1 && theta<barrel2) return 4;
    else if(theta>=barrel2 && theta<=bwd1) return 5;//gap
    else if(theta>bwd1 && theta<bwd2) return 6;
    else if(theta>=bwd2) return 7;

    return -1;
  }

// HepLorentzVector firstquark(0,0,0,0);

Hep3Vector LambdaAna::readMC( Hep3Vector &m_qqaxis, int &numquark)
{

   numquark =0;
   int numLambda_c=0;//recored number of Lambda_c in this event
   int numLambda_c_mode0=0; //lambda_c --> lambda pi
   int numSigma0=0;//recored number of sigma0 in this event
   int numLambda=0;//recored number of lambda in this event
   int numLambda_b=0;//recored number of lambda in this event
   vector<int>  VquarkID; VquarkID.clear();
   vector<Hep3Vector> VquarkVec; VquarkVec.clear();
   vector<float> VquarkAng; VquarkAng.clear();
    HepLorentzVector firstquark(0,0,0,0);
   Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();
   //general propety thrust... 
   vector<Hep3Vector> allMCParticlesLab; allMCParticlesLab.clear();
   vector<Hep3Vector> allMCParticlesCMS; allMCParticlesCMS.clear();
   HepLorentzVector  qmomentum(0,0,0,0);
   HepLorentzVector  qmomentum_lab(0,0,0,0);
   HepLorentzVector  qbarmomentum(0,0,0,0);
   HepLorentzVector  qbarmomentum_lab(0,0,0,0);
   if(m_debug) cout<<" how many particles in MC " <<gen_hep_Mgr.size()<<endl;
   for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
    {
    bool leafbool=false; bool mleaf=true;
    int cnt_lm=0;
    bool interestL=false;
    int cntproton =0, cntpion=0;
    Gen_hepevt_Manager& gen_hepevt_mgr = Gen_hepevt_Manager::get_manager();
    if(abs(gen_it->idhep())==3122)
    {
    int ndaug = (gen_it->daLast()-gen_it->daFirst()) + 1;
    for (int i = 0; i < ndaug; i++) {
        Panther_ID ID(gen_it->daFirst()+i);
     if(ID==0)
       {
         if(m_debug)
           cout <<"daughter ID is wrong!!!" <<endl;
         break;
       }
     Gen_hepevt& daughter = gen_hepevt_mgr(ID);
     if(abs(daughter.idhep())==211) { cntpion++; }
     if(abs(daughter.idhep())==2212){ cntproton++;}
     }
     if(cntpion==1 &&cntproton==1) interestL =true;
     }

    //if(interestL) numLambda_b++; 
    if(abs(gen_it->idhep())==3122&&gen_it->mother().idhep()==10022) numLambda_b++;
    //if(gen_it->isthep()<0) continue; //???
    if(abs(gen_it->idhep())==3122) numLambda++;
    
    if(abs(gen_it->idhep())==4122) 
    {
    numLambda_c++;
    int Lcpion=0;
    int LcLambda=0;
    int Lcdaug = (gen_it->daLast()-gen_it->daFirst()) + 1;
     for (int i = 0; i < Lcdaug; i++) {
     Panther_ID ID(gen_it->daFirst()+i);
     if(ID==0)
       {
         if(m_debug)
           cout <<"daughter ID is wrong!!!" <<endl;
         break;
      }
     Gen_hepevt& daughter = gen_hepevt_mgr(ID);
     if(abs(daughter.idhep())==211) { Lcpion++; }
     if(abs(daughter.idhep())==3122){ LcLambda++;}
     }
    //count lambda_c --> lambda pi, we don't require that lambda --> p pi because this br is known well
    if (Lcdaug==2&&Lcpion==1&&LcLambda==1) numLambda_c_mode0++;
    }
    if(abs(gen_it->idhep())==3212) numSigma0++;
    int ndaug = (gen_it->daLast()-gen_it->daFirst()) + 1;
    //cout<<"mcParticle: "<<gen_it->idhep()<<endl;
    bool isquark=false;

    HepLorentzVector v4tmp (gen_it->PX(),gen_it->PY(),gen_it->PZ(),gen_it->E()); 
    HepLorentzVector v4tmp_boost = v4tmp;
    v4tmp_boost.boost(kinematics::CMBoost);

    if(abs(gen_it->idhep())==1||abs(gen_it->idhep())==2||abs(gen_it->idhep())==3||abs(gen_it->idhep())==4) { isquark=true; }
    //if(isquark&&gen_it->mother().idhep()==10022)
    //HepLorentzVector firstquark(0,0,0,0);
    if(isquark)
    { 
     //HepLorentzVector tmp(gen_it->mother().PX(),gen_it->mother().PY(),gen_it->mother().PZ(),gen_it->mother().E());
     //tmp.boost(kinematics::CMBoost);
     //cout<<gen_it->idhep() <<" mother "<< gen_it->mother()<<" mother p4 "<< tmp.px()<<", "<<tmp.py()<<", "<<tmp.pz()<<", "<<tmp.e()<<endl;
     //cout<<" mother mass "<<tmp.m()<<endl;
     numquark++;
     if(numquark==1) firstquark=v4tmp_boost;
     if(gen_it->idhep()>0) { qmomentum =  v4tmp_boost;  qmomentum_lab = v4tmp; }
     else { qbarmomentum = v4tmp_boost;  qbarmomentum_lab = v4tmp; }
     VquarkID.push_back(gen_it->idhep());
     VquarkVec.push_back(v4tmp_boost.vect());
     VquarkAng.push_back(v4tmp_boost.vect().angle(firstquark.vect()));
     //cout<<" numquark "<<numquark << " px "<<firstquark.py()<<" angle "<<v4tmp_boost.vect().angle(firstquark.vect())<<endl;
    } 
   

    for(int ii=0;ii<7;ii++)
    {
   if(abs(gen_it->idhep())==leafParticleID[ii]) leafbool=true;
    //if(gen_it->mother()) cout<<" MC mother "<<abs(gen_it->mother().idhep())<<endl;
   if(gen_it->mother()&&abs(gen_it->mother().idhep())==leafParticleID[ii]) { mleaf=false;}
    }
   if(!(leafbool&&mleaf)) continue; 
   
   //HepLorentzVector v4tmp (gen_it->PX(),gen_it->PY(),gen_it->PZ(),gen_it->E()); 
   allMCParticlesLab.push_back(v4tmp.vect());
   allMCParticlesCMS.push_back(v4tmp_boost.vect());
    }

//  cout<<"numquark "<<numquark<<endl;

   m_qqaxis = qmomentum.vect();
   Hep3Vector  m_qqaxis_lab = qmomentum_lab.vect(); 
   if(rand() % 100 <50) { m_qqaxis= qbarmomentum.vect();  m_qqaxis_lab = qbarmomentum_lab.vect();} 
   kinematics::qqDirCMS_MC = m_qqaxis;
   kinematics::qqDirLab_MC = m_qqaxis_lab;

   if(m_debug) cout<<" how many leaf particles in MC " <<allMCParticlesCMS.size()<<endl;

   Thrust t_MC=thrustall(allMCParticlesCMS.begin(),allMCParticlesCMS.end(),retSelf);
   Thrust t_MC_lab=thrustall(allMCParticlesLab.begin(),allMCParticlesLab.end(),retSelf);
   kinematics::thrustDirCMS_MC=t_MC.axis;
   kinematics::thrustDirLab_MC=t_MC_lab.axis;
   if(rand() % 100 <50)
     {
      kinematics::thrustDirCMS_MC=-1*kinematics::thrustDirCMS_MC;
      kinematics::thrustDirLab_MC=-1*kinematics::thrustDirLab_MC;
     }
   double thrust_mc=t_MC.thru;
   int nleaf_mc= allMCParticlesCMS.size();

   if(m_debug) cout<<" thrust value in MC " <<thrust_mc<<endl;

   int cnt_Lambda_MC=0;
   int cnt_Lambda_decayToPPi_MC=0;
   int index=0;
   for( Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++,index++)
    {
     //cout<<"ii "<<ii<<endl;
     HepLorentzVector pp4_MC, pip4_MC;
     HepLorentzVector pp4_MC_CMS, pip4_MC_CMS;
     Gen_hepevt_Manager& gen_hepevt_mgr = Gen_hepevt_Manager::get_manager();
     //Gen_hepevt gParticle = Gen_hepevt(*gen_it);
     //if (gParticle.isthep()<0) continue;
     //if(gen_it->isthep()<0  ) continue;
     if(gen_it->idhep()!=3122&&gen_it->idhep()!=-3122) continue; 
     cnt_Lambda_MC++;
     int cnt_proton =0, cnt_pion=0;
     //find daughters
     int ndaug = (gen_it->daLast()-gen_it->daFirst()) + 1;
     if(m_debug) cout<<"Lambda ndaughters "<< ndaug <<endl;
     int pcharge=0;
     if(ndaug!=2) continue;
     int LambdaID = gen_it->get_ID();
     int ProtonID =-1; 
     int PionID =-1; 
     for (int i = 0; i < ndaug; i++) {
        Panther_ID ID(gen_it->daFirst()+i);
        if(m_debug) cout<<"Lambda ndaughters in list "<< ID <<endl;
       if(ID==0)
       {
         if(m_debug)
           cout <<"daughter ID is wrong!!!" <<endl;
         break;
       }
     Gen_hepevt& daughter = gen_hepevt_mgr(ID);
     if(m_debug) cout<<"Lambda daughters hepid "<< i << " :  "<<daughter.idhep() <<endl;
     if(abs(daughter.idhep())==211) { HepLorentzVector pip4(daughter.PX(),daughter.PY(),daughter.PZ(),daughter.E()); cnt_pion++; pip4_MC=pip4; PionID = ID;} 
     if(abs(daughter.idhep())==2212){
     HepLorentzVector pp4(daughter.PX(),daughter.PY(),daughter.PZ(),daughter.E()); cnt_proton++; pp4_MC=pp4; 
     if (daughter.idhep()>0) pcharge=1; else pcharge=-1;
     ProtonID =ID;
     }    
      }

   if(cnt_pion!=1 || cnt_proton!=1) continue;
   int m_Lmother=gen_it->mother().idhep();
   int m_Listhep=gen_it->isthep();
   int m_Lmothermother=-1;
   Gen_hepevt& Itmother=gen_hepevt_mgr(Panther_ID(gen_it->mother_ID()));
   if(Itmother.mother()!=NULL) m_Lmothermother=Itmother.mother().idhep();
   //if(m_Lmother==3312||m_Lmother==-3312) cout<<" gen_it->isthep()<0  "<<gen_it->isthep()<<endl;
   //cout<<"mother "<< m_Lmother<<"check "<<m_Lmothermother<<endl;
   //if(pcharge>0) cout<<" Lambda its mother "<< gen_it->mother().idhep() <<endl;
   //if(pcharge<0) cout<<" antiLambda its mother "<< gen_it->mother().idhep() <<endl;
   cnt_Lambda_decayToPPi_MC++;
   HepLorentzVector Lp4_lab(gen_it->PX(),gen_it->PY(),gen_it->PZ(),gen_it->E()); 
   HepLorentzVector Lp4_CMS = Lp4_lab; Lp4_CMS.boost(kinematics::CMBoost);
   double m_pmag=Lp4_lab.vect().mag();
   pp4_MC_CMS = pp4_MC; pp4_MC_CMS.boost(kinematics::CMBoost);
   pip4_MC_CMS = pip4_MC; pip4_MC_CMS.boost(kinematics::CMBoost);

    float m_Lm_mc=Lp4_lab.m();
    float m_Lm_check=(pp4_MC+pip4_MC).m();
    float m_z_mc=2*Lp4_CMS.e()/kinematics::Q;
    int zbinID=-1;
    for(int dd=0;dd<5;dd++)
    { if(m_z_mc>zbinRange[dd]&&m_z_mc<zbinRange[dd+1]){ zbinID=dd; break;}}
    if(zbinID<0) { continue; cout<<" zbinID "<<zbinID<<endl; }
    float hrustProj=kinematics::thrustDirCMS_MC.dot(Lp4_CMS.vect())/(kinematics::thrustDirCMS_MC.mag()*Lp4_CMS.vect().mag());
    double pt_mc= Lp4_CMS.vect().perp(kinematics::thrustDirCMS_MC);
    int ptbinID=-1;
    for(int bb=0;bb<5;bb++)
    { if(pt_mc>ptbinRange[bb]&&pt_mc<ptbinRange[bb+1]){ ptbinID=bb; break;}}
    if(ptbinID<0) { continue; cout<<" ptbinID "<<ptbinID<<endl; }

    double tangle = Lp4_CMS.vect().angle(kinematics::thrustDirCMS_MC);
    int hemiflag=0;
    if(hrustProj>0) hemiflag=1;
    else { hemiflag=-1;
    tangle = Lp4_CMS.vect().angle(-kinematics::thrustDirCMS_MC);
    }

    ///calculate the costheta_p 
    double lptheta_mc =-1;
    double thetap_mc_fake=-1;
    double costheta_mc_fake=-1;
    double costheta_mc_fake_b=-1;
    double costheta_mc_fake_c=-1;
    double costheta_mc_d=-1;
    double costheta_mc_fake_d=-1;
    double costheta_mc_test=-1;
    Hep3Vector refdir_mc(0,0,0);
    double thetap_mc_qq = CalTheta(m_qqaxis,m_qqaxis_lab,Lp4_lab,pp4_MC, pip4_MC, lptheta_mc, costheta_mc_test, thetap_mc_fake, costheta_mc_fake, costheta_mc_fake_b, costheta_mc_fake_c, costheta_mc_d, costheta_mc_fake_d, refdir_mc);
    double costheta_mc_qq = cos(thetap_mc_qq);
    double thetap_mc = CalTheta(kinematics::thrustDirCMS_MC,kinematics::thrustDirLab_MC,Lp4_lab,pp4_MC, pip4_MC, lptheta_mc, costheta_mc_test, thetap_mc_fake, costheta_mc_fake, costheta_mc_fake_b, costheta_mc_fake_c, costheta_mc_d, costheta_mc_fake_d, refdir_mc);
    double costheta_mc = cos(thetap_mc);
    double sintheta_mc = sin(thetap_mc);
    double cosl_mc = cos(lptheta_mc);
    double protonmag_mc= pp4_MC.vect().mag();//in the lab frame
    double pionmag_mc = pip4_MC.vect().mag();
    double protonTheta_mc= pp4_MC.vect().cosTheta();//in the lab frame
    double pionTheta_mc = pip4_MC.vect().cosTheta();
    double LTheta_mc = Lp4_lab.vect().cosTheta();
    double protonTheta_CM_mc= pp4_MC_CMS.vect().cosTheta();//in the cms frame
    double pionTheta_CM_mc = pip4_MC_CMS.vect().cosTheta();
    double LTheta_CM_mc = Lp4_CMS.vect().cosTheta();

    //Weighted input Pt=0.05;
    int m_drop=0;
    double decayPar=0;
    if(gen_it->idhep()>0) decayPar = m_LamDecayPar;
    else decayPar = m_antiLamDecayPar;

    //double random = r1.Uniform(0, 1+m_inputPolarization*fabs(decayPar));
    //double var = 1+decayPar*m_inputPolarization*costheta_mc;
    double random = r1.Uniform(0, 1+fabs(m_inputPolarizationArray[zbinID]*decayPar));
    double var = 1+decayPar*m_inputPolarizationArray[zbinID]*costheta_mc;
    if(random>var) {  m_drop=1;  weightedProton.push_back(ProtonID); weightedPion.push_back(PionID);}

    //random = r5.Uniform(0, 1+fabs(m_inputZPtcombinedArray[zbinID][ptbinID]*decayPar));
    //var = 1+decayPar*m_inputZPtcombinedArray[zbinID][ptbinID]*costheta_mc;
    int m_drop_zptcom=0;
    int m_drop_zptcom_b=0;
    int m_drop_zptcom_c=0;
    double inputlinear=-0.1*m_z_mc*pt_mc;
    random = r5.Uniform(0, 1+fabs(inputlinear*decayPar));
    var = 1+decayPar*inputlinear*costheta_mc;
    if(random>var) {  m_drop_zptcom=1;  weightedProton_zptcom.push_back(ProtonID); weightedPion_zptcom.push_back(PionID);}
    inputlinear=-0.1*m_z_mc*pt_mc*4;
    random = r8.Uniform(0, 1+fabs(inputlinear*decayPar));
    var = 1+decayPar*inputlinear*costheta_mc;
    if(random>var) {  m_drop_zptcom_b=1;  weightedProton_zptcom_b.push_back(ProtonID); weightedPion_zptcom_b.push_back(PionID);}
 
    random = r9.Uniform(0, 1+fabs(m_inputZPtcombinedArray[zbinID][ptbinID]*decayPar));
    var = 1+decayPar*m_inputZPtcombinedArray[zbinID][ptbinID]*costheta_mc;
    if(random>var) {  m_drop_zptcom_c=1;  weightedProton_zptcom_c.push_back(ProtonID); weightedPion_zptcom_c.push_back(PionID);}

    float weight_factor = weight_zpt[zbinID][ptbinID];
    weight_factor*=var;
    //weightedFactor_proton.push_back(ProtonID);
    //weightedFactor_pion.push_back(PionID);
    weightedFactor_lambdaId.push_back(LambdaID);
    weightedFactor_factor.push_back(weight_factor);
    //find a hadron(pion/kaon)  in the other hemi
   int hindex=0;
   int nkp_mc=0, nkm_mc=0, npip_mc=0, npim_mc=0; 
   for( Gen_hepevt_Manager::iterator gen_jt=gen_hep_Mgr.begin();gen_jt!=gen_hep_Mgr.end();gen_jt++,hindex++)
    {
     if(abs(gen_jt->idhep())!=211&&abs(gen_jt->idhep())!=321) continue; 
     HepLorentzVector p4(gen_jt->PX(), gen_jt->PY(),gen_jt->PZ(),gen_jt->E());
     HepLorentzVector p4_CMS = p4; p4_CMS.boost(kinematics::CMBoost);
     int hadron_hemiflag=0;
     float hadron_hrustProj=kinematics::thrustDirCMS_MC.dot(p4_CMS.vect())/(kinematics::thrustDirCMS_MC.mag()*p4_CMS.vect().mag());
     if(hadron_hrustProj>0) hadron_hemiflag=1; else { hadron_hemiflag=-1;}
     if(hemiflag*hadron_hemiflag>0) continue;
     int hID = gen_jt->get_ID();
     if(hID==ProtonID||hID==PionID) continue;
     double z_hadron_mc=2*p4_CMS.e()/kinematics::Q;
     int zhID=-1;
     for(int dd=0;dd<5;dd++)
     {if(z_hadron_mc>HadronZRange[dd]&&z_hadron_mc<HadronZRange[dd+1]) {zhID=dd; break;} }
     if(zhID<0) { continue; cout<<" zhID "<<zhID<<endl; }
     //if() continue;  //remove daughter of the Lambda

    //int m_drop_comh=0; //combine hadron thrust frame
    //int m_drop_comh_href=0; //combine hadron hadron frame
    //double random_com = r2.Uniform(0, 1+m_inputPolarization*fabs(decayPar));
    //double var_com = 1+decayPar*m_inputPolarization*costheta_mc;
    double inputPol = 0;
    //if(pcharge*(gen_jt->idhep())<0) inputPol = m_inputCombinedPolArray_B[zbinID][zhID]; 
    //else {inputPol = m_inputCombinedPolArray[zbinID][zhID]; }
    //double random_com = r2.Uniform(0, 1+fabs(inputPol*decayPar));
    //double var_com = 1+decayPar*inputPol*costheta_mc;
    //int hID = gen_jt->get_ID();
    //if(hID==ProtonID||hID==PionID) continue;
    //if(random_com>var_com) { m_drop_comh=1; comweightedProton.push_back(ProtonID); comweightedPion.push_back(PionID); comweighted_Hadron.push_back(hID); }
    float pt_href =0;
    double costheta_href = Cal_HadronFrame( p4_CMS, Lp4_lab,pp4_MC,pt_href );
    double theta_GJ = Cal_HadronGJFrame( p4_CMS, Lp4_lab );
    /*double random_com_href = r4.Uniform(0, 1+fabs(inputPol*decayPar));
    double var_com_href = 1+decayPar*inputPol*costheta_href;
    if(random_com_href>var_com_href) { m_drop_comh_href=1; comweightedProton_href.push_back(ProtonID); comweightedPion_href.push_back(PionID); comweighted_Hadron_href.push_back(hID); }
    */

    int ptbin_href_ID=-1;
    int m_drop_zptcomh_href=0;
    for(int bb=0;bb<5;bb++)
    {if(pt_href>ptbinRange_href[bb]&&pt_href<ptbinRange_href[bb+1]) {ptbin_href_ID=bb; break;} }
    if(ptbin_href_ID<0) { continue; cout<<" ptbin_href_ID "<<ptbin_href_ID<<endl; }
    //inputPol = m_inputZPtcombinedArray[zbinID][ptbin_href_ID];
    inputPol = -0.05*m_z_mc*pt_href;
    double random_zptcom_href = r6.Uniform(0, 1+fabs(inputPol*decayPar));
    double var_zptcom_href = 1+decayPar*inputPol*costheta_href;
    if(random_zptcom_href>var_zptcom_href) { m_drop_zptcomh_href=1; comweightedProton_zptcom_href.push_back(ProtonID); comweightedPion_zptcom_href.push_back(PionID); comweighted_Hadron_zptcom_href.push_back(hID); }

    int m_drop_zptcomh_href_b=0;
    inputPol = -0.05*m_z_mc*pt_href*4;
    random_zptcom_href = r7.Uniform(0, 1+fabs(inputPol*decayPar));
    var_zptcom_href = 1+decayPar*inputPol*costheta_href;
    if(random_zptcom_href>var_zptcom_href) { m_drop_zptcomh_href_b=1; comweightedProton_zptcom_href_b.push_back(ProtonID); comweightedPion_zptcom_href_b.push_back(PionID); comweighted_Hadron_zptcom_href_b.push_back(hID); }

     if((gen_jt->idhep())==321) 
     	{kpz_mc[nkp_mc]=z_hadron_mc; /*kp_cost_cms_mc[nkp_mc]=p4_CMS.vect().cosTheta(); comkp_drop[nkp_mc]=m_drop_comh;comkp_drop_href[nkp_mc]=m_drop_comh_href;*/comkp_zpt_drop_href[nkp_mc]=m_drop_zptcomh_href;comkp_zpt_drop_href_b[nkp_mc]=m_drop_zptcomh_href_b; comkp_cost_href_mc[nkp_mc]=costheta_href; comkp_pt_href_mc[nkp_mc]=pt_href;  comkp_theta_GJ_mc[nkp_mc]=theta_GJ; nkp_mc++; }
     if((gen_jt->idhep())==-321)
	 {kmz_mc[nkm_mc]=z_hadron_mc; /*km_cost_cms_mc[nkm_mc]=p4_CMS.vect().cosTheta(); comkm_drop[nkm_mc]=m_drop_comh;comkm_drop_href[nkm_mc]=m_drop_comh_href;*/ comkm_zpt_drop_href[nkm_mc]=m_drop_zptcomh_href;comkm_zpt_drop_href_b[nkm_mc]=m_drop_zptcomh_href_b;comkm_cost_href_mc[nkm_mc]=costheta_href; comkm_pt_href_mc[nkm_mc]=pt_href;  comkm_theta_GJ_mc[nkm_mc]=theta_GJ; nkm_mc++; }
     if((gen_jt->idhep())==211)  
	{pipz_mc[npip_mc]=z_hadron_mc; /*pip_cost_cms_mc[npip_mc]=p4_CMS.vect().cosTheta(); compip_drop[npip_mc]=m_drop_comh; compip_drop_href[npip_mc]=m_drop_comh_href;*/ compip_zpt_drop_href[npip_mc]=m_drop_zptcomh_href;compip_zpt_drop_href_b[npip_mc]=m_drop_zptcomh_href_b; compip_cost_href_mc[npip_mc]=costheta_href; compip_pt_href_mc[npip_mc]=pt_href; compip_theta_GJ_mc[npip_mc]=theta_GJ; npip_mc++; }
     if((gen_jt->idhep())==-211) 
	{pimz_mc[npim_mc]=z_hadron_mc; /*pim_cost_cms_mc[npim_mc]=p4_CMS.vect().cosTheta(); compim_drop[npim_mc]=m_drop_comh; compim_drop_href[npim_mc]=m_drop_comh_href;*/ compim_zpt_drop_href[npim_mc]=m_drop_zptcomh_href;compim_zpt_drop_href_b[npim_mc]=m_drop_zptcomh_href_b;compim_cost_href_mc[npim_mc]=costheta_href; compim_pt_href_mc[npim_mc]=pt_href; compim_theta_GJ_mc[npim_mc]=theta_GJ; npim_mc++; }

    }//end of loop hadrons

    m_mc_info.hflag=hemiflag;//hemiflag of Lambda
    m_mc_info.pcharge=pcharge;
    m_mc_info.Lmag_mc=m_pmag;
    m_mc_info.thr_mc=thrust_mc;
    //angle_qqthr=m_qqaxis.angle(thrustDirCMS_MC);
    angle_qqthr = m_qqaxis.angle(kinematics::thrustDirCMS_MC);
    if(m_qqaxis.dot(kinematics::thrustDirCMS_MC)<0) angle_qqthr = m_qqaxis.angle(-1 * kinematics::thrustDirCMS_MC);;
    m_mc_info.Nleaf=nleaf_mc;
    m_mc_info.tan_mc=tangle;
    m_mc_info.drop=m_drop; // alone Lambda or anti-Lambda
    m_mc_info.drop_zpt=m_drop_zptcom; // alone Lambda or anti-Lambda
    m_mc_info.drop_zpt_b=m_drop_zptcom_b; // alone Lambda or anti-Lambda
    m_mc_info.drop_zpt_c=m_drop_zptcom_c; // alone Lambda or anti-Lambda
    m_mc_info.costheta=costheta_mc;
    m_costheta_mc_qq = costheta_mc_qq;
    m_mc_info.sintheta=sintheta_mc;
    m_mc_info.cosl_mc=cosl_mc;
    m_mc_info.thetap_fake=thetap_mc_fake;
    m_mc_info.costheta_fake=costheta_mc_fake;
    m_mc_info.costheta_fake_b=costheta_mc_fake_b;
    m_mc_info.costheta_fake_c=costheta_mc_fake_c;
    m_mc_info.costheta_d=costheta_mc_d;
    m_mc_info.costheta_fake_d=costheta_mc_fake_d;
    //m_mc_info.lp_new=cosl_mc_new;
    //m_mc_info.lp_check=cosl_mc_check;
    //m_mc_info.tp_test=cost_mc_test;
    //m_mc_info.tp_new=cosl_mc_new;
    m_mc_info.thep_mc=thetap_mc;
    m_mc_info.z_mc=m_z_mc;
    m_mc_info.Lmother=m_Lmother;
    m_mc_info.Listhep=m_Listhep;
    m_mc_info.Lmothermother=m_Lmothermother;
    m_mc_info.Lm_mc=m_Lm_mc;
    m_mc_info.Lm_ch=m_Lm_check;
    m_mc_info.pt_mc=pt_mc;
    m_mc_info.prmag_mc=protonmag_mc;
    m_mc_info.pimag_mc=pionmag_mc;
    m_mc_info.prthe_mc=protonTheta_mc;
    m_mc_info.pithe_mc=pionTheta_mc;
    m_mc_info.lthe_mc=LTheta_mc;
    m_mc_info.prthe_cm=protonTheta_CM_mc;
    m_mc_info.pithe_cm=pionTheta_CM_mc;
    m_mc_info.lthe_cm=LTheta_CM_mc;
    m_mc_info.nkp_mc=nkp_mc;//in the opposite hemisphere
    m_mc_info.nkm_mc=nkm_mc;//in the opposite hemisphere
    m_mc_info.npip_mc=npip_mc;
    m_mc_info.npim_mc=npim_mc;
    m_mc_info.nquark = numquark;
    m_weight_mc = weight_factor;
    for(int hh = 0; hh<numquark; hh++)
    {
    int tmphemi = 0; 
    if(VquarkVec[hh].angle(kinematics::thrustDirCMS_MC)<0.5*pi) tmphemi=1;
    else tmphemi= -1;
    quarkhemi[hh] =  tmphemi;
    quarkId[hh] = VquarkID[hh];
    m_quarkAng[hh] = VquarkAng[hh];
    //cout<<" angle "<<hh<<" , "<<m_quarkAng[hh]<<endl;
    }

    for(int mm=0;mm<nkp_mc;mm++)   comkp_cost_mc[mm] = costheta_mc;
    for(int mm=0;mm<nkm_mc;mm++)   comkm_cost_mc[mm] = costheta_mc;
    for(int mm=0;mm<npip_mc;mm++)  compip_cost_mc[mm] = costheta_mc;
    for(int mm=0;mm<npim_mc;mm++)  compim_cost_mc[mm] = costheta_mc;

    if(m_rm_mctree!=1) m_mctree->Fill();
    
   } 

  //m_tup_mc_evt->column("NLam",cnt_Lambda_decayToPPi_MC);
  //m_tup_mc_evt->column("thr_mc",thrust_mc);
  //m_tup_mc_evt->dumpData();

 m_nLambda_b=numLambda_b;
 m_nLambda=numLambda;
 m_nLc=numLambda_c;
 m_nLc_mode0=numLambda_c_mode0;
 m_nSigma0=numSigma0;
 m_evttree->Fill();

return kinematics::thrustDirCMS_MC; 


}

void LambdaAna::other(int*, BelleEvent*, int*)
{

}

double LambdaAna::CalTheta(Hep3Vector thrustDirCMS, Hep3Vector thrustDirLab,  HepLorentzVector Lambda_p4_Lab , HepLorentzVector Proton_p4_Lab, HepLorentzVector Pion_p4_Lab, double &lptheta, double &costheta_test, double &theta_fake, double &cost_fake, double &cost_fake_b, double &cost_fake_c, double &cost_d, double &cost_fake_d, Hep3Vector &refdir_record)
{
    HepLorentzVector Lambda_p4_CMS = Lambda_p4_Lab;
    Lambda_p4_CMS.boost(kinematics::CMBoost);
    if(m_debug) cout<<"CalTheta:: CMS boost vector "<<kinematics::CMBoost<<endl;
    float hrustProj=thrustDirCMS.dot(Lambda_p4_CMS.vect())/(thrustDirCMS.mag()*Lambda_p4_CMS.vect().mag());
    //double pt= Lambda_p4_CMS.vect().perp(thrustDirCMS);

    ///calculate the costheta_p 
    //the reference direction
    //Hep3Vector refdir =  Lambda_p4_CMS.vect().cross(thrustDirCMS);//.cross(Lambda_p4_CMS.vect());
    //if(hrustProj<0) refdir = -1*Lambda_p4_CMS.vect().cross(thrustDirCMS);

    Hep3Vector refdirLab = Lambda_p4_Lab.vect().cross(thrustDirLab);//.cross(Lambda_p4_Lab.vect());
    if(hrustProj<0) refdirLab = -1*Lambda_p4_Lab.vect().cross(thrustDirLab);
    HepLorentzVector refdir_test_c(1,refdirLab.x()/refdirLab.mag(),refdirLab.y()/refdirLab.mag(),refdirLab.z()/refdirLab.mag());

     Hep3Vector refdir = thrustDirCMS.cross(Lambda_p4_CMS.vect());
     if(hrustProj<0)  refdir = -1*thrustDirCMS.cross(Lambda_p4_CMS.vect());
     refdir_record = refdir;
     HepLorentzVector refdir_test_b(1,refdir.x()/refdir.mag(),refdir.y()/refdir.mag(),refdir.z()/refdir.mag());

    // Hep3Vector refdir_fake = thrustDirCMS;
     //if(hrustProj<0)  refdir_fake= -1*thrustDirCMS;

     //double alpha=r3.Uniform(0, m_PI);
     //double beta=r3.Uniform(0, 2*m_PI);
     //Hep3Vector refdir_fake(sin(alpha)*sin(beta),sin(alpha)*cos(beta),cos(alpha));
     //if(hrustProj<0)  refdir_fake= -1*refdir_fake;

     //Hep3Vector refdir_tmp = thrustDirCMS.cross(Lambda_p4_CMS.vect());//randomly
     //Hep3Vector refdir_fake_b = Lambda_p4_CMS.vect().cross(refdir_tmp);
     Hep3Vector refdir_fake_b = Lambda_p4_CMS.vect().cross(refdir);
     //if(hrustProj<0)  refdir_fake_b= -1*refdir_fake_b;

     Hep3Vector refdir_fake= Lambda_p4_CMS.vect().cross(refdir);
     HepLorentzVector refdir_fake_tmp(1, refdir_fake.x()/refdir_fake.mag(),refdir_fake.y()/refdir_fake.mag(),refdir_fake.z()/refdir_fake.mag()); 

    // Hep3Vector refdir_fake_c = refdir_fake.cross(refdir);
     //if(hrustProj<0)  refdir_fake_c= -1*refdir_fake_c;

     
     //Hep3Vector beamdir(0,0,1);/// ???
     Hep3Vector beamdir(0,0,0);
     beamdir= kinematics::firstElectronCM.vect();
     //cout<<"beam dir in cms: "<<beamdir.x()<<", " <<beamdir.y()<<" , "<<beamdir.z()<<endl;
     //
     Hep3Vector refdir_d = beamdir.cross(Lambda_p4_CMS.vect());
     if(beamdir.dot(Lambda_p4_CMS.vect())<0)  refdir_d = -1*beamdir.cross(Lambda_p4_CMS.vect());

     Hep3Vector refdir_fake_d = Lambda_p4_CMS.vect().cross(refdir_d);
     //if(beamdir.dot(Lambda_p4_CMS.vect())<0)  refdir_fake_d = -1*refdir_d.cross(Lambda_p4_CMS.vect());

     ///asume direction of reference  
    //HepLorentzVector refdirNew(refdir,refdir.mag());
    //boost proton and pion into the CMS of Lambda
    Hep3Vector Lboostvector =  Lambda_p4_Lab.boostVector();
    Hep3Vector Lboostvector_2 =  Lambda_p4_CMS.boostVector();
    HepLorentzVector Lambda_boost = Lambda_p4_Lab;
    HepLorentzVector Lambda_CMS_boost = Lambda_p4_CMS;
    HepLorentzVector p4Proton_boost = Proton_p4_Lab;
    HepLorentzVector p4Proton_CMS= Proton_p4_Lab;
    p4Proton_CMS.boost(kinematics::CMBoost);
    //HepLorentzVector p4PionLab_boost = Pion_p4_Lab;
    Lambda_boost.boost(-Lboostvector);
    //if(m_debug) cout<<"Lambda in Lambda CMS "<< Lambda_boost<<endl;
    Lambda_CMS_boost.boost(-Lboostvector_2);
    //if(m_debug) cout<<"Lambda in Lambda CMS "<< Lambda_CMS_boost<<endl;
    p4Proton_boost.boost(-Lboostvector);
    HepLorentzVector p4Proton_CMS_boost = p4Proton_CMS;
    p4Proton_CMS_boost.boost(-Lboostvector_2);
    refdir_test_b.boost(-Lboostvector_2);//to Lambda rest frame
    refdir_test_c.boost(-Lboostvector);//to Lambda rest frame
    //cout<<"proton in Lambda rest frame "<<Proton_p4_Lab<<", "<<p4Proton_boost<<", "<<p4Proton_CMS_boost<<endl;
    refdir_fake_tmp.boost(-Lboostvector_2);
    //p4PionLab_boost.boost(-Lboostvector);
    //refdirNew.boost(-Lboostvector_2);// to Lambda rest frame
    //double costheta = p4Proton_boost.vect().dot(refdir)/refdir.mag()/p4Proton_boost.vect().mag();
    HepLorentzVector thrusttmp(1,thrustDirCMS.x()/thrustDirCMS.mag(),thrustDirCMS.y()/thrustDirCMS.mag(),thrustDirCMS.z()/thrustDirCMS.mag());
    thrusttmp.boost(-Lboostvector_2);//boost from CMS to Lambda rest frame;
    HepLorentzVector beamdirtmp(1,beamdir.x()/beamdir.mag(),beamdir.y()/beamdir.mag(),beamdir.z()/beamdir.mag());
    beamdirtmp.boost(-Lboostvector_2);
    Hep3Vector refdir_test=thrusttmp.vect().cross(-1*beamdirtmp.vect());
    if(hrustProj<0)  refdir_test = thrusttmp.vect().cross(beamdirtmp.vect());
    double thetap = p4Proton_boost.vect().angle(refdir);
    //costheta_test = cos(p4Proton_boost.vect().angle(refdir_test));
    //costheta_test = cos(p4Proton_CMS_boost.vect().angle(refdir_test_b.vect()));
    //costheta_test = cos(p4Proton_CMS_boost.vect().angle(refdir_test_c.vect()));
    //cout<<"test "<< costheta_test <<", "<<cos(p4Proton_CMS_boost.vect().angle(refdir_test))<<endl;
    costheta_test = cos(p4Proton_boost.vect().angle(refdirLab));
    theta_fake= p4Proton_boost.vect().angle(refdir_fake_tmp.vect());
    double  theta_fake_b= p4Proton_boost.vect().angle(refdir_fake_b);
    //double  theta_fake_c= p4Proton_boost.vect().angle(refdir_fake_c);
    double  theta_d= p4Proton_boost.vect().angle(refdir_d);
    double  theta_fake_d= p4Proton_boost.vect().angle(refdir_fake_d);
    cost_fake=cos(theta_fake);
    cost_fake_b=cos(theta_fake_b);
    //cost_fake_c=cos(theta_fake_c);
    cost_d=cos(theta_d);
    cost_fake_d=cos(theta_fake_d);
    lptheta = p4Proton_boost.vect().angle(Lambda_p4_CMS.vect());
    //lptheta_new = p4Proton_boost.vect().angle(Lambda_p4_Lab.vect());//longitudinal
    //lptheta_check = p4Proton_CMS_boost.vect().angle(Lambda_p4_CMS.vect());//longitudinal wrong
    //tptheta = p4Proton_boost.vect().angle(refdirNew.vect());//transverse test
    //tptheta_new = p4Proton_boost.vect().angle(refdirLab);//transverse new check 
    return thetap;
}

double LambdaAna::Cal_HadronFrame(HepLorentzVector refHadron_p4_CMS, HepLorentzVector Lambda_p4_Lab, HepLorentzVector Proton_p4_Lab, float & pt_href)
{
 HepLorentzVector Lambda_p4_CMS = Lambda_p4_Lab;
 Lambda_p4_CMS.boost(kinematics::CMBoost);
 //Hep3Vector hrdir = refHadron_p4_CMS.vect().cross(Lambda_p4_CMS.vect());
 //if(refHadron_p4_CMS.vect().dot(Lambda_p4_CMS.vect())<0)  hrdir = -1*refHadron_p4_CMS.vect().cross(Lambda_p4_CMS.vect());
 Hep3Vector hrdir = -1*refHadron_p4_CMS.vect().cross(Lambda_p4_CMS.vect());

 pt_href = Lambda_p4_CMS.vect().perp(refHadron_p4_CMS.vect());

  HepLorentzVector p4Proton_boost = Proton_p4_Lab;
  Hep3Vector Lboostvector =  Lambda_p4_Lab.boostVector();
  p4Proton_boost.boost(-Lboostvector);

 double costheta_href = p4Proton_boost.vect().angle(hrdir);

 return cos(costheta_href);

}

double LambdaAna::Cal_HadronGJFrame(HepLorentzVector refHadron_p4_CMS, HepLorentzVector Lambda_p4_Lab)
{
 HepLorentzVector Lambda_p4_CMS = Lambda_p4_Lab;
 Lambda_p4_CMS.boost(kinematics::CMBoost);
 Hep3Vector lhplane = kinematics::firstElectronCM.vect().cross(refHadron_p4_CMS.vect());
 //if(refHadron_p4_CMS.vect().dot(kinematics::firstElectronCM.vect())<0)  lhplane = -1*refHadron_p4_CMS.vect().cross(kinematics::firstElectronCM.vect());

 Hep3Vector hplane = refHadron_p4_CMS.vect().cross(Lambda_p4_CMS.vect());
 //if(refHadron_p4_CMS.vect().dot(Lambda_p4_CMS.vect())<0)  hplane = -1*refHadron_p4_CMS.vect().cross(Lambda_p4_CMS.vect());

 double theta_GJ =hplane.angle(lhplane);
 Hep3Vector tmp = lhplane.cross(hplane);
 double sign = refHadron_p4_CMS.angle(tmp) ;
 if(sign>(3.1415926/2)) theta_GJ*=-1;
 theta_GJ*=2;
 if(theta_GJ<-1*3.1415926)  theta_GJ =  theta_GJ + 2*3.1415926;
 else if (theta_GJ>3.1415926) theta_GJ = theta_GJ - 2*3.1415926;

//return cos(theta_GJ);
return theta_GJ;
//return lhplane.cosTheta();
//return refHadron_p4_CMS.vect().angle(kinematics::firstElectronCM.vect());
//return refHadron_p4_CMS.vect().cosTheta();
//return hplane.theta();
//return lhplane.phi();
//return hplane.phi();

}

void LambdaAna::GetNhadron(int LheimFlag,Vp4 vp4kaon,Vp4 vp4pion, Vp4 vp4kaon_lab, Vp4 vp4pion_lab, Vint kaonchr, Vint pionchr, int LprotonTrkID, int LpionTrkID, Vint pionTrkID,Vint kaonTrkID, Vfloat& kpz, Vfloat& kmz, Vfloat& pipz, Vfloat& pimz, Vfloat& kp_oat, Vfloat& km_oat, Vfloat& pip_oat, Vfloat& pim_oat, Vint &kp_nMatch, Vint &km_nMatch, Vint &pip_nMatch, Vint &pim_nMatch, Vint &kp_trackIndex, Vint &km_trackIndex, Vint &pip_trackIndex, Vint &pim_trackIndex, Vfloat &vkp_mom_diff, Vfloat &vkm_mom_diff,Vfloat &vpip_mom_diff,Vfloat &vpim_mom_diff,Vfloat &vkp_angle_diff, Vfloat &vkm_angle_diff, Vfloat &vpip_angle_diff, Vfloat &vpim_angle_diff, Vp4 &vkp_oh_p4, Vp4 &vkm_oh_p4,Vp4 &vpip_oh_p4,Vp4 &vpim_oh_p4, Vp4 &vkp_oh_p4_matched, Vp4 &vkm_oh_p4_matched,Vp4 &vpip_oh_p4_matched,Vp4 &vpim_oh_p4_matched, Vint& kp_matchedType, Vint& km_matchedType, Vint& pip_matchedType, Vint& pim_matchedType, Vint& Vkaonindex, Vint& Vpionindex)
{
kpz.clear();
kmz.clear();
pipz.clear();
pimz.clear();
kp_oat.clear();
km_oat.clear();
pip_oat.clear();
pim_oat.clear();
Vkaonindex.clear();
Vpionindex.clear();
if(LheimFlag!=1&&LheimFlag!=-1) { cout<<"somthing is wrong hemi of the Lambda is not decided "<<endl; }
if(LpionTrkID<0) { cout<<"somthing is wrong TrkID of pion "<< LpionTrkID<<endl; return;}
if(LprotonTrkID<0) { cout<<"somthing is wrong TrkID of proton "<< LprotonTrkID<<endl; return;}

for(int kk=0;kk<vp4kaon.size();kk++)
{
if(LpionTrkID==kaonTrkID[kk]) continue;//remove the track from Lambda candidate
if(LprotonTrkID==kaonTrkID[kk]) continue;//remove the track from Lambda candidate
float hrustProj=kinematics::thrustDirCM.dot(vp4kaon[kk].vect());
//if(hrustProj>0) hemiflag=1;
float opangle;
if(LheimFlag==-1) opangle = vp4kaon[kk].vect().angle(kinematics::thrustDirCM);//using the reversed thrust axis
else opangle = vp4kaon[kk].vect().angle(-1*kinematics::thrustDirCM);

if(hrustProj*LheimFlag<0) 
{
double kaonz = 2*vp4kaon[kk].e()/kinematics::Q;
if(kaonz<m_zhcut || kaonz > m_zhcut_up) continue;
Vkaonindex.push_back(kk);
int nMatch=0, trackIndex=-1;
double mom_diff=100; double angle_diff=100;
HepLorentzVector p4_cms_match(0,0,0,0);
int htype=-1;
if(m_mc==1)TrackMatch(vp4kaon_lab[kk], pkaon[kk], nMatch, trackIndex,mom_diff, angle_diff, p4_cms_match, htype);
if(kaonchr[kk]>0) {  kpz.push_back(kaonz); kp_oat.push_back(opangle); kp_nMatch.push_back(nMatch); kp_trackIndex.push_back(trackIndex); vkp_mom_diff.push_back(mom_diff); vkp_angle_diff.push_back(angle_diff); vkp_oh_p4.push_back(vp4kaon[kk]); vkp_oh_p4_matched.push_back(p4_cms_match); kp_matchedType.push_back(htype);
}//nkaonp_otherHemi++;
else { kmz.push_back(kaonz); km_oat.push_back(opangle); km_nMatch.push_back(nMatch); km_trackIndex.push_back(trackIndex); vkm_mom_diff.push_back(mom_diff); vkm_angle_diff.push_back(angle_diff);vkm_oh_p4.push_back(vp4kaon[kk]);vkm_oh_p4_matched.push_back(p4_cms_match); km_matchedType.push_back(htype);}
}
}
///pions

for(int pp=0;pp<vp4pion.size();pp++)
{
if(LpionTrkID==pionTrkID[pp]) continue;//remove the track from Lambda candidate
if(LprotonTrkID==pionTrkID[pp]) continue;
float hrustProj=kinematics::thrustDirCM.dot(vp4pion[pp].vect());
float opangle;
if(LheimFlag==-1) opangle = vp4pion[pp].vect().angle(kinematics::thrustDirCM);//using the reversed thrust axis
else opangle = vp4pion[pp].vect().angle(-1*kinematics::thrustDirCM);
if(hrustProj*LheimFlag<0)
{
float pionz = 2*vp4pion[pp].e()/kinematics::Q;
if(pionz<m_zhcut || pionz > m_zhcut_up) continue;
Vpionindex.push_back(pp);
int nMatch=0, trackIndex=-1;
HepLorentzVector p4_cms_match(0,0,0,0);
double mom_diff=100; double angle_diff=100;
int htype=-1;
if(m_mc==1)TrackMatch(vp4pion_lab[pp], ppion[pp], nMatch, trackIndex,mom_diff,angle_diff, p4_cms_match, htype);
if(pionchr[pp]>0) { pipz.push_back(pionz); pip_oat.push_back(opangle); pip_nMatch.push_back(nMatch); pip_trackIndex.push_back(trackIndex); vpip_mom_diff.push_back(mom_diff); vpip_angle_diff.push_back(angle_diff); vpip_oh_p4.push_back(vp4pion[pp]); vpip_oh_p4_matched.push_back(p4_cms_match); pip_matchedType.push_back(htype);}
else { pimz.push_back(pionz); pim_oat.push_back(opangle); pim_nMatch.push_back(nMatch); pim_trackIndex.push_back(trackIndex); vpim_mom_diff.push_back(mom_diff); vpim_angle_diff.push_back(angle_diff);vpim_oh_p4.push_back(vp4pion[pp]);vpim_oh_p4_matched.push_back(p4_cms_match); pim_matchedType.push_back(htype);
//cout<<"charge check "<<pionchr[pp]<<endl;
}
}
}

}

int LambdaAna::TrackMatch( HepLorentzVector p4_rec_lab, Particle p, int &nMatch, int& trackIndex, double &mom_diff, double &angle_diff, HepLorentzVector &p4_matched_cms, int& htype)
{
//nMatch=0; 
trackIndex=-1;
mom_diff=10000;
angle_diff=10000;
int index=0;
nMatch=-1;

//const Gen_hepevt& hep_lp(gen_level(get_hepevt(p.mdstCharged())));
//const Gen_hepevt& hep(get_hepevt(mdstCharged()));
const Gen_hepevt &h_hadron = gen_level(get_hepevt(p.mdstCharged()));
//if (abs(h_hadron.idhep())!=211 && abs(h_hadron.idhep())!=321) { nMatch=0; return 0; }
if(h_hadron)
{
htype = h_hadron.idhep();
trackIndex= h_hadron.get_ID();
nMatch=1;
HepLorentzVector p4_lab(h_hadron.PX(),h_hadron.PY(),h_hadron.PZ(),h_hadron.E());
mom_diff= p4_rec_lab.vect().mag()-p4_lab.vect().mag();
//if(abs(mom_diff)>2)cout<<" gen. "<< p4_lab.vect().mag() <<"  rec. "<< p4_rec_lab.vect().mag() <<endl;
angle_diff=p4_rec_lab.vect().angle(p4_lab.vect());
p4_matched_cms=p4_lab;
p4_matched_cms.boost(kinematics::CMBoost);
return 1;
}
return 0;
}

int LambdaAna::MatchToMCTruth_new(int pcharge, int& mother, int &mothermother, int &isthep, HepLorentzVector Lp4_rec_lab,  Particle proton, Particle pion, Particle Lambda, double &cosl_true, double &cost_true, double &cost_true_qq, double &cost_d_true, int &Lindex, int &Prindex, int &Piindex, double &mom_diff, double &angle_diff, HepLorentzVector& Lp4_lab_true, HepLorentzVector& pp4_lab_true,HepLorentzVector& pip4_lab_true, Hep3Vector &refdir_true, int& numMatchquark, int& MatchFalvor)
{
Lindex=-1;
Prindex=-1;
Piindex=-1;

  int nMatch=0; 
  angle_diff=100000;
  mom_diff=100000;
  //const Gen_hepevt &h_proton=proton.genHepevt();
  //const Gen_hepevt &h_pion=pion.genHepevt();
  const Gen_hepevt &h_proton = gen_level(get_hepevt(proton.mdstCharged()));
  const Gen_hepevt &h_pion = gen_level(get_hepevt(pion.mdstCharged()));
  Gen_hepevt_Manager& gen_hepevt_mgr = Gen_hepevt_Manager::get_manager();

  if(h_proton && h_pion && h_proton.mother() && h_pion.mother() && h_proton.mother().get_ID()==h_pion.mother().get_ID() &&abs(h_pion.mother().idhep())==3122) 
    { 
     Gen_hepevt &h_Lambda= h_proton.mother();
    //if(abs(h_Lambda.idhep())!=3122) {/* cout<<" mother of proton/pion is " << h_Lambda.idhep()<<endl;*/ continue; }
     Lambda.relation().genHepevt(h_proton.mother()); 
     HepLorentzVector Lp4_lab(h_Lambda.PX(),h_Lambda.PY(),h_Lambda.PZ(),h_Lambda.E());
     HepLorentzVector pp4_lab(h_proton.PX(),h_proton.PY(),h_proton.PZ(),h_proton.E());
     HepLorentzVector pip4_lab(h_pion.PX(),h_pion.PY(),h_pion.PZ(),h_pion.E());
     mom_diff = Lp4_rec_lab.vect().mag()-Lp4_lab.vect().mag();
     angle_diff = Lp4_rec_lab.vect().angle(Lp4_lab.vect());
     nMatch=1;
     Lindex=h_Lambda.get_ID();
     Prindex=h_proton.get_ID();
     Piindex=h_pion.get_ID();
     mother = h_Lambda.mother().idhep();
     isthep = h_Lambda.isthep();
     Gen_hepevt& Itmother=gen_hepevt_mgr(Panther_ID(h_Lambda.mother_ID()));
     if(Itmother.mother()!=NULL) mothermother=Itmother.mother().idhep();
     else mothermother=-1;
     //cout<<"mother "<< mother<<" grandmother "<<mothermother<<endl;
     //mothermother = h_Lambda.mother().mother().idhep();

 //cout<<" Lindex "<< Lindex <<" Prindex "<<Prindex<<endl;

//calculate angles using MCTruth
  double costheta_test=-1;
 double thetap_fake=-1; double costheta_fake=-1;  double costheta_fake_b=-1;  double costheta_fake_c=-1;
 double costheta_d=-1;  double costheta_fake_d=-1;
 Hep3Vector refdir_mc;
 double thetap_mc_qq = CalTheta(kinematics::qqDirCMS_MC,kinematics::qqDirLab_MC,Lp4_lab, pp4_lab, pip4_lab, cosl_true,costheta_test, thetap_fake, costheta_fake, costheta_fake_b,  costheta_fake_c, costheta_d, costheta_fake_d, refdir_mc);
 double thetap_mc = CalTheta(kinematics::thrustDirCMS_MC,kinematics::thrustDirLab_MC,Lp4_lab, pp4_lab, pip4_lab, cosl_true,costheta_test, thetap_fake, costheta_fake, costheta_fake_b,  costheta_fake_c, costheta_d, costheta_fake_d, refdir_mc);

 refdir_true=refdir_mc;
 cosl_true = cos(cosl_true);
 cost_true = cos(thetap_mc);
 cost_true_qq = cos(thetap_mc_qq);
 cost_d_true = costheta_d;
 Lp4_lab_true=Lp4_lab;
 pp4_lab_true=pp4_lab;
 pip4_lab_true=pip4_lab;


//get the quark flavor that is in same hameishpere with Lambda
vector<Hep3Vector> VquarkVecCM; VquarkVecCM.clear();
Vint VquarkFlavor; VquarkFlavor.clear();
Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();
for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
{
bool isquark=false;
if(abs(gen_it->idhep())==1||abs(gen_it->idhep())==2||abs(gen_it->idhep())==3||abs(gen_it->idhep())==4) { isquark=true; }
if(isquark)
{
  HepLorentzVector v4tmp_boost (gen_it->PX(),gen_it->PY(),gen_it->PZ(),gen_it->E());
  v4tmp_boost.boost(kinematics::CMBoost);
  VquarkVecCM.push_back(v4tmp_boost.vect());
  VquarkFlavor.push_back(gen_it->idhep());
}
}
HepLorentzVector Lp4_CM = Lp4_lab;
Lp4_CM.boost(kinematics::CMBoost);
int lambdaHemei =0;
if(Lp4_CM.vect().angle(kinematics::thrustDirCMS_MC)<0.5*pi) lambdaHemei =1;
else lambdaHemei =-1;
 numMatchquark=0;
 MatchFalvor = 0;
for(int hh = 0; hh<VquarkFlavor.size(); hh++)
    {
    int quarkhemi= 0;
    if(VquarkVecCM[hh].angle(kinematics::thrustDirCMS_MC)<0.5*pi) quarkhemi=1;
    else quarkhemi= -1;
    if(lambdaHemei*quarkhemi<0) continue;
    numMatchquark++;
    MatchFalvor=VquarkFlavor[hh];
    }
}
return nMatch;

}

//try to match to Xi
int LambdaAna::MatchToMCTruth(int pcharge, int& mother, int &mothermother, int& isthep, HepLorentzVector pp4, HepLorentzVector pip4, HepLorentzVector &pp4Matched_lab, HepLorentzVector &pip4Matched_lab, HepLorentzVector& Lp4Matched_lab, int &trkindex, double &mom_diff, double &angle_diff)
{
trkindex=-1;
Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();

  int nMatch=0;
  mother=0;
  isthep=0;
  angle_diff=100000;
  mom_diff=100000;
  int index=0;
  for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++,index++)
    {
     HepLorentzVector pp4_MC, pip4_MC;
     int cnt_npion=0; int cnt_nproton=0;
      Gen_hepevt_Manager& gen_hepevt_mgr = Gen_hepevt_Manager::get_manager();
     if(gen_it->idhep()!=3122&&gen_it->idhep()!=-3122) continue;
     int ndaug = (gen_it->daLast()-gen_it->daFirst()) + 1;
     if(m_debug) cout<<"Lambda ndaughters "<< ndaug <<endl;
     int pcharge_MC=0;
     for (int i = 0; i < ndaug; i++) {
        Panther_ID ID(gen_it->daFirst()+i);
        if(m_debug) cout<<"Lambda ndaughters in list "<< ID <<endl;
       if(ID==0)
       {
         if(m_debug)
           cout <<"daughter ID is wrong!!!" <<endl;
         break;
       }
     Gen_hepevt& daughter = gen_hepevt_mgr(ID);
     if(m_debug) cout<<"Lambda daughters hepid "<< i << " :  "<<daughter.idhep() <<endl;
     if(abs(daughter.idhep())==211) {
      HepLorentzVector pip4_MC_tmp(daughter.PX(),daughter.PY(),daughter.PZ(),daughter.E());
      pip4_MC = pip4_MC_tmp;
      cnt_npion++;
      }
      if(abs(daughter.idhep())==2212){
     HepLorentzVector pp4_MC_tmp(daughter.PX(),daughter.PY(),daughter.PZ(),daughter.E());
      pp4_MC = pp4_MC_tmp;
      cnt_nproton++;
     if (daughter.idhep()>0) pcharge_MC=1; else pcharge_MC=-1;
     }
     }

  if(ndaug!=2) continue;
  if(cnt_nproton!=1||cnt_npion!=1) { if(m_debug)cout<<" cannot pass number of pront, pion "<<cnt_nproton<<" ,  "<<cnt_npion<<endl; continue;}
  if(pcharge_MC*pcharge<0) { if(m_debug)cout<<" cannot pass charge  of proton "<<pcharge_MC<<" , "<<pcharge<<endl; continue; }
  if(pp4.vect().angle(pp4_MC.vect())>20/180.*3.1415926) continue;
  if(pip4.vect().angle(pip4_MC.vect())>20/180.*3.1415926) continue;
  if(fabs(pip4.vect().mag()-pip4_MC.vect().mag())>0.2) continue;
  if(fabs(pp4.vect().mag()-pp4_MC.vect().mag())>0.2) continue;
   HepLorentzVector Lp4_lab(gen_it->PX(),gen_it->PY(),gen_it->PZ(),gen_it->E());
   HepLorentzVector Lp4_rec_lab = pp4 +pip4;
   nMatch++;

 if(Lp4_rec_lab.vect().angle(Lp4_lab.vect())<angle_diff)
   {
   if(m_debug) cout<<"find one matched Lambda "<<endl;
   mom_diff = Lp4_rec_lab.vect().mag()-Lp4_lab.vect().mag();
   angle_diff = Lp4_rec_lab.vect().angle(Lp4_lab.vect());
   pp4Matched_lab = pp4_MC;
   pip4Matched_lab = pip4_MC;
   Lp4Matched_lab = Lp4_lab;
   mother= gen_it->mother().idhep();
   isthep= gen_it->isthep();
   //mothermother= (gen_it->mother().mother()).idhep();
   trkindex = index;
   }
 }
 if(nMatch<1) return 0;
 return nMatch;
}

/*genhep_vec* LambdaAna::getDaughters(const Gen_hepevt &mother)
{

   Gen_hepevt_Manager& gen_hepevt_mgr = Gen_hepevt_Manager::get_manager();

   int n_children = mother.daLast() - mother.daFirst() + 1;

   genhep_vec *children = new genhep_vec();

   for(int i=0; i<n_children; i++) {

     Panther_ID ID0(mother.daFirst()+i);
     if(ID0==0)
       {
         if(m_debug)
           cout <<"wrong!!!" <<endl;
         break;
       }
     Gen_hepevt& temp = gen_hepevt_mgr(ID0);

     if (temp)
   {
         children->push_back(&temp);
   }
   }

   return children;
 }
*/


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
