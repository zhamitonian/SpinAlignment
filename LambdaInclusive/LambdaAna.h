#include "belle.h"
#include <cmath>
//#include "TH1F.h"
#include "particle/Particle.h"
#include "particle/utility.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "tables/evtcls.h"
#include "basf/module_descr.h"
#include "TFile.h"
#include "TH1F.h"
#include <TTree.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TRandom3.h>

#include <panther/panther.h>
#include "tables/evtcls.h"
//#include "toolbox/Thrust.h"
#include "toolbox/FuncPtr.h"
//#include "TMath.h"

#include BELLETDF_H
#include HEPEVT_H
#include MDST_H


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<float> Vfloat;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<Hep3Vector> Vp3;
typedef std::vector<Particle> VParticle;

class LambdaAna : public Module {
public:
  LambdaAna();
  ~LambdaAna() {}
  void init(int *dummy);
  void term();
  void disp_stat(const char*);
  void hist_def(void);
  void event(BelleEvent*, int*);
  void begin_run(BelleEvent*, int *);
  void end_run(BelleEvent*, int *);
  void other(int*, BelleEvent*, int*);
  Hep3Vector readMC( Hep3Vector &m_qqaxis, int &numquark);
  void getDrDz(Mdst_charged_Manager::iterator chr_it, int masshyp, double& dr, double& dz, double& refitPx, double& refitPy, double& refitPz);

  double CalTheta(Hep3Vector thrustDirCM,Hep3Vector thrustDirLab, HepLorentzVector Lambda_p4_Lab , HepLorentzVector Proton_p4_Lab, HepLorentzVector Pion_p4_Lab, double &lptheta,double &costheta_test, double &theta_fake, double &cost_fake, double &cost_fake_b, double &cost_fake_c, double &cost_d, double &cost_fake_d, Hep3Vector &refdir);

  double Cal_HadronFrame( HepLorentzVector refHadron_p4_CMS, HepLorentzVector Lambda_p4_Lab, HepLorentzVector Proton_p4_Lab, float & pt_href);
  double Cal_HadronGJFrame(HepLorentzVector refHadron_p4_CMS, HepLorentzVector Lambda_p4_Lab);
 void GetNhadron(int LheimFlag, Vp4 vp4kaon, Vp4 vp4pion, Vp4 vp4kaon_lab, Vp4 vp4pion_lab,Vint kaonchr, Vint pionchr,int LprotonTrkID, int LpionTrkID,Vint pionTrkID, Vint kaonTrkID, Vfloat& kpz, Vfloat& kmz, Vfloat& pipz, Vfloat& pimz, Vfloat& kp_oat, Vfloat& km_oat, Vfloat& pip_oat, Vfloat& pim_oat, Vint &kp_nMatch, Vint &km_nMatch, Vint &pip_nMatch, Vint &pim_nMatch, Vint &kp_trackIndex, Vint &km_trackIndex, Vint &pip_trackIndex, Vint &pim_trackIndex, Vfloat &vkp_mom_diff, Vfloat &vkm_mom_diff, Vfloat &vpip_mom_diff, Vfloat &vpim_mom_diff, Vfloat &vkp_angle_diff,  Vfloat &vkm_angle_diff, Vfloat &vpip_angle_diff, Vfloat &vpim_angle_diff, Vp4 &vkp_oh_p4, Vp4 &vkm_oh_p4, Vp4 &vpip_oh_p4, Vp4 &vpim_oh_p4, Vp4 &vkp_oh_p4_matched, Vp4 &vkm_oh_p4_matched, Vp4 &vpip_oh_p4_matched, Vp4 &vpim_oh_p4_matched, Vint& kp_matchedType, Vint& km_matchedType, Vint& pip_matchedType, Vint& pim_matchedType, Vint& Vkaonindex, Vint& Vpionindex);

  int GetECLSector(double theta);

  //int MatchToMCTruth(int pcharge, int &mother, HepLorentzVector pp4, HepLorentzVector pip4, HepLorentzVector &pp4Matched_lab, HepLorentzVector &pip4Matched_lab, HepLorentzVector&  Lp4Matched_lab, double &cosl_true, double &cost_true, int &matchedTrkIndex, double &mom_diff, double &angle_diff);
  int MatchToMCTruth(int pcharge, int& mother, int &mothermother, int& isthep, HepLorentzVector pp4, HepLorentzVector pip4, HepLorentzVector &pp4Matched_lab, HepLorentzVector &pip4Matched_lab, HepLorentzVector& Lp4Matched_lab, int &trkindex, double &mom_diff, double &angle_diff);

  int MatchToMCTruth_new(int pcharge, int &mother, int &mothermother, int& isthep, HepLorentzVector p4_lab ,Particle proton, Particle pion, Particle Lambda, double &cosl_true, double &cost_true, double &cost_true_qq, double &cost_d_true, int &Lmatched_trkIndex, int &Prmatched_trkIndex, int &Pimatched_trkIndex, double &mom_diff, double &angle_diff, HepLorentzVector&  Lp4Matched_lab, HepLorentzVector &pp4Matched_lab, HepLorentzVector &pip4Matched_lab, Hep3Vector &ref_dir, int& numMatchquark, int& MatchFalvor);

  //int TrackMatch(HepLorentzVector p4_lab, int charge, int type, int &nMatch, int &trackIndex, double &mom_diff, double &angle_diff);
  int TrackMatch(HepLorentzVector p4_lab, Particle p, int &nMatch, int &trackIndex, double &mom_diff, double &angle_diff, HepLorentzVector& p4_cms_matched, int& htype);
  //genhep_vec* getDaughters(const Gen_hepevt &mother);
public: // BASF parameter

   bool m_debug;
   int m_mc;
   bool m_PID;
   int m_rm_mctree;
   int m_rm_tree;
   int m_rm_mixtree;
   double m_thrustMagCut;
   double m_zhcut;
   double m_zhcut_up;
   char *m_output_filename;
   char *m_output_filename_genTree;
   double m_inputPolarization;
   //double m_inputPolarizationArray[5];
   //double m_inputCombinedPolArray[5][5];
   //double m_inputCombinedPolArray_B[5][5];
   //double zbinRange[6];
   //double HadronZRange[6];
   double m_LamDecayPar;
   double m_antiLamDecayPar;

private:
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

   //struct m_mix{
   float m_Lmass_mix;
   float m_pimag_mix;
   //float m_oa_mix;
   float m_prmag_mix;
   float m_lamag_mix;
   float m_lamag_mix_cms;
   float m_pt_mix;
   float m_z_mix;
   float m_cost_mix;
   float m_cost_mix_b;
   int m_kind_mix;
   //}
   

   int countEvt;
   Vint weightedProton;
   Vint weightedPion;

   Vint weightedProton_zptcom;
   Vint weightedPion_zptcom;
   Vint weightedProton_zptcom_b;//enlarge
   Vint weightedPion_zptcom_b;//enlarge
   Vint weightedProton_zptcom_c;
   Vint weightedPion_zptcom_c;

   Vint weightedFactor_lambdaId;
   Vfloat weightedFactor_factor;

   Vint comweightedProton;
   Vint comweightedPion;
   Vint comweighted_Hadron;
   Vint comweightedProton_href;
   Vint comweightedPion_href;
   Vint comweighted_Hadron_href;
   Vint comweightedProton_zptcom_href;
   Vint comweightedPion_zptcom_href;
   Vint comweighted_Hadron_zptcom_href;
   Vint comweightedProton_zptcom_href_b;//enlarge
   Vint comweightedPion_zptcom_href_b;//enlarge
   Vint comweighted_Hadron_zptcom_href_b;//enlarge
   VParticle pkaon;
   VParticle ppion;
   TRandom3 r1;
   //TRandom3 r2;
   //TRandom3 r3;
   //TRandom3 r4;
   TRandom3 r5;
   TRandom3 r6;
   TRandom3 r7;
   TRandom3 r8;
   TRandom3 r9;
   double m_PI;

   //BelleTuple *m_tup;
   //BelleTuple *m_tup_evt;
   //BelleTuple *m_tup_mc;
   //BelleTuple *m_tup_mc_evt;
   //BelleHistogram *hmass;
   //BelleHistogram *hthr;
   TFile* m_file;
   TTree* m_tree;
   TTree* m_mixtree;
   TTree* m_mctree;
   TTree* m_evttree;
    int m_nLambda_b;
    int m_nLambda;
    int m_nLc;
    int m_nLc_mode0;
    int m_nSigma0;

   struct lambda_info{
	int  evtNo;
	int  runNo;
	int  expNo;
	float Q;
        float visEnergy;//calculate by myself
        float Evis;//calculated by thrust 
	int ngood;
	float eler;//beam
	float eher;
	float thrust;
	float thrust_dir[3];//direction of thrust axis
	float thrust_dir_mc[3];//direction of thrust axis from MC truth
	float qq_axis_mc[3];//direction of thrust axis from MC truth
        float thrust_phi;
        float thrust_theta;
        float thrust_phi_lab;
        float thrust_theta_lab;
        float thrust_mc_phi;
        float thrust_mc_theta;
        float qq_axis_phi;
        float qq_axis_theta;
        float thrust_mc_phi_lab;
        float thrust_mc_theta_lab;
	float tangle;//angle between lambda and thrust axis 
        float thrdiff;//differcen between thrust axis and truth axis 
        float throff;//differcen between thrust axis and original q-qbar axis 
        float thrshift;
        int   numquark;
        //int   numu;
        //int   numd;
        //int   nums;
        //int   numc;
        float refdirdiff;//angle between reference direction cal. from rec. thrust axis and mc. thrust axis

	float Lmass;
	float rawmass;
	float masscheck;
	float z;
	float pt;
	float z_mc;//true value in MCTruth
	float pt_mcThr;
	float pt_mc;
	int   hemiflag;
	int   kind;//Lamda or anti-Lambda
	int   Isgood;//type from FindL
	float  chisq;//chisq of vertex fit
	float Lp4[4];//Lambda in CMS frame
	float Prop4[4];//proton in CMS frame
	float Piop4[4];//pion in CMS frame
	float Lp4_lab[4];//Lambda in Lab frame
	float Prop4_lab[4];//proton in Lab frame
	float Piop4_lab[4];//pion in Lab frame
        float Lp4_lab_phi;
        float Lp4_lab_theta;
        float Prop4_lab_phi;
        float Prop4_lab_theta;
        float Pionp4_lab_phi;
        float Pionp4_lab_theta;
        float Lp4_raw_phi;
        float Lp4_raw_theta;
        float Prop4_raw_phi;
        float Prop4_raw_theta;
        float Pionp4_raw_phi;
        float Pionp4_raw_theta;
        float Lp4_match_phi;
        float Lp4_match_theta;
        float Prop4_match_phi;
        float Prop4_match_theta;
        float Pionp4_match_phi;
        float Pionp4_match_theta;
	float prmag;//momentum of proton in Lab frame
	float pimag;//momentum of pion in Lab frame
	float lamag;//momentum of L in Lab frame
	float prtheta;//costheta of proton in Lab frame
	float pitheta;//costheta of pion in Lab frame
	float latheta;//costheta of L in Lab frame
	float prtheta_cms;//costheta of proton in cms 
	float pitheta_cms;//costheta of pion in cms 
	float lacostheta_cms;//costheta of L in cms 
	float y;//
	float latheta_cms;//theta of L in cms 
	float thetap;//transversed angle
	float costheta;//transversed angle
	float sintheta;//transversed angle
	float costheta_mcThrust;
        float cost_truemom_recThr;
	float cosl;//longitudinal 
  float weight;
        int    drop; //drop this lambda if match to the weighed Lambda in MCTruth
        //int    drop_pt;
        int    drop_zpt;
        int    drop_zpt_b;//enlarge
        int    drop_zpt_c;//enlarge
        int    matchedTrkIndex; 
	float thetap_fake;//transversed angle
	float costheta_test;//
	float costheta_fake;//
	float costheta_fake_b;//
	float costheta_fake_b_mcT;//
	float costheta_fake_c;//
	float costheta_d;//
	float costheta_fake_d;
	
	int  nkpoh;//number of kaon+ in the opposite hemisiphere
	int  nkmoh;//number of kaon- in the opposite hemisiphere
	int  npipoh;//number of pion+ in the opposite hemisiphere
	int  npimoh;//number of pion- in the opposite hemisiphere

	int nMatch;//match to MCTruth
	int mother;//mother in MCTruth
	int isthep;//mother in MCTruth
	int mothermother;//grandmother in MCTruth
	int nMatch_b;//check
	int mother_b;//search Xi
	int isthep_b;//search Xi
	int mothermother_b;
	float cosl_true;//match to MCTruth
	float cost_true;//match to MCTruth
	float refdif_true;//match to MCTruth
	float cost_d_true;//match to MCTruth
	float mom_diff;//match to MCTruth
	float angle_diff;//match to MCTruth
 
     } m_info;

        int m_nMatchedQuark;
        int m_matchedFlavor;
        float m_prob_proton;
        float m_prob_pk;
        float m_fl;
        //float m_flerr;
        float m_kp_prob[30];
        float m_km_prob[30];
        float kpz[30];
        float kmz[30];
        float pipz[30];
        float pimz[30];
        float m_kp_theta[30];
        float m_km_theta[30];
        float m_pip_theta[30];
        float m_pim_theta[30];
        float kp_oat[30];
        float km_oat[30];
        float pip_oat[30];
        float pim_oat[30];
        float kp_p4[30][4];//kp in the opposite hemi, p4 in the e+e- cms
        float km_p4[30][4];//
        float pip_p4[30][4];//
        float pim_p4[30][4];//

        int kp_nmatch[30];
        //int kp_drop[30];
        int km_nmatch[30];
        //int km_drop[30];
        int pip_nmatch[30];
        //int pip_drop[30];
        int pim_nmatch[30];
        float kpz_true[30];
        float kmz_true[30];
        float pipz_true[30];
        float pimz_true[30];
        int m_kp_matchedType[30];
        int m_km_matchedType[30];
        int m_pip_matchedType[30];
        int m_pim_matchedType[30];
        //int kp_drop[30];
        //int pim_drop[30];
        //int km_drop_href[30]; //repalced by zpt combined input
        //int kp_drop_href[30];
        //int pip_drop_href[30];
        //int pim_drop_href[30];
        int km_drop_zpt_href[30];
        int kp_drop_zpt_href[30];
        int pip_drop_zpt_href[30];
        int pim_drop_zpt_href[30];
        int km_drop_zpt_href_b[30];
        int kp_drop_zpt_href_b[30];
        int pip_drop_zpt_href_b[30];
        int pim_drop_zpt_href_b[30];
        float kp_mom_diff[30];
        float km_mom_diff[30];
        float pip_mom_diff[30];
        float pim_mom_diff[30];
        float kp_angle_diff[30];
        float km_angle_diff[30];
        float pip_angle_diff[30];
        float pim_angle_diff[30];
        float kp_oa_lam_cms[30];
        float km_oa_lam_cms[30];
        float pip_oa_lam_cms[30];
        float pim_oa_lam_cms[30];
        float kp_lam_mass[30];
        float km_lam_mass[30];
        float pip_lam_mass[30];
        float pim_lam_mass[30];
        float comkp_cost[30];
        float comkm_cost[30];
        float compip_cost[30];
        float compim_cost[30];
        float comkp_cost_href[30];
        float comkm_cost_href[30];
        float compip_cost_href[30];
        float compim_cost_href[30];
        float comkp_pt_href[30];
        float comkm_pt_href[30];
        float compip_pt_href[30];
        float compim_pt_href[30];
        float comkp_cost_href_true[30];
        float comkm_cost_href_true[30];
        float compip_cost_href_true[30];
        float compim_cost_href_true[30];
        float comkp_pt_href_true[30];
        float comkm_pt_href_true[30];
        float compip_pt_href_true[30];
        float compim_pt_href_true[30];
        float comkp_theta_GJ[30];
        float comkm_theta_GJ[30];
        float compip_theta_GJ[30];
        float compim_theta_GJ[30];
        float comkp_theta_GJ_true[30];
        float comkm_theta_GJ_true[30];
        float compip_theta_GJ_true[30];
        float compim_theta_GJ_true[30];

        float m_cost_qq;
				int m_index;
        int   nphoton;
        float Ephoton[50];
	float theta_photon[50];
	float phi_photon[50];
	float inv_lamgam[50];
	//float helicity[50];

        int m_nlambdac; 
        float m_lc_mass[100];
        int m_lc_charge[100];
        int m_lc_mode[100];
        int m_nD; 
        float m_D_mass[100];
        int m_D_charge[100];
        int m_D_mode[100];

   float eler;
   float eher;

   struct mc_info{ 
   float thr_mc;
   float thep_mc;
   float z_mc;
   int  Lmother;
   int  Listhep;
   int  Lmothermother;
   int Nleaf;
   float Lm_mc;
   float Lm_ch;
   float pt_mc;
   float costheta;
   float sintheta;
   float cosl_mc;
   int pcharge;
   int  hflag;
   float Lmag_mc;
   float pimag_mc;
   float prmag_mc;
   float pithe_mc;
   float lthe_mc;
   float prthe_mc;
   float pithe_cm;
   float prthe_cm;
   float lthe_cm;
   float tan_mc;
   int drop;
   int drop_pt;
   int drop_zpt;
   int drop_zpt_b;
   int drop_zpt_c;
   int nkp_mc;//number of k+ in the opposite hemi
   int nkm_mc;
   int npip_mc;
   int npim_mc;
   float thetap_fake;
   float costheta_fake;
   float costheta_fake_b;
   float costheta_fake_c;
   float costheta_d;
   float costheta_fake_d;
   
   int nquark;
   //int quarkhemi[nquark];/// 1 same with lambda -1 opposite hemi
   //int quarkId[nquark]; // u 1, d 2, s 3, c 4

   }m_mc_info;

   float m_weight_mc;
   float angle_qqthr;
   float m_costheta_mc_qq; 
   
  int quarkhemi[20];
  int quarkId[20];
  float m_quarkAng[20];
  float kpz_mc[30];
  float kmz_mc[30];
  float pipz_mc[30];
  float pimz_mc[30];
  //float kp_cost_cms_mc[30];
  //float km_cost_cms_mc[30];
  //float pip_cost_cms_mc[30];
  //float pim_cost_cms_mc[30];
  float comkp_cost_mc[30];
  float comkm_cost_mc[30];
  float compip_cost_mc[30];
  float compim_cost_mc[30];
  float comkp_cost_href_mc[30];
  float comkm_cost_href_mc[30];
  float compip_cost_href_mc[30];
  float compim_cost_href_mc[30];
  float comkp_pt_href_mc[30];
  float comkm_pt_href_mc[30];
  float compip_pt_href_mc[30];
  float compim_pt_href_mc[30];
  float comkp_theta_GJ_mc[30];
  float comkm_theta_GJ_mc[30];
  float compip_theta_GJ_mc[30];
  float compim_theta_GJ_mc[30];
  //int comkp_drop[30];
  //int comkm_drop[30];
  //int compip_drop[30];
  //int compim_drop[30];

  //int comkp_drop_href[30];
  //int comkm_drop_href[30];
  //int compip_drop_href[30];
  //int compim_drop_href[30];

  int comkp_zpt_drop_href[30];
  int comkm_zpt_drop_href[30];
  int compip_zpt_drop_href[30];
  int compim_zpt_drop_href[30];

  int comkp_zpt_drop_href_b[30];
  int comkm_zpt_drop_href_b[30];
  int compip_zpt_drop_href_b[30];
  int compim_zpt_drop_href_b[30];

///save data
   //double Lmass;


};

extern "C" Module_descr *mdcl_LambdaAna() {
  LambdaAna *module = new LambdaAna;
  Module_descr *dscr = new Module_descr( "LambdaAna", module );
  dscr->define_param( "output_filename", "output root file name",255 ,(char*)module->m_output_filename );
  //dscr->define_param( "isMCSample", "is data or MC ",(int)module->m_mc);
  dscr->define_param("isMCSample", "is data or MC ", &module->m_mc);
  dscr->define_param("rmMCTree", "keep the mcTree or not ", &module->m_rm_mctree);
  dscr->define_param("rmTree", "keep the Tree or not ", &module->m_rm_tree);
  dscr->define_param("rmMixTree", "keep the mixed event Tree or not ", &module->m_rm_mixtree);
  //dscr->define_param( "isMCSample", "is data or MC ",256,(char*)module->m_mc);
  //dscr->define_param( "output_filename_mcTree", "output root file name",255 ,(char*)module->m_output_filename_genTree );
  //dscr->define_param( "Debug", "open the debug mode",256,(int)module->m_debug);
  // dscr->define_param ( "my_filename", "filename", 256,
  //                      module->my_filename );
  return dscr;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
