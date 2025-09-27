//*-----
//  File : EventClassify.cc
//  Description : event classification routines. 
//
//  Date  : 10 - Feb - 1999
//  Author : Ichiro Adachi, Physics Division, KEK
//*-----

// for user's include file.
#include "belle.h"
#include MDST_H
#include EVTCLS_H
#include BELLETDF_H
#include TRG_H
#include EVTVTX_H
#include MDST_OBS_H

#include "trg/TrgBit.h"
#include "event/BelleEvent.h"

//*---
// add on '99-Aug-20. for gamma+phi selection. IA
//*---
#include "eid/eid.h"

//*---
// add on '99-Oct-01. for gamma+phi selection. IA
//*---
#include "kid/atc_pid.h"


// for thrust calculation.
#include "toolbox/Thrust.h"
#include "toolbox/FoxWolfr.h"
#include "toolbox/FuncPtr.h"

//*---
// add on '00-Mar-12. for HadronB selection from Brendan. IA.
//*---
//
// This is not a class header. Just code
//
#include "evtcls/hadsel.h"

#include "evtcls/EventClassify.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif



//
// static constants. 
//
static const double ELER(3.5);
static const double EHER(7.9965);
static const double theta(0.022);

//
// user's function for sorting.
//   descending order.
//
#if defined(__GNUG__)
int SortByEnergy( const HepLorentzVector **a,  
		  const HepLorentzVector **b )
{
#else
extern "C" int SortByEnergy( const void *av,  
			     const void *bv )
{
  const HepLorentzVector **a((const HepLorentzVector **)av);
  const HepLorentzVector **b((const HepLorentzVector **)bv); 
#endif
  if( (*a)->e() < (*b)->e() ) return 1;
  else if( (*a)->e() == (*b)->e() ) return 0;
  else return -1;
}

//   descending order.
#if defined(__GNUG__)
int SortByP( const Mdst_charged **a,  
	     const Mdst_charged **b )
{
#else
extern "C" int SortByP( const void *av,  
			const void *bv )
{
  const Mdst_charged **a((const Mdst_charged **)av);
  const Mdst_charged **b((const Mdst_charged **)bv);
#endif
  Hep3Vector v0( (*a)->px(), (*a)->py(), (*a)->pz() );
  Hep3Vector v1( (*b)->px(), (*b)->py(), (*b)->pz() );
  if( v0.mag() < v1.mag() ) return 1;
  else if( v0.mag() == v1.mag() ) return 0;
  else return -1;
}

//   descending order.
#if defined(__GNUG__)
int SortByPrest( const Mdst_charged **a,  
		 const Mdst_charged **b )
{
#else
extern "C" int SortByPrest( const void *av,  
			    const void *bv )
{
  const Mdst_charged **a((const Mdst_charged **)av);
  const Mdst_charged **b((const Mdst_charged **)bv);
#endif
  Hep3Vector v0( (*a)->px(), (*a)->py(), (*a)->pz() );
  Hep3Vector v1( (*b)->px(), (*b)->py(), (*b)->pz() );
  HepLorentzVector v0rest4( v0, sqrt(v0.mag2()+(*a)->mass()*(*a)->mass()) );
  HepLorentzVector v1rest4( v1, sqrt(v1.mag2()+(*b)->mass()*(*b)->mass()) );

  //  double ELER =3.5;
  //  double EHER =7.9965;
  //  double theta=0.022;
  HepLorentzVector CMrest(-EHER*sin(theta), 0., -EHER*cos(theta)+ELER,
			  EHER+ELER );
  Hep3Vector CMboostV=CMrest.boostVector();

  v0rest4.boost(CMboostV);
  v1rest4.boost(CMboostV);
  Hep3Vector v0rest3 = v0rest4.vect();
  Hep3Vector v1rest3 = v1rest4.vect();
  if( v0rest3.mag() < v1rest3.mag() ) return 1;
  else if( v0rest3.mag() == v1rest3.mag() ) return 0;
  else return -1;
}

#if defined(__GNUG__)
int SortByE( const Mdst_gamma **a,  
	     const Mdst_gamma **b )
{
#else
extern "C" int SortByE( const void *av,  
			const void *bv )
{
  const Mdst_gamma **a((const Mdst_gamma**)av);
  const Mdst_gamma **b((const Mdst_gamma**)bv);
#endif
  Hep3Vector v0( (*a)->px(), (*a)->py(), (*a)->pz() );
  Hep3Vector v1( (*b)->px(), (*b)->py(), (*b)->pz() );
  if( v0.mag() < v1.mag() ) return 1;
  else if( v0.mag() == v1.mag() ) return 0;
  else return -1;
}

//   descending order.
#if defined(__GNUG__)
int SortByErest( const Mdst_gamma **a,  
		 const Mdst_gamma **b )
{
#else
extern "C" int SortByErest( const void *av,  
			    const void *bv )
{
  const Mdst_gamma **a((const Mdst_gamma**)av);
  const Mdst_gamma **b((const Mdst_gamma**)bv);
#endif
  Hep3Vector v0( (*a)->px(), (*a)->py(), (*a)->pz() );
  Hep3Vector v1( (*b)->px(), (*b)->py(), (*b)->pz() );
  HepLorentzVector v0rest4( v0, v0.mag() );
  HepLorentzVector v1rest4( v1, v1.mag() );

  //  double ELER =3.5;
  //  double EHER =7.9965;
  //  double theta=0.022;
  HepLorentzVector CMrest(-EHER*sin(theta), 0., -EHER*cos(theta)+ELER,
			  EHER+ELER );
  Hep3Vector CMboostV=CMrest.boostVector();

  v0rest4.boost(CMboostV);
  v1rest4.boost(CMboostV);
  Hep3Vector v0rest3 = v0rest4.vect();
  Hep3Vector v1rest3 = v1rest4.vect();
  if( v0rest3.mag() < v1rest3.mag() ) return 1;
  else if( v0rest3.mag() == v1rest3.mag() ) return 0;
  else return -1;
}

//   decending order.
#if defined(__GNUG__)
int SortByECL( const Mdst_ecl **a,  
	       const Mdst_ecl **b )
{
#else
extern "C" int SortByECL( const void *av,  
			  const void *bv )
{
  const Mdst_ecl **a((const Mdst_ecl**)av);
  const Mdst_ecl **b((const Mdst_ecl**)bv);
#endif
  float v0 = (*a)->energy();
  float v1 = (*b)->energy();
  if( v0 < v1 ) return 1;
  else if( v0 == v1 ) return 0;
  else return -1;
}

// Constructor
EventClassify::EventClassify ( void ) 
{

  //*---
  // good track selection criteria.
  //   all events except for cosmic rays are selected
  //   using good tracks defined here.
  //*---
  //  R_cut = 30.0;
  //  Z_cut_min = -40.0;
  //  Z_cut_max =  40.0;
  // updated on 19-March-1999. I.Adachi
  //  cos_trk_min = -0.85;
  //  cos_trk_max =  0.95;
  //
  //  R_cut = 1.0;
  //  Z_cut_min = -5.0;
  //  Z_cut_max =  5.0;
  //  cos_trk_min =  -1.0;
  //  cos_trk_max =   1.0;
  //  P_t_cut = 0.10;

  R_cut = 2.0;
  //*---- 
  //  -3cm < Z < 3cm ---> -4cm < Z < 4cm.
  //
  // good track selection on Z widened... since IP shifted in range of
  // -0.8 - 1.0 cm in Z.
  //                          I.Adachi   '99-Aug-24.
  //*----
  Z_cut_min = -4.0;
  Z_cut_max =  4.0;
  //  Z_cut_min = -3.0;
  //  Z_cut_max =  3.0;
  //  Z_cut_min = -5.0;
  //  Z_cut_max =  5.0;

  // Z-cut for beam gas monitoring.
  Z_cut_min_bg = -10.0;
  Z_cut_max_bg =  10.0;

  // Z-cut for tau selection.
  Z_cut_min_tau = -5.0;
  Z_cut_max_tau =  5.0;

  cos_trk_min =  -1.0;
  cos_trk_max =   1.0;
  P_t_cut = 0.10;
  Nsvd_cut  = 0;

  // default( in case of good_track_type == 1 ).
  conf_cut  = 1.0e-25;

  // cut for # of associated wire hits(instead of CL cut).
  //   effective in good_track_type == 2.
  N_wirehits_cut = 16;

  //*---
  // track selections for two photon events.
  //*---
  drmax1 = 5.0;
  drmax2 = 2.0;
  drmax3 = 1.0;
  dzmax  = 5.0;
  ptmin1 = 0.1;
  ptmin2 = 0.3;
  costhmin = -0.8660;
  costhmax =  0.9563;

  //*---
  // good cluster selection criteria.
  // E > 0.1GeV/c && cluster_quality should be good.
  //*---
  E_min_cut = 0.1; 

  //*---
  // select good acceptance region for ECL. add on '99-Sep-28.
  //*---
  ECLtheta_min =  17.0;   // min. angle(deg).
  ECLtheta_max = 150.0;   // max. angle(deg).

  //*---
  // select acceptance for ECL cluster counting. add on '00-Mar-10.
  //*---
  ECLcos_min = -0.7;   // min. cosine.
  ECLcos_max =  0.9;   // max. cosine.

  //*---
  // hadronic event selection criteria.
  //*---
  Ntrk_cut = 3;
  Ev_cut = 0.4;   // changed '99-May-14.
  Pz_cut = 1.0;
  Esum_min_cut = 0.05;
  Esum_max_cut = 1.8 ;

  //  PVr_cut = 1.0;   //primary vertex cut '99-Jul-01.
  //  PVz_cut = 2.0;
  //
  // prmary vertex cut value also widened because of
  // large shift of IP as -0.6 - 1.0 cm.
  //                      I.Adachi   '99-Aug-24 
  PVr_cut = 1.5;
  PVz_cut = 3.5;

  // new cuts for Ntrk=2 events. add on '99-May-14.
  Ev_low_cut = 0.4;  
  Pz_low_cut = 1.0;
  Psum_low_cut = 1.8;  // add on '99-May-17.

  // add for updating HadronB cuts. IA. '00-Feb-08.
  // barrel track definition.
  cos_barrel_min = -0.643;
  cos_barrel_max =  0.866;
  ECLsumB_min  = 0.2 ; // in unit of 0.5*CM-energy.
  ECLsumB_max  = 1.6 ; // in unit of 0.5*CM-energy.
  //*obsolute   Nbarrel_cut = 2   ; // # of tracks in barrel.
  Necl_cut    = 2   ; // # of ECL clusters.

  //*---
  // Jpsi candidates selection criteria. add on '99-Oct-19 IA.
  //*---
  Plepton_cut = 0.8;      // lepton caididate momentum >0.8 GeV/c
  JpsiMass_min_cut = 2.5; // min. value of J/psi mass region.
  JpsiMass_max_cut = 4.0; // max. value of J/psi mass region.

  //*---
  // Bhabha event selection criteria.
  //*---
  // for rest frame.
  N_bhabha_trk_cut =  2;
  N_bhabha_cls_cut =  2;
  P_bhabha_cut     =  0.5;
  E_bhabha_1st_cut   =  2.0; //  2.0GeV for 1st E cut.
  E_bhabha_total_cut =  4.0; //  4.0GeV for total energy cut.
  acol_angle_cut   = 10.0;

  // regions are defined in negatively charged track.
  cos_barrel_rest1   =  46.7419; // barrel region cut.
  cos_barrel_rest2   = 145.715;
  // restricted Barrel region for Luminosity.
   cos_barrel_rest11   =  50.5; // barrel region cut.
   cos_barrel_rest12   =  129.5;
  //cos_barrel_rest11   =  46.7419; // barrel region cut.
  //cos_barrel_rest12   = 145.715;
  cos_foreward_rest1 =  23.0462; // foreward region cut.
  cos_foreward_rest2 =  45.262;
  cos_backward_rest1 = 146.804;  // backward region cut.
  cos_backward_rest2 = 160.594;

  // for lab frame
  P_bhabha_lab_cut1 =  4.0;  // Pe- for 4.0GeV/c
  P_bhabha_lab_cut2 =  3.0;  // Pe+ for 3.0GeV/c
  acop_angle_cut    = 10.0;

  //*---
  // gamma-pair event selection criteria.
  //*---
  // for rest frame.
  N_gamma_trk_cut =  2;  // max. # of good tracks allowed. modified at "001011.
  N_gamma_cls_cut =  2;  // min. # of good clusters required.  
  E_gamma_cut     =  0.5; // energy deposit cut in unit of W/2. 
  acol_gamma_cut   = 15.0;  // cluster acollinearity cut in deg. modified at "001011.
  E_gamma_1st_cut  =  2.0;   // 2.0GeV for 1st E cut.
  E_gamma_total_cut =  4.0;  //  4.0GeV for total energy cut.

  //*---
  // mu-pair event selection criteria.
  //*---
  // rest frame.
  N_mu_trk_cut = 2;
  P_mu_cut     = 0.5; //  track momentum in units of W/2.
  //  E_mu_cut     = 1.0; //  GeV unit. Maximun energy deposit.
  T_mu_cut     = 8.0; //  nsec. Time difference allowed.
  E_mu_sum_cut = 2.0; //  total deposit energy cut in GeV.
  E_mu_trk_sum_cut =  2.0; //GM 01/03/05; 1.0; //GM 12/11/2004, total track associated ECL energysum  
  acol_mu_cut  = 10.0; //  acollinearily angle cut in deg.

  //*---
  // calorimeter Bhabha selections.
  //*---
  N_CalBhabha_cls_cut = 2;   // min. # of clusters.
  E_CalBhabha_tot_cut = 5.0; // total energy in ECL(GeV).
  E_CalBhabha_1st_cut = 1.0; // max. cluster energy cut(GeV).
  E_CalBhabha_2nd_cut = 1.0; // 2nd max. cluster energy cut(GeV).
  E_CalBhabha_3rd_cut = 0.5; // 3rd max. cluster energy cut(GeV). '99-Oct-19.
  E_CalBhabha_sum12_cut = 14.0; // sum of two energitic cluster energy(GeV).

  //*---
  // tau-pair selection criteria.
  //*---
  Nmin_tau_cut = 2; // min. # of good tracks.  
  Nmax_tau_cut = 8; // max. # of good tracks.  

  Psum_tau_cut = 10; // max. allowed momentum sum(GeV).
  Esum_tau_cut = 10; // max. allowed ECL energy sum(GeV).

  // add on '99-Nov-12. IA. 
  Ptmax_tau_cut         = 0.5;
  Erec_tau_2Dcut        = 3.0;
  //  Ptmax_tau_2Dcut       = 0.8;
  Ptmax_tau_2Dcut       = 1.0;
  Etot_tau_2_4track_cut = 9.0;
  ChargeSum_tau_cut     = 2;

  Vr_tau_cut = 0.5;  // add Primary Vtx cut. '99-Jul-01
  Vz_tau_cut = 3.0;

  Psum_tau_2track_cut = 9.0;  // add further cuts for 2 track case. '99-Aug-23
  Esum_tau_2track_cut = 9.0;
  ThetaMiss_tau_min_cut =   5.0;
  ThetaMiss_tau_max_cut = 175.0;

  //*---
  // cosmic event selection criteria.
  //*---
  //  R_SVD_cut         = 6.05; // SVD outer radius(cm).
  //  R_cos_min_cut     = 2.0;  // beampipe outer radius(cm).
  R_cos_min_cut     = 3.0;  // beampipe outer radius(cm).
  //  R_cos_max_cut     = 8.5;  // CDC inner radius(cm).
  Rimpact_cut       = 0.15; // 3sigma from cosmic study.
  Zimpact_cut       = 3.0;  // 3sigma from cosmic study.
  //  Pz_cosmic_cut     = 0.2;  // this is temporaly. +/-20% allowed.
  Pz_cosmic_cut     = 0.05;  // modified. 99-Mar-29. IA.
  tof_time_diff_cut = 7.0;  // unit of nsec.
  Esum_MIP_cut      = 2.0;  // unit of GeV.

  //*---
  // radiative Bhabha events selection criteria.  '99-May-26. IA
  //*---
  Ntrk_eeg_cut =  1;     // # of min. tracks.
  Ncls_eeg_cut =  3;     // # of clusters in ECL.
  E3sum_eeg_cut = 8.0;   // sum of the 3 bigest cluster energy.
  E_eeg_min_cut = 0.1;   // min. energy deposit. 
  E_eeg_total_cut = 11.0;   // total ECL energy. add on '99-Dec-13.
  Q_eeg_cut = 0.0;          // charge sum of tracks. add on '99-Dec-13.

  //*---
  // eeee event selection criteria.
  //*---
  Vz_4e_cut_level     =1;  // switch for vertex cut method
  //
  // modified on '99-Aug-24. 
  //  Vr_4e_cut           =0.2; // Vertex R in cm
  //  Vz_4e_cut           =2.5; // Vertex Z in cm
  //  Vz_4e_cut2          =0.6; // 2nd level cut of Vertex Z
  // 
  Vr_4e_cut           =0.3; // Vertex R in cm
  Vz_4e_cut           =3.2; // Vertex Z in cm
  Vz_4e_cut2          =0.8; // 2nd level cut of Vertex Z
  Pt_4e_cut           =0.25; // Pt cut
  S_Pz_4e_cut         =-1.; // Pz balance in CMS flame
  S_Pt_4e_cut         =0.2; // Pt balance in CMS flame
  W_4e_cut            =5.0; // Invariant mass of 2-trk
  E_tot_4e_max        =6.0; // Max. total energy deposit at ECL
  E_tot_4e_min        =0.6; // Min. total energy deposit at ECL

  //*---
  // gamma+phi event selection criteria.
  //*---
  N_gphi_cls_cut  = 1;
  N_gphi_trk_cut  = 2;
  E_gphi_min_cut  = 4.5;
  E_gphi_max_cut  = 5.5;
  eprob_gphi_tight_cut  = 0.1;
  kid_gphi_cut    = 0.6;

  //*---
  // radiative mu-pair selecion criteria. add on '99-Oct-18 IA.
  //*---
  N_radmu_trk_cut = 1;
  P_radmu_cm_cut  = 5.0;  // 4.8
  tof_radmu_time_diff_cut = 5.0;  // time difference (5nsec)
  // ---- gamma case.
  N_radmu_ecl_case1       = 4;
  ECL_radmu_min_case1     = 0.8;  // 1.5
  ECL_radmu_max_case1     = 8.0; // 7.5;  // 6.0;
  Gamma_radmu_min_case1   = 0.5;  // 1.0;
  Gamma_radmu_max_case1   = 7.5; // 6.5;  // 5.5;
  E_radmu_total_cut_case1 = 9.5;  // 9.0;
  // ---- no gamma case.
  N_radmu_ecl_case2       = 2;
  ECL_radmu_max_case2     = 1.0;
  acol_radmu_cut          = 20.0;
  acop_radmu_cut          = 170.0;
  P_radmu_lab_cut         = 0.5;

  //*---
  // Bhabha selections initial values. by Adachi 2005-April-11 from Kakuno san's notice.
  //*---
  pres_bh_other_trg = 0;
  pres_gg_other_trg = 0;

  //*---
  // beam BG event selection criteria.
  //*---
  Esum_BG_cut = 1.0;

  //*---
  // particle masses.
  //*---
  pi_mass = 0.13956995;
  mu_mass = 0.105658389;
  e_mass  = 0.51099907e-3;

  //*---
  // machine parameters.
  //*---
  //  E_LER =3.5;
  //  E_HER =7.9965;
  //  double theta=0.022;
  cm=HepLorentzVector(-EHER*sin(theta), 0., -EHER*cos(theta)+ELER,
		      EHER+ELER );
  ebeam=cm.mag()/2.;
  CMBoost=cm.boostVector();

}

// init function
void EventClassify::init ( int *status ) 
{
  *status = 0; 
  //*---
  // print out selection cuts.
  //*---

  dout(Debugout::INFO,"EventClassify") << " " <<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "   << good track criteria used in selections >>          "<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "             Rclosest < "<< R_cut        <<std::endl;
  dout(Debugout::INFO,"EventClassify") << "        "<< Z_cut_min << " < Zclosest < " << Z_cut_max   <<std::endl;
  dout(Debugout::INFO,"EventClassify") << "                           Pt > "<< P_t_cut << " GeV/c "   <<std::endl;
  dout(Debugout::INFO,"EventClassify") << "             confidence level < "<< conf_cut               <<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "   << good cluster criteria used in selections >>          "<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "                            E > "<< E_min_cut<<"GeV"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "                   cluster quality should be good        "<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "   << hadronic event selections criteria >>              "<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "              # of good tracks >= " << Ntrk_cut-1    <<std::endl;
  dout(Debugout::INFO,"EventClassify") << "  visible energy(CMS)/(0.5Ecm) >= " << Ev_cut        <<std::endl;
  dout(Debugout::INFO,"EventClassify") << "               abs(Pzsum(CMS)) <= " << Pz_cut        <<std::endl;
  dout(Debugout::INFO,"EventClassify") << "          Esum in ECL/(0.5Ecm) >= " << Esum_min_cut  <<std::endl;
  dout(Debugout::INFO,"EventClassify") << "          Esum in ECL/(0.5Ecm) <= " << Esum_max_cut  <<std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << "                 jet mass(CMS) >= " << Mjet_cut<<"GeV"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "   << Bhabha event selections criteria >>              "<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "              # of good tracks >= " << N_bhabha_trk_cut   <<std::endl;
  dout(Debugout::INFO,"EventClassify") << "            # of good clusters >= " << N_bhabha_cls_cut   <<std::endl;
  dout(Debugout::INFO,"EventClassify") << "                track momentum  > " << P_bhabha_cut<<"*W/2"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "           max. cluster energy  > " << E_bhabha_1st_cut<<"*W/2"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "                 ECL energy sum > " << E_bhabha_total_cut<<"*W/2"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "            acollinearity angle < " << acol_angle_cut<<"deg"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "   << gamma-pair event selections criteria >>            "<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "              # of good tracks <= " << N_gamma_trk_cut   <<std::endl;
  dout(Debugout::INFO,"EventClassify") << "            # of good clusters >= " << N_gamma_cls_cut   <<std::endl;
  dout(Debugout::INFO,"EventClassify") << "           max. cluster energy  > " << E_gamma_1st_cut<<"*W/2"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "                 ECL energy sum > " << E_gamma_total_cut<<"*W/2"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "           acollinearity angle < " << acol_gamma_cut<<"deg"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "   << mu-pair event selections criteria >>              "<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "      # of high momentum tracks = " << N_mu_trk_cut   <<std::endl;
  dout(Debugout::INFO,"EventClassify") << "         -high momentum means P > " << P_mu_cut<<"*W/2"<<std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << "           Ecluster associated  < " << E_mu_cut<<"GeV"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "     ECL total energy deposite  < " << E_mu_sum_cut<<"GeV"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " Trk associated ECL total energy deposite  < " << E_mu_trk_sum_cut<<"GeV"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "   TOF difference(if available) < " << T_mu_cut<<"nsec"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "            acollinearity angle < " << acol_mu_cut<<"deg"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "   << eeee event selections criteria >>              "<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << "             Vertex_z cut Level = " << Vz_4e_cut_level<< std::endl;
  dout(Debugout::INFO,"EventClassify") << "             Vertex_r           < " << Vr_4e_cut<<"cm"<< std::endl;
  dout(Debugout::INFO,"EventClassify") << "             Vertex_z           < " << Vz_4e_cut<<"cm"<< std::endl;
  dout(Debugout::INFO,"EventClassify") << "             Vertex_z_tight     < " << Vz_4e_cut2<<"cm"<< std::endl;
  dout(Debugout::INFO,"EventClassify") << "             Pt                 > " << Pt_4e_cut<<"GeV/c"<< std::endl;
  dout(Debugout::INFO,"EventClassify") << "             Sigma Pz*          < " << S_Pz_4e_cut<<"GeV/c"<< std::endl;
  dout(Debugout::INFO,"EventClassify") << "             Sigma Pt*          < " << S_Pt_4e_cut<<"GeV/c"<< std::endl;
  dout(Debugout::INFO,"EventClassify") << "             invariant mass (W) < " << W_4e_cut<<"GeV/c^2"<< std::endl;
  dout(Debugout::INFO,"EventClassify") << "         max total E deposit    < " << E_tot_4e_max<<"GeV"<< std::endl;
  dout(Debugout::INFO,"EventClassify") << "         min total E deposit    > " << E_tot_4e_min<<"GeV"<< std::endl;
  dout(Debugout::INFO,"EventClassify") << " --------------------------------------------------------"<<std::endl;
  dout(Debugout::INFO,"EventClassify") << " " <<std::endl;

}

// initialize event data.
void EventClassify::init_event ( int* evtcls_flag, int* hadronic_flag )
{

  //*--
  // initialize evtcls flags.
  //*---
  for( int i=0;i<20;i++){
    evtcls_flag[i]=0;
    hadronic_flag[i]=0;
  }

  //*---
  // delete all panther tables containing already in the event
  //   except for evtcls_flag.
  // Now evtcls_flag also reset at the beginning. '99-Oct-07 I.Adachi
  //*---
  // evtcls_flag.
  Evtcls_flag_Manager& evtflag = Evtcls_flag_Manager::get_manager();
  evtflag.remove();
  // evtcls_flag2.
  Evtcls_flag2_Manager& evtflag2 = Evtcls_flag2_Manager::get_manager();
  evtflag2.remove();

  // hadronic_flag.
  Evtcls_hadronic_flag_Manager& hadflag = Evtcls_hadronic_flag_Manager::get_manager();
  hadflag.remove();

  // hadronic events.
  Evtcls_hadron_info_Manager& hadinfo = Evtcls_hadron_info_Manager::get_manager();
  hadinfo.remove();
  Evtcls_hadron_charged_Manager& hadchg = Evtcls_hadron_charged_Manager::get_manager();
  hadchg.remove();
  Evtcls_hadron_neutral_Manager& hadcls = Evtcls_hadron_neutral_Manager::get_manager();
  hadcls.remove();

  // Bhabha events.
  Evtcls_bhabha_info_Manager& bhainfo = Evtcls_bhabha_info_Manager::get_manager();
  bhainfo.remove();
  Evtcls_bhabha_charged_Manager& bhachg = Evtcls_bhabha_charged_Manager::get_manager();
  bhachg.remove();
  Evtcls_bhabha_cluster_Manager& bhacls = Evtcls_bhabha_cluster_Manager::get_manager();
  bhacls.remove();

  // Gamma_pair events.
  Evtcls_gamma_info_Manager& gaminfo = Evtcls_gamma_info_Manager::get_manager();
  gaminfo.remove();
  Evtcls_gamma_cluster_Manager& gamcls = Evtcls_gamma_cluster_Manager::get_manager();
  gamcls.remove();

  // mu_pair events.
  Evtcls_mupair_info_Manager& muinfo = Evtcls_mupair_info_Manager::get_manager();
  muinfo.remove();
  Evtcls_mupair_charged_Manager& muchg = Evtcls_mupair_charged_Manager::get_manager();
  muchg.remove();

  // radiative Bhabha_pair events.
  Evtcls_radbha_info_Manager& radinfo = Evtcls_radbha_info_Manager::get_manager();
  radinfo.remove();
  Evtcls_radbha_cluster_Manager& radcls = Evtcls_radbha_cluster_Manager::get_manager();
  radcls.remove();
  
  //two-photons
  Evtcls_twophoton_info_Manager& twoinfo 
    = Evtcls_twophoton_info_Manager::get_manager();
  twoinfo.remove();
  Evtcls_twophoton_level1_Manager& twolv1 = 
    Evtcls_twophoton_level1_Manager::get_manager();
  twolv1.remove();
  Evtcls_twophoton_level2_Manager& twolv2 =
    Evtcls_twophoton_level2_Manager::get_manager();
  twolv2.remove();
  Evtcls_twophoton_level3_Manager& twolv3 = 
    Evtcls_twophoton_level3_Manager::get_manager();
  twolv3.remove();

  //cosmic
  Evtcls_cosmic_info_Manager& cosinfo 
    = Evtcls_cosmic_info_Manager::get_manager();
  cosinfo.remove();
  Evtcls_cosmic_charged_Manager& coschg 
    = Evtcls_cosmic_charged_Manager::get_manager();
  coschg.remove();

  //*---
  // get manager for belle_event & check data_type.
  //*---
  Belle_event_Manager& RunHeadMgr = 
    Belle_event_Manager::get_manager();
  std::vector<Belle_event>::iterator it = RunHeadMgr.begin();
  if( it != RunHeadMgr.end() ) {
    data_type = (*it).ExpMC();
    //add 100908 
    expno = (*it).ExpNo();
    runno = (*it).RunNo();

    // dout(Debugout::DDEBUG,"EventClassify") << " data_type : " << data_type <<std::endl;
  }else{
    // dout(Debugout::DDEBUG,"EventClassify") << "%-evtcls undefined data_type. assume MC data..."<<std::endl;
    data_type = 2;
  }

  //*---
  // set initial values for trigger info.
  //*---
  // add on "011001
   bhabha_trg     = -99;
   bhabha_pre_trg = -99;
   bhabha_cdcbb_trg = -99;
   bhabha_klm_trg   = -99;
   bhabha_lum_trg   = -99;
   bhabha_bar_trg=-99;
   gg_bar_trg=-99;
   gg_cal_trg=-99;
   bhabha_flat_trg=-99;
   bhabha_cal_trg=-99;
   bhabha_other_trg=-99; 
   bhabha_dimu_trg=-99; 
 //*---
  // check Bhabha trigger info. if real data.
  //*---
  if( data_type >= 2 ) return;


    Rectrg_summary3_Manager& RecTrgMgr =Rectrg_summary3_Manager::get_manager();
  // change from end()-1 -> begin().


    //add at "021002  begin

    TrgBit TrgBits;
    std::vector<Rectrg_summary3>::iterator it1 = RecTrgMgr.begin();

    if( it1 != RecTrgMgr.end() && *it1 ){

    bhabha_pre_trg = TrgBits.getPSNM(TrgBit::csi_bhabha_p);
    bhabha_trg = TrgBits.getPSNM(TrgBit::csi_bhabha);
    bhabha_lum_trg = TrgBits.getPSNM(TrgBit::csi_lum_e);

      bhabha_cdcbb_trg =  TrgBits.get(TrgBit::cdc_bb);
      bhabha_dimu_trg = TrgBits.getPSNM(TrgBit::dimu) ;
     
     bhabha_klm_trg = (int)TrgBits.get(TrgBit::klm_fwd)+2*(int)TrgBits.get(TrgBit::klm_bwd)+4*(int)TrgBits.get(TrgBit::klm_brl);
     
    bhabha_flat_trg= TrgBits.getPSNM(TrgBit::csi_bhabha_p);

    //add 100908 
    if(expno>13 ){
      bhabha_bar_trg=TrgBits.getPSNM(TrgBit::brl_bhabha);
      gg_bar_trg=  TrgBits.getPSNM(TrgBit::brl_2gamma);
    }   
    else{ 
      bhabha_bar_trg= bhabha_trg;
      gg_bar_trg= 0; 
    }

    //add at "021002  end 
    bhabha_cal_trg= (int) (  bhabha_trg==1 || bhabha_lum_trg==1 || bhabha_bar_trg==1 || (bhabha_dimu_trg==1 && bhabha_klm_trg==0) );     
    if( bhabha_cal_trg != 1 ) bhabha_other_trg=1; 

    gg_cal_trg= (int) ( bhabha_trg==1 || bhabha_lum_trg==1 || bhabha_bar_trg==1 ||  gg_bar_trg==1 );     
    if(gg_cal_trg != 1) gg_other_trg=1;   
    } 

    if(bhabha_other_trg==1) pres_bh_other_trg++;
    if(gg_other_trg==1) pres_gg_other_trg++;
  
  // dout(Debugout::DDEBUG,"EventClassify") << " Bhabha trg : " << bhabha_trg <<std::endl;

}

// fill primary_vertex info.
void EventClassify::fill_vertex ( void )
{

  //*---
  // Get managers for Primary_Vertex.
  //*---
  Evtvtx_primary_vertex_Manager& VtxMgr 
    = Evtvtx_primary_vertex_Manager::get_manager();
  std::vector<Evtvtx_primary_vertex>::iterator it = VtxMgr.begin();

  if( it != VtxMgr.end() && *it ){
    VtxQuality = (*it).quality();
    if( VtxQuality >= 2 ){
      PrimaryVtx = HepPoint3D( (*it).PV_x(), (*it).PV_y(), (*it).PV_z() );
    }
  }else{
    VtxQuality = -99;
  }

}

// fill MDST_Charged info.
void EventClassify::fill_charged ( int good_track_type )
{

  //*---
  // Get managers for MDST_Charged and MDST_Trk.
  //*---
  Mdst_charged_Manager& ChgMgr = Mdst_charged_Manager::get_manager();
  Mdst_trk_Manager& TrkMgr = Mdst_trk_Manager::get_manager();
  Mdst_trk_fit_Manager& TrkfitMgr = Mdst_trk_fit_Manager::get_manager();

  //*---
  // Pick up good tracks from charged data.
  //*---
  // Loop over charged tracks
  for(std::vector<Mdst_charged>::iterator it = ChgMgr.begin();
      it != ChgMgr.end(); ++it){
    Mdst_charged& Charged = *it;           // Get one track
    Mdst_trk& Trk = Charged.trk();         // Obtain pointer to MDST_trk
    Mdst_trk_fit& Trkfit  = Trk.mhyp(2);   // get pointer to MDST_TRK_fit.

    Hep3Vector P ( Charged.px(), Charged.py(), Charged.pz() ); 
    float cosine_trk = P.z()/P.mag();
    float dr       = Trkfit.helix(0);
    float dz       = Trkfit.helix(3);
     int   wirehits = Trkfit.nhits(0);
    float conf_level;
    if( Trkfit.ndf() ==0 ){
      conf_level = -999.0;
    }else{
      conf_level = prob_(Trkfit.chisq(),Trkfit.ndf());
    }      

    // track selection here.
    // P_t cut
    if ( P.perp() >= P_t_cut ) {
      if( cosine_trk > cos_trk_min && cosine_trk < cos_trk_max ){
	if( Trkfit.nhits(3) >= Nsvd_cut ){
	  if( abs(dr)< R_cut ){
	    if( ( dz> Z_cut_min ) && (dz< Z_cut_max) ){
	      if( good_track_type == 1 ){
		if( conf_level > conf_cut ){

		  // Register good track from CL cut.
		  TrkList.append ( Charged ); 
		}
	      }else if( good_track_type == 2 ){
		if( wirehits >= N_wirehits_cut ){

		  // Register good track from # of wire hit cut.
		  TrkList.append ( Charged ); 
		}
	      }else{
		
		// Register good track w/o further cuts.
		TrkList.append ( Charged ); 
	      }
	    }

	    // track candidate for tau selection.
	    if( (dz> Z_cut_min_tau) && (dz< Z_cut_max_tau) ){
	      TrkListTau.append( Charged );
	    }

	    // track candidate for monitoring beam gas events.
	    if( (dz> Z_cut_min_bg) && (dz< Z_cut_max_bg) ){
	      TrkListBG.append( Charged );
	    }

	  }
	  //    	  }else if( abs(dr) > R_cos_min_cut ) {
    	  if( abs(dr) < R_cos_min_cut ) {
	    // Register a track coming from off-vertex as cosmic candidate.
	    //	      TrkList_offvtx.push_back ( &Charged );
	    TrkList_offvtx.append( Charged );
	  }
	}
      }
    }
    // this is track selection for two photon events.
    //  from Uehara_san.
    if( Trk.quality() == 0 ){
      TrkList_TwoPho.append( Charged );
      
      // classify tracks into three levels.
      if( P.perp() > ptmin1 && abs(dr) <= drmax1 && abs(dz) <=dzmax ) 
    	TrkList_TwoPhoLevel1.append( Charged );  // TL=1
      if( P.perp() > ptmin2 && abs(dr) <= drmax2 && abs(dz) <= dzmax 
    	  &&   // for V0
	  P.cosTheta() >= costhmin && P.cosTheta() <= costhmax ) 
    	TrkList_TwoPhoLevel2.append( Charged );
      if( P.perp() > ptmin2 && abs(dr) <= drmax3 && abs(dz) <= dzmax &&
	  P.cosTheta() >= costhmin && P.cosTheta() <= costhmax )
    	TrkList_TwoPhoLevel3.append( Charged );
    }
    
    // this is track selection for eeee events.
    if( abs(dr)<= Vr_4e_cut ){
      if( abs(dz) <= Vz_4e_cut ){
	TrkList_eeee.append ( Charged ); 
      }                 
    }

    // this is track candidates for gamma+phi events. '99-Aug-20.
    TrkList_gphi.append( Charged );

  }
  

  //*---
  // separate track list by looking at charge.
  //*---
  for( int itrk=0; itrk< TrkList.length(); itrk++ ){
    if( TrkList[itrk]->charge() > 0 ){
      TrkList_pos_rest.append( TrkList[itrk] );
      TrkList_pos.append( TrkList[itrk] );
    }else{
      TrkList_neg_rest.append( TrkList[itrk] );
      TrkList_neg.append( TrkList[itrk] );
    }
  }

  //*---
  // sort track list by track's magnitude of momentum.
  //*---
  TrkList_neg_rest.sort( SortByPrest );
  TrkList_pos_rest.sort( SortByPrest );
  TrkList_neg.sort( SortByP );
  TrkList_pos.sort( SortByP );
  
  //*---
  // fill two highest momentum tracks
  //     as QED event candidate/rest frame.
  //*---
  if( TrkList_neg_rest.length() && TrkList_pos_rest.length() ) {
    TrkList_pair_rest.append( TrkList_neg_rest[0] );
    TrkList_pair_rest.append( TrkList_pos_rest[0] );
  }else if( ! TrkList_neg_rest.length() && TrkList_pos_rest.length() ){
    TrkList_pair_rest.append( TrkList_pos_rest[0] );
    TrkList_pair_rest.append( TrkList_pos_rest[1] );
  }else if( TrkList_neg_rest.length() && ! TrkList_pos_rest.length() ){
    TrkList_pair_rest.append( TrkList_neg_rest[0] );
    TrkList_pair_rest.append( TrkList_neg_rest[1] );
  }

  //*---
  // fill two highest momentum tracks
  //     as QED event candidate/lab frame.
  //*---
  if( TrkList_neg.length() && TrkList_pos.length() ) {
    TrkList_pair.append( TrkList_neg[0] );
    TrkList_pair.append( TrkList_pos[0] );
  }else if( ! TrkList_neg.length() && TrkList_pos.length() ){
    TrkList_pair.append( TrkList_pos[0] );
    TrkList_pair.append( TrkList_pos[1] );
  }else if( TrkList_neg.length() && ! TrkList_pos.length() ){
    TrkList_pair.append( TrkList_neg[0] );
    TrkList_pair.append( TrkList_neg[1] );
  }

  
  //*---
  // that's it.
  //*---
  
}


// fill MDST_Gamma info.
void EventClassify::fill_gamma ( void )
{
  //*---
  // Get Gamma Managers
  //*---
  Mdst_gamma_Manager& GamMgr = Mdst_gamma_Manager::get_manager();
  Mdst_ecl_Manager& EclMgr = Mdst_ecl_Manager::get_manager();

  //*---
  // Select good clusters.
  //*---
  // Loop over gamma tracks
  for(std::vector<Mdst_gamma>::iterator it = GamMgr.begin();
      it != GamMgr.end(); it++){
    Mdst_gamma& Gamma = *it;           // Get one track
    Mdst_ecl& Ecl = Gamma.ecl();         // Obtain pointer to MDST_ecl
    Hep3Vector P ( Gamma.px(), Gamma.py(), Gamma.pz() ); 
    // Define momentum vector
    if ( Ecl.quality() == 0 ) {       // identified as a good cluster
      if ( P.mag() >= E_min_cut ) {         // E cut
	//	GamList.push_back ( &Gamma );
 	GamList.append ( Gamma ); 
 	GamList_rest.append ( Gamma ); 
      }
    }

    // pick up all gamma clusters for phi_gamma selection.
    GamList_gphi.append ( Gamma ); 

  }

  //*---
  // Sort by energy at rest/lab frame.
  //*---
  GamList.sort(SortByE);
  GamList_rest.sort(SortByErest);
  GamList_gphi.sort(SortByErest);

  //*---
  // that's it.
  //*---  

}

// fill MDST_ECL info.
void EventClassify::fill_ecl ( void )
{

  //*---
  // Get ECL Manager.
  //*---
  Mdst_ecl_Manager& EclMgr = Mdst_ecl_Manager::get_manager();

  //*---
  // Check cluster info.
  //*---
  Nlarge      = 0;
  Necl_barrel = 0;
  for(std::vector<Mdst_ecl>::iterator it = EclMgr.begin();
      it != EclMgr.end(); it++){
    Mdst_ecl& cls = *it;

    //append all cluster info.
    EclListAll.append( cls );

    //boosted into the rest frame.
    double Ex = cls.energy()*sin( cls.theta() )*cos( cls.phi() );
    double Ey = cls.energy()*sin( cls.theta() )*sin( cls.phi() );
    double Ez = cls.energy()*cos( cls.theta() );
    HepLorentzVector P4(Ex, Ey, Ez, cls.energy());
    HepLorentzVector Prest=P4;
    // rest-frame. register AList.
    Prest.boost(CMBoost);

    //cluster position.
    float cls_theta = cls.theta()*180.0/acos(-1.0);

    //check cluster quality.
    if( cls.quality() == 0 ){
      //save in the list.
      EclList.append ( cls ); 
      // lab-frame.
      EclListLab.append( new HepLorentzVector(P4) );
      if( Prest.e() >= 0.5 * ebeam ) Nlarge++;
      EclListCM.append( new HepLorentzVector(Prest) );
      //here save "track associated clusters for Mupair selection criteria
      if (cls.match() !=0) EclListTrack.append(new HepLorentzVector(P4) );

      // here save "good" clusters for hadronic selections.
      if( cls.energy() > E_min_cut ){
	if( cls_theta >= ECLtheta_min && cls_theta <= ECLtheta_max )
	  EclListCMSel.append( new HepLorentzVector(Prest) );
	if( cos( cls.theta() ) > ECLcos_min 
	    && cos( cls.theta() ) < ECLcos_max ) Necl_barrel++;
      }
      
    }

    //save all cluster info.
    EclListCMAll.append( new HepLorentzVector(Prest) );

  }

  // dout(Debugout::DDEBUG,"EventClassify") << " before sort " <<std::endl;
  //  for( int i=0; i<EclListCM.length();i++){
  // dout(Debugout::DDEBUG,"EventClassify") << " cluster # : "<<i<<" energy : "<<EclListCM[i]->e()<<std::endl;
  //  } 

  //*---
  // Sort HepAList by energy.
  //*---
  EclListTrack.sort( SortByEnergy );
  EclListLab.sort( SortByEnergy );
  EclListCM.sort( SortByEnergy );
  EclListCMSel.sort( SortByEnergy );
  EclListCMAll.sort( SortByEnergy );

  EclList.sort( SortByECL );
  EclListAll.sort( SortByECL );

  // dout(Debugout::DDEBUG,"EventClassify") << " after sort " <<std::endl;
  //  for( int i=0; i<EclListCM.length();i++){
  // dout(Debugout::DDEBUG,"EventClassify") << " cluster # : "<<i<<" energy : "<<EclListCM[i]->e()<<std::endl;
  //  } 


  //*---
  // that's it.
  //*---  

}

// hadronic event selection.
void EventClassify::hadsel ( int* evtcls_flag, int* hadronic_flag, int output_level )
{

  //*---
  // reset global variables here.
  //*---
  //  for charged tracks.
  float Psum     = 0.0;
  float Esum_chg = 0.0;
  float Psum_cms     = 0.0;
  float Esum_chg_cms = 0.0;
  float Pzsum    = 0.0;
  //  for neutral clusters.
  float Esum_cls     = 0.0;
  float Esum_cls_cms = 0.0;
  //  for ECL clusters.
  float ECL_sum = 0.0;


  //*---
  // calculate charged observables.
  //*---
  int Ntrk = TrkList.length();
  int Nbarrel = 0;

  // loop over # of good tracks.  
  for( int itrk=0; itrk< Ntrk; itrk++ ){

    Hep3Vector P (TrkList[itrk]->px(),TrkList[itrk]->py(),
		  TrkList[itrk]->pz());
    double E =sqrt(P.mag2()+ pow(pi_mass,2.));
    HepLorentzVector P_cms(P,E);
    P_cms.boost(CMBoost);
    Psum += P.mag();
    Esum_chg += E;
    float ppcms = sqrt( pow(P_cms.e(),2.0) - P_cms.mag2() );
    Psum_cms += sqrt( pow(P_cms.e(),2.0) - P_cms.mag2() );
    Esum_chg_cms += P_cms.e();
    Pzsum += P_cms.pz();

    //check if track contains in the barrel region.
    if( P.cosTheta() > cos_barrel_min && P.cosTheta() < 0.866 ) Nbarrel++;

    //store boosted good tracks for thrust calculation.
    //    Hep3Vector P3( P_cms.px(), P_cms.py(), P_cms.pz() );
    //    vlist.append(new Hep3Vector(P3));
    Hep3Vector P3( P_cms.px(), P_cms.py(), P_cms.pz() );
    vlist.push_back( P3 );
    HepLorentzVector P4 = P_cms;
    v4list.append(new HepLorentzVector(P4));

  }

  //*---
  // here look at J/psi candidates. add on '99-Oct-19 IA.
  //*---
  int JpsiCand = 0;
  for( int itrk=0; itrk<Ntrk; itrk++ ){
    if( TrkList[itrk]->charge() > 0 ){
      Hep3Vector P3plus (TrkList[itrk]->px(),TrkList[itrk]->py(),
			 TrkList[itrk]->pz());
      // momentum cut...
      if( P3plus.mag() > Plepton_cut ){

	HepLorentzVector P4plus( P3plus, P3plus.mag() );
      
	for( int jtrk=0; jtrk<Ntrk; jtrk++ ){
	  if( TrkList[jtrk]->charge() < 0 ){
	    Hep3Vector P3minus (TrkList[jtrk]->px(),TrkList[jtrk]->py(),
				TrkList[jtrk]->pz());
	    // momentum cut...
	    if( P3minus.mag() > Plepton_cut ){

	      HepLorentzVector P4minus( P3minus, P3minus.mag() );

	      HepLorentzVector Pjpsi = P4plus + P4minus;
	      if( Pjpsi.mag() > JpsiMass_min_cut &&
		  Pjpsi.mag() < JpsiMass_max_cut )   JpsiCand++;

	    }
	  }
	}
      } 
    }
  }

  //*---
  // calculate neutral observables.
  //*---
  int Ncls = GamList.length();

  // dout(Debugout::DDEBUG,"EventClassify") << " #good clusters: " << Ncls <<std::endl;

  // loop over # of good clusters.
  for( int icls =0; icls<Ncls; icls++ ){

    Hep3Vector P (GamList[icls]->px(),GamList[icls]->py(),
		  GamList[icls]->pz());
    double E = P.mag();
    Esum_cls += E;
    HepLorentzVector P_cms(P,E);
    P_cms.boost(CMBoost);
    Esum_cls_cms += P_cms.e();
    Pzsum += P_cms.pz();
 
    //store boosted good clusters for thrust calculation.
    //    Hep3Vector P3( P_cms.px(), P_cms.py(), P_cms.pz() );
    Hep3Vector P3( P_cms.px(), P_cms.py(), P_cms.pz() );

    vlist.push_back( P3 );
    HepLorentzVector P4 = P_cms;
    v4list.append(new HepLorentzVector(P4));

  }
  
  //*---
  // compute visible energy. 
  //*---
  float Evis     = Esum_cls + Esum_chg;
  float Evis_cms = Esum_cls_cms + Esum_chg_cms;

  //*---
  // ECL deposit energy. add on '99-May-14. I.Adachi
  //                     update to use "good" cluster. '99-Sep-28.
  //*---
  // # of good ECL clusters.
  int Necl = EclListCMSel.length();

  for( int i=0;i<EclListCMSel.length();i++){
    Hep3Vector Ecls = EclListCMSel[i]->vect();
    ECL_sum += Ecls.mag();
  }


  //*---
  //for debug.
  //*---
  // dout(Debugout::DDEBUG,"EventClassify") << "                            Ntrk : " << Ntrk <<std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << "                            Ncls : " << Ncls <<std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << " size of v4list defined HepAList : " <<v4list.length()<<std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << " size of vlist defined list      : " <<vlist.size()<<std::endl;
  
  //*---
  // get thrust for invariant jet mass calculation.
  //*---
  double Hjetmass = 0.0;
  double Ljetmass = 0.0;
  double ThrParam = 0.0;
  double FWParam  = 0.0;

  if( (Ntrk+Ncls) >= 2 ){
    Hep3Vector thr = thrust( vlist.begin(),vlist.end(),SelfFunc(Hep3Vector()));
    Hep3Vector tn  = thr.unit();
    ThrParam = 1.0-thr.mag();

    //jet mass in each hemisphere.
    HepLorentzVector Vpos ;
    HepLorentzVector Vneg ;
    // dout(Debugout::DDEBUG,"EventClassify") << " size of AList : " << v4list.length()<< std::endl;
    // dout(Debugout::DDEBUG,"EventClassify") << " track+cluster : " << Ntrk+Ncls<< std::endl;
    //    for( int ip=0; ip < Ntrk; ip++ ){
    for( int ip=0; ip < Ntrk+Ncls; ip++ ){
      Hep3Vector Vsign(v4list[ip]->px(),v4list[ip]->py(),v4list[ip]->pz());
      double scalar = Vsign*tn;
      if( scalar >= 0.0 ) {
	Vpos += *v4list[ip];
      }else{
	Vneg += *v4list[ip];
      }
    }
    if( Vpos.mag() >= Vneg.mag() ){
      Hjetmass = Vpos.mag();
      Ljetmass = Vneg.mag();
    }else{
      Hjetmass = Vneg.mag();
      Ljetmass = Vpos.mag();
    }
    
    //*---
    // Fox-Wolfram variable.
    //*---
    FoxWolfram fw =
      foxwolfram( vlist.begin(),vlist.end(),SelfFunc(Hep3Vector()));
    FWParam = fw.R(2);
  }
  
  //*---
  // hadronic event classification.
  //*---
  // remove heavy jet mass cut... '99-May-7  I.Adachi
  // add R2 collinear cut.(->now removed)
  // add low multiplicity cut.    '99-May-14 I.Adachi
  // add primary vertex cut if quality >=2
  //                              '99-Jul-01 I.Adachi
  if( Ntrk >= Ntrk_cut ){
    if( (Evis_cms/ebeam) >= Ev_cut ){
      if( abs(Pzsum/ebeam) <= Pz_cut ) {
	if( (ECL_sum/ebeam) >= Esum_min_cut &&
	    (ECL_sum/ebeam) <= Esum_max_cut ) {

	  // if quality >= 2 case.
	  if( VtxQuality >= 2 ){
	    float PVr = sqrt( PrimaryVtx.x()*PrimaryVtx.x()+
			      PrimaryVtx.y()*PrimaryVtx.y() );
	    float PVz = PrimaryVtx.z();
	    if( abs(PVr) <= PVr_cut && abs(PVz)<= PVz_cut ){
	  
	      if( FWParam < 0.2 ){
		evtcls_flag[0]   = 10;
		hadronic_flag[0] = 10;
	      }else{
		evtcls_flag[0]   = 20;
		hadronic_flag[0] = 20;
	      }
	  
	      //*---
	      // more tight cut.
	      //*---
	      if( Ntrk >= 5 ){
		if( (Evis_cms/ebeam) >= 1.0 ){
		  if( abs(Pzsum/ebeam) <= 0.6 ) {

		    if( FWParam < 0.2 ){
		      evtcls_flag[0]   = 30;
		      hadronic_flag[5] = 10;
		    }else{
		      evtcls_flag[0]   = 40;
		      hadronic_flag[5] = 20;
		    }

		  }
		}
	      }

	    }

	  // if quality < 2 case. not PV cut applied.
	  }else{

	    if( FWParam < 0.2 ){
	      evtcls_flag[0]   = 10;
	      hadronic_flag[0] = 10;
	    }else{
	      evtcls_flag[0]   = 20;
	      hadronic_flag[0] = 20;
	    }

	    //*---
	    // more tight cut.
	    //*---
	    if( Ntrk >= 5 ){
	      if( (Evis_cms/ebeam) >= 1.0 ){
		if( abs(Pzsum/ebeam) <= 0.6 ) {
		  
		  if( FWParam < 0.2 ){
		    evtcls_flag[0]   = 30;
		    hadronic_flag[5] = 10;
		  }else{
		    evtcls_flag[0]   = 40;
		    hadronic_flag[5] = 20;
		  }
		  
		}
	      }
	    }
	    
	  }
	  
	}
      }
    }

  }else if( Ntrk == 2 ){
    if( (Evis_cms/ebeam) >= Ev_low_cut ){
      if( abs(Pzsum/ebeam) <= Pz_low_cut ) {
	if( (ECL_sum/ebeam) >= Esum_min_cut &&
	    (ECL_sum/ebeam) <= Esum_max_cut ) {
	  if( (Psum_cms/ebeam) <= Psum_low_cut ){

	    if( FWParam < 0.2 ){
	      evtcls_flag[0] = 1;
	    }else{
	      evtcls_flag[0] = 2;
	    }
	  
	  }
	}
      }
    }
  }

  //*---
  // new HadronA selection. '00-Mar-10. IA. 
  //*---
  if( evtcls_flag[0] >=10 ){
    if( (ECL_sum/ebeam) > ECLsumB_min &&
	(ECL_sum/ebeam) < ECLsumB_max ) {
      if( Necl_barrel >= Necl_cut ){
	if( FWParam < 0.2 ){
	  hadronic_flag[1] = 10;
	  evtcls_flag[13]  = 10;
	}else{
	  hadronic_flag[1] = 20;
	  evtcls_flag[13]  = 20;
	}
      }
    }
  }
  
  //*---
  // HadronB from Brendan Casey's selection. '00-Mar-12 IA.
  //*---
  if( hadronic_flag[1] ){
    Mdst_charged_Manager& ChgMgr = Mdst_charged_Manager::get_manager();
    Mdst_ecl_Manager& EclMgr = Mdst_ecl_Manager::get_manager();

    float EcmDefault = 2.0*ebeam;
    HadB HadronB = good_hadronB(ChgMgr, EclMgr, EcmDefault);
    int cut = HadronB.type;

    if( cut == 0 ){
      if( FWParam < 0.2 ){
	hadronic_flag[2] = 10;
	if( Ntrk>= 5 ) 	hadronic_flag[4] = 10;
      }else{
	hadronic_flag[2] = 20;
	if( Ntrk>= 5 ) 	hadronic_flag[4] = 20;
      }
    }

  }

  //*---
  // # tracks >= 5 added on new HadronA. '00-Mar-10. IA.
  //*---
  if( hadronic_flag[1] ) {
    if( Ntrk >= 5 ) {
      if( FWParam < 0.2 ){
	hadronic_flag[3] = 10;
      }else{
	hadronic_flag[3] = 20;
      }
    }
  }

  //*---
  // flag for J/psi candidates. add on '99-Oct-19 IA.
  //*---
  if( evtcls_flag[0] >= 10 ){
    if( JpsiCand ) evtcls_flag[12] = JpsiCand;
  }
  
  //*---
  // Fill panther table if this event is identified as a hadronic event.
  //*---

  //output for debugging.
  // dout(Debugout::DDEBUG,"EventClassify") << " evtcls_flag(hadsel) : " << evtcls_flag[0]  <<std::endl;

  //++++
  //(1) fill event_flag.
  //++++
  //  Evtcls_flag_Manager& EvtFlgMgr = Evtcls_flag_Manager::get_manager();
  //  Evtcls_flag& evtflag = EvtFlgMgr.add();

  //  evtflag.flag(0,evtcls_hadron_flag);


  //*----
  // here return if this event is not identified as a hadron.
  //*----
  if( output_level == 0 ){ 
    if( evtcls_flag[0] == 0 ) return;
  }

  //++++  
  //(2) fill hadronic event global informations.
  //++++  
  Evtcls_hadron_info_Manager & EvtInfoMgr 
    = Evtcls_hadron_info_Manager::get_manager();
  Evtcls_hadron_info& evtinfo = EvtInfoMgr.add();
  
  evtinfo.Ntrk(Ntrk);
  evtinfo.Ncls(Ncls);
  evtinfo.Psum(Psum_cms);
  evtinfo.Esum(ECL_sum);
  //  evtinfo.Esum(Esum_cls_cms);
  //  evtinfo.ECLsum(ECL_sum);
  evtinfo.Evis(Evis_cms);
  evtinfo.Pz(Pzsum);
  evtinfo.HeavyJetMass(Hjetmass);
  evtinfo.Thrust(1.0-ThrParam);
  evtinfo.R2(FWParam);

  // dout(Debugout::DDEBUG,"EventClassify") << " Ebeam : "<< ebeam <<std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << " Psum : " << Psum_cms << " Esum : " << Esum_cls_cms 
  //       << " Evis : " << Evis_cms << " Pz   : " << Pzsum <<std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << " HJM  : " << Hjetmass << " T    : " << 1.0-ThrParam
  //       << " R2   : " << FWParam  <<std::endl;


  //++++  
  //(3) fill pointers for good tracks to MDST_Charged.
  //++++  
  Evtcls_hadron_charged_Manager & EvtChgMgr 
    = Evtcls_hadron_charged_Manager::get_manager();
  if( Ntrk ){
    for( int i = 0; i<Ntrk; ++i ){
      Evtcls_hadron_charged& evtchg = EvtChgMgr.add();
      evtchg.charged( *TrkList[i] );
    }
  }else{
    Evtcls_hadron_charged& evtchg = EvtChgMgr.add();
    evtchg.reset_charged();
  }

  //++++  
  //(4) fill pointers for good neutrals to MDST_Gamma.
  //++++  
  Evtcls_hadron_neutral_Manager & EvtNeuMgr 
    = Evtcls_hadron_neutral_Manager::get_manager();
  if( Ncls ){
    for( int i = 0; i<Ncls; ++i ){
      Evtcls_hadron_neutral& evtgam = EvtNeuMgr.add();
      evtgam.gamma( *GamList[i] );
    }
  }else{
    Evtcls_hadron_neutral& evtgam = EvtNeuMgr.add();
    evtgam.reset_gamma();
  }

  // dout(Debugout::DDEBUG,"EventClassify") << " panther output for evtcls completed " << std::endl;

}


// Bhabha selection.
void EventClassify::bhasel ( int* evtcls_flag, int output_level )
{

  //*---
  // reset local flags.
  //*---
  int bhabha_rest_flag = 0;
  int bhabha_lab_flag  = 0;

  //*---
  // const.
  //*---
  float const pi = acos(-1.0);

  //*---
  // return if # good tracks < N_bhabha_trk_cut(=2).
  //*---
  int Ntrk = TrkList.length();
  if( Ntrk < N_bhabha_trk_cut ) return;


  // dout(Debugout::DDEBUG,"EventClassify") << " # tracks : "<< Ntrk <<std::endl;

  //*---
  // get momentum for candidate pair.
  //*---
  Hep3Vector P3;
  P3=Hep3Vector( TrkList_pair_rest[0]->px(), TrkList_pair_rest[0]->py(),
		 TrkList_pair_rest[0]->pz() );
  HepLorentzVector PeleCM( P3, sqrt(P3.mag2()+
	       TrkList_pair_rest[0]->mass()*TrkList_pair_rest[0]->mass()) );
  PeleCM.boost(CMBoost);
  Hep3Vector P3ele = PeleCM.vect();
  P3=Hep3Vector( TrkList_pair_rest[1]->px(), TrkList_pair_rest[1]->py(),
		 TrkList_pair_rest[1]->pz() );
  HepLorentzVector PposCM( P3, sqrt(P3.mag2()+
	       TrkList_pair_rest[1]->mass()*TrkList_pair_rest[1]->mass()) );
  PposCM.boost(CMBoost);
  Hep3Vector P3pos = PposCM.vect();

  Hep3Vector P3eleLab( TrkList_pair[0]->px(), TrkList_pair[0]->py(),
		       TrkList_pair[0]->pz() );
  Hep3Vector P3posLab( TrkList_pair[1]->px(), TrkList_pair[1]->py(),
		       TrkList_pair[1]->pz() );

  //*---
  // sum of total energy observed in ECL.
  //*---
  float  ECLsum = 0.0;
  for( int i=0;i<EclListCMAll.length();i++ ){
    ECLsum += EclListCMAll[i] -> e();
  }
  
  //*---
  // calculate acollinearity angle.
  //*---
  float acoll = P3ele.angle( -P3pos )*180.0/pi;

  // add at "001011
  float dthe = P3ele.theta()+P3pos.theta()-pi;
         dthe=asin(sin(dthe))*180.0/pi;
  float dphi = P3ele.phi()-P3pos.phi();
         dphi=asin(sin(dphi))*180.0/pi;
  
  //*---
  // calculate acoplanarity angle.
  //*---
  Hep3Vector P3eleLabxy( P3eleLab.x(), P3eleLab.y(), 0.0 );
  Hep3Vector P3posLabxy( P3posLab.x(), P3posLab.y(), 0.0 );
  float acopl = P3eleLabxy.angle( -P3posLabxy )*180.0/pi;
  
  //*---
  //* event selections here.
  //*---
  // using variables in the rest frame.
  if(  ( (P3ele.mag()/ebeam) > P_bhabha_cut ) && 
       ( (P3pos.mag()/ebeam) > P_bhabha_cut ) ){


    //add on "001012
    if( (bhabha_dimu_trg  && !bhabha_klm_trg   ) ||  data_type != 1  ){

      if( acoll < acol_angle_cut ) {
        bhabha_rest_flag =40;
      } else {
        bhabha_rest_flag =41;
      }
    } 

   //add on 20031021
  //*---
  // return if # good clusters < N_bhabha_cls_cut(=2).
  //*---
  if( EclListCMAll.length() < N_bhabha_cls_cut ) return;
  // dout(Debugout::DDEBUG,"EventClassify") << " # clusters : "<< EclList.length() <<std::endl;




    // new cuts from June-3
    //    if( ( (EclListCM[0]->e()/ebeam) > -99.0 ) &&
    //	( (EclListCM[1]->e()/ebeam) > -99.0 ) ){

    if( ( EclListCMAll[0]->e() > E_bhabha_1st_cut ) &&
	( ECLsum > E_bhabha_total_cut ) ){


      // '99-Jun-01
      //    if( ( (EclListCM[0]->e()/ebeam) > P_bhabha_cut ) &&
      //	( (EclListCM[1]->e()/ebeam) > P_bhabha_cut ) ){

      // now events pass through cuts.
      // assign rough flag here.

      // this flag will be over-written if accepted by angle cuts.

      //add 011001 
      if(data_type != 1 || bhabha_cal_trg==1 || (bhabha_other_trg==1 &&  pres_bh_other_trg%10==1) ) {  	  
      if( acoll < acol_angle_cut ){
	bhabha_rest_flag =40;
      } else {
	bhabha_rest_flag =41;
      }
      }

      //      dout(Debugout::INFO,"EventClassify") << "My_check:: "<<  bhabha_cal_trg<<" "<< bhabha_flat_trg<<" "<< bhabha_other_trg<<":"<< pres_bh_other_trg<<" "<< pres_cal_trg<<" "<< pres_flat_trg<< " "<< pres_dimu_trg<< " "<< bhabha_rest_flag <<" "<<std::endl ;

      if(bhabha_rest_flag == 40) {
      float theta_ele =  P3ele.theta()*180.0/pi;
      float theta_pos =  P3pos.theta()*180.0/pi;

      // this is berrel region.
      // both e-/e+ tracks should be in the barrel region.
      if( theta_ele >=  cos_barrel_rest11 &&
	  theta_ele <=  cos_barrel_rest12 ){
	if( theta_pos >= cos_barrel_rest11 &&
	    theta_pos <= cos_barrel_rest12 ){

	  // collinear Bhabha events.
	  if( acoll < acol_angle_cut ) {
	    bhabha_rest_flag =10;
	    // loose Bhabha.
	  } else {
	    bhabha_rest_flag =11;
	  }

	}

      // this is forward region.
      } else if( theta_ele >=  cos_foreward_rest1 &&
		theta_ele <=  cos_foreward_rest2 ){
	if( theta_pos >= cos_backward_rest1 &&
	    theta_pos <= cos_backward_rest2 ){

	  // collinear Bhabha events.
	  if( acoll < acol_angle_cut ) {
	    bhabha_rest_flag =20;
	    // loose Bhabha.
	  } else {
	    bhabha_rest_flag =21;
	  }
	  
	}
	
      // this is forward region.
      }else if( theta_ele >=  cos_backward_rest1 &&
		theta_ele <=  cos_backward_rest2 ){
	if( theta_pos >=  cos_foreward_rest1 &&
	    theta_pos <=  cos_foreward_rest2 ){

	  // collinear Bhabha events.
	  if( acoll < acol_angle_cut ) {
	    bhabha_rest_flag =30;
	    // loose Bhabha.
	  } else {
	    bhabha_rest_flag =31;
	  }
	  
	}

      // other regions.
      }else{
	// collinear Bhabha events.
	if( acoll < acol_angle_cut ) {
	  bhabha_rest_flag =40;
	  // loose Bhabha.
	} else {
	  bhabha_rest_flag =41;
	}
      }
      }
    }
  }


  
  //*---
  // check if Bhabha trigger is fired.
  //*---
  if( data_type == 1 ){
 //add at "021002  
    if( bhabha_bar_trg == 1 ){
      if( bhabha_rest_flag ){
	bhabha_rest_flag = bhabha_rest_flag + 100;
      }
    }
  }

  evtcls_flag[1] = bhabha_rest_flag;
  
  //*---
  // fill panther table.
  //*---
  if( output_level == 0 ){
    if( evtcls_flag[1] == 0 ) return;
  }  

  //*---
  // global informations.
  //*---
  Evtcls_bhabha_info_Manager & BhaMgr 
    = Evtcls_bhabha_info_Manager::get_manager();
  Evtcls_bhabha_info& evtinfo = BhaMgr.add();

  evtinfo.Ntrk(Ntrk);
  evtinfo.Ncls(EclListCMAll.length());
  evtinfo.E1st(EclListCMAll[0]->e());
  evtinfo.E2nd(EclListCMAll[1]->e());
  if( P3pos.mag() > P3ele.mag() ){
    evtinfo.P1st(P3pos.mag());
    evtinfo.P2nd(P3ele.mag());
  }else{
    evtinfo.P1st(P3ele.mag());
    evtinfo.P2nd(P3pos.mag());
  }
  evtinfo.acollinearity(acoll);
  
  //*---  
  // fill pointers for good tracks to MDST_Charged.
  //*---  
  Evtcls_bhabha_charged_Manager & EvtChgMgr 
    = Evtcls_bhabha_charged_Manager::get_manager();

  if( Ntrk ){
    for( int i=0; i< 2 ; i++ ){
      Evtcls_bhabha_charged& evtchg = EvtChgMgr.add();
      evtchg.charged( *TrkList_pair_rest[i] );
    }
  }else{
    Evtcls_bhabha_charged& evtchg = EvtChgMgr.add();
    evtchg.reset_charged();
  }
  
  //*---  
  // fill pointers for good clusters to MDST_Ecl.
  //*---  
  Evtcls_bhabha_cluster_Manager & EvtClsMgr 
    = Evtcls_bhabha_cluster_Manager::get_manager();

  if( EclList.length() ){
    for( int i=0; i< EclList.length(); i++ ){
      Evtcls_bhabha_cluster& evtcls = EvtClsMgr.add();
      evtcls.ecl( *EclList[i] );
    }
  }else{
    Evtcls_bhabha_cluster& evtcls = EvtClsMgr.add();
    evtcls.reset_ecl();
  }

}

// gamma-pair selection.
void EventClassify::ggpair_sel ( int* evtcls_flag, int output_level )
{

  //*---
  // const.
  //*---
  float const pi = acos(-1.0);

  //*---
  // reset local variables.
  //*---
  int gamma_rest_flag = 0;
  
  // dout(Debugout::DDEBUG,"EventClassify") << " #tracks   : "<<TrkList.length()<<std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << " #clusters : "<<GamList.length()<<std::endl;

  //*---
  //* check # of good tracks.
  //*---

  //Allow to pass ECL selected Bhabha
   if( TrkList.length() > N_gamma_trk_cut+3 ) return;

  //*---
  // check # of clusters.
  //*---
  if( GamList.length() < N_gamma_cls_cut ) return;

  //*---
  // get the 1st & 2nd largest energy clusters.
  //*--- 
  Hep3Vector Erest[2]; 
  Hep3Vector Elab[2]; 
  for( int i=0;i<2;i++){
    Hep3Vector P ( GamList_rest[i]->px(), GamList_rest[i]->py(),
		   GamList_rest[i]->pz()  );
    HepLorentzVector EposCM( P, P.mag() );
    EposCM.boost(CMBoost);
    Erest[i]=EposCM.vect();
    Elab[i] =Hep3Vector( GamList[i]->px(), GamList[i]->py(),
			 GamList[i]->pz() );

  }

  //*---
  // sum of total energy observed in ECL.
  //*---
  float  ECLsum = 0.0;
  for( int i=0;i<EclListCM.length();i++ ){
    ECLsum += EclListCM[i] -> e();
  }
  
  //*---
  //* selections.
  //*---
  float acol = Erest[0].angle( -Erest[1] )*180.0/pi;
  float acop =  Elab[0].angle(  -Elab[1] )*180.0/pi;

  // add on "001012
  float dthe =Erest[0].theta()+Erest[1].theta()-pi;
  dthe=asin(sin(dthe))*180.0/pi;
  float dphi = Erest[0].phi()-Erest[1].phi();
  dphi=asin(sin(dphi))*180.0/pi;

  if( ( Erest[0].mag() > E_gamma_1st_cut ) &&
      ( ECLsum > E_gamma_total_cut ) ){

    // add on "001012
 
   if( TrkList.length() == 2 ) {   

      Hep3Vector P3eleLab( TrkList_pair[0]->px(), TrkList_pair[0]->py(),
			   TrkList_pair[0]->pz() );
      Hep3Vector P3posLab( TrkList_pair[1]->px(), TrkList_pair[1]->py(),
			   TrkList_pair[1]->pz() );
      HepLorentzVector PeleCM( P3eleLab, P3eleLab.mag() );
      HepLorentzVector PposCM( P3posLab, P3posLab.mag() );  
      PeleCM.boost(CMBoost);
      PposCM.boost(CMBoost);

      float dthc=P3eleLab.angle(P3posLab); 
      float minv2=2*P3eleLab.mag()*P3posLab.mag()*(1-cos(dthc));
      dthc = dthc*180.0/pi;




	// Select converted GG->ee 

      if( dthc < 10 &&  minv2 < 0.04 ) {
	gamma_rest_flag =40;
	//	dout(Debugout::INFO,"EventClassify") << "My_checkGG::2_1 "<< acol <<" "<< dphi<<" "<<dthc<<std::endl;
      }  

//add 011001 
    
      // Select Bhabha with bad tracking     
      else if( abs(acol) < 5 &&  abs(dphi) > 5 &&  abs(dphi) < 15 && ((PeleCM.mag() < P_bhabha_cut  || PposCM.mag() < P_bhabha_cut) || ((PeleCM.mag() > P_bhabha_cut  || PposCM.mag() > P_bhabha_cut) &&  abs(dthc)>5)) ) { 
	gamma_rest_flag =40;
	//	dout(Debugout::INFO,"EventClassify") << "My_checkGG::2_2 "<< acol <<" "<<dphi <<" "<<PeleCM.mag()<<" "<<dthc<<std::endl;
      }
   }

   //   }	  

 
  //add 011001 Select Bhabha with bad tracking 
   if(data_type != 1 || gg_cal_trg==1 || 
      (gg_other_trg==1 &&  pres_gg_other_trg%10==1) ){
    
    if( TrkList.length() < 2) {
      if( acol < acol_gamma_cut ){
	gamma_rest_flag = 40;
      }else{
	gamma_rest_flag = 41;
      }
    }
   
      float cos_gamma0 = Erest[0].theta()*180.0/pi;
      float cos_gamma1 = Erest[1].theta()*180.0/pi;
    

      //add 100908 Selected Bhabha  events  with ECL only
      if(data_type != 1 || bhabha_bar_trg==1 || bhabha_lum_trg==1 ){
 	if( cos_gamma0 >= cos_barrel_rest11-1 &&
	    cos_gamma0 <= cos_barrel_rest12+1 ){
	  if( cos_gamma1 >= cos_barrel_rest11-1 &&
	      cos_gamma1 <= cos_barrel_rest12+1 ){
	    
	    if( acol < acol_gamma_cut-5 && abs(dphi) > 2. && abs(dphi) < 15. ){
	      gamma_rest_flag = 50;
	    }
	  }
	}
      }
   

      if( TrkList.length() < 2) {
   
      if( cos_gamma0 >= cos_barrel_rest11 &&
	  cos_gamma0 <= cos_barrel_rest12 ){
	if( cos_gamma1 >= cos_barrel_rest11 &&
	    cos_gamma1 <= cos_barrel_rest12 ){
	
	  if( acol < acol_gamma_cut-5 && abs(dphi) < 2.3 ){
	    gamma_rest_flag = 10;
	  }
	  else{ 
	    gamma_rest_flag = 11;
	  }
	}
      } else if( cos_gamma0 >= 29 + 0*cos_foreward_rest1 &&
		 cos_gamma0 <= cos_foreward_rest2 ){
	if( cos_gamma1 >= cos_backward_rest1 &&
	    cos_gamma1 <= 151 + 0*cos_backward_rest2 ){

	  if( acol < acol_gamma_cut && abs(dphi)< 2.3 ){
	    gamma_rest_flag = 20;
	  }
	  else{
	    gamma_rest_flag = 21;
	  }
	}

      } else if( cos_gamma0 >= cos_backward_rest1 &&
		 cos_gamma0 <= 151 + 0*cos_backward_rest2 ){
	if( cos_gamma1 >= 29 + 0*cos_foreward_rest1 &&
	    cos_gamma1 <= cos_foreward_rest2 ){
	  if( acol < acol_gamma_cut && abs(dphi)< 2.3 ){
	    gamma_rest_flag = 30;

	  }else{
	    gamma_rest_flag = 31;
	  }
	}
      }else{
	if( acol < acol_gamma_cut ){
	  gamma_rest_flag = 40;
	}else{
	  gamma_rest_flag = 41;
	}
      }
   
      //add 011001 - select Bhabha events with bad tracking

      if( gamma_rest_flag != 10 && gamma_rest_flag  != 20 &&  gamma_rest_flag  != 30  && gamma_rest_flag != 40 ) { 

	if(  acol < 5  && abs(dphi) > 2.0 && abs(dphi) < 15 ){

	  if( cos_gamma0 >= 18 &&  cos_gamma0 <= 152 && cos_gamma1 >= 18 &&  cos_gamma1 <= 152 ){ 
	    //	    dout(Debugout::INFO,"EventClassify") << "My_checkGG::1_1 "<< acol <<" "<< dphi<<" "<< bhabha_trg<<" "<< bhabha_lum_trg<< " fl "<< gamma_rest_flag <<  std::endl;
	    gamma_rest_flag = 40;
	  }

	}
      }

    }
    
   }
  }

  // dout(Debugout::DDEBUG,"EventClassify") << " gamma_rest_flag : " << gamma_rest_flag << std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << "  E1 : "<< Erest[0].mag()<< " E2 : "<<Erest[1].mag()<<std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << "  acolinearity : " << acol << std::endl;
  
  //*---
  // check if Bhabha trigger is fired.
  //*---

  if( data_type == 1 ){
    //add at "021002
    if( bhabha_bar_trg == 1 ){
      if( gamma_rest_flag ){
	gamma_rest_flag = gamma_rest_flag + 100;
      }
    }
  }

  //*---
  // put panther data.
  //*---
  evtcls_flag[2] = gamma_rest_flag;
  if( ! output_level ){
    if( ! evtcls_flag[2] ) return;
  }

  // dout(Debugout::DDEBUG,"EventClassify") << " evtcls_flag[2] : "<< evtcls_flag[2]<<std::endl;

  //*---
  // global informations.
  //*---
  Evtcls_gamma_info_Manager & GamMgr 
    = Evtcls_gamma_info_Manager::get_manager();
  Evtcls_gamma_info& evtinfo = GamMgr.add();

  evtinfo.Ntrk(TrkList.length());
  evtinfo.Ncls(GamList.length());
  evtinfo.E1st(Erest[0].mag());
  evtinfo.E2nd(Erest[1].mag());
  evtinfo.acollinearity(acol);

  // dout(Debugout::DDEBUG,"EventClassify") << " add info completed.. " <<std::endl;
  
  //*---  
  // fill pointers for good clusters to MDST_gamma.
  //*---  
  Evtcls_gamma_cluster_Manager & EvtClsMgr 
    = Evtcls_gamma_cluster_Manager::get_manager();

  if( GamList_rest.length() ){
    for( int i=0; i< 2 ; i++ ){
      Evtcls_gamma_cluster& evtcls = EvtClsMgr.add();
      evtcls.gamma( *GamList_rest[i] );
    }
  }else{
    Evtcls_gamma_cluster& evtcls = EvtClsMgr.add();
    evtcls.reset_gamma();
  }

  //*---
  // that's it.
  //*---

  }

// mu-pair selection.
void EventClassify::musel ( int* evtcls_flag, int trkecl_crit, int output_level )
{

  //*---
  // const.
  //*---
  float const pi = acos(-1.0);

  //*---
  // reset local variables.
  //*---
  int Nmu_cand = 0;
  int mu_rest_flag = 0;

  //*---
  // # good tracks.
  //*---
  int Ntrk = TrkList.length();

  //*---
  // check # of high momentum tracks.
  //*---
  for( int itrk=0;itrk<Ntrk; itrk++){
    Hep3Vector P(TrkList[itrk]->px(),TrkList[itrk]->py(),
		 TrkList[itrk]->pz());
    HepLorentzVector P4(P,sqrt(P.mag2()+
		    TrkList[itrk]->mass()*TrkList[itrk]->mass()));
    P4.boost(CMBoost);
    Hep3Vector Pcm = P4.vect();
    if( (Pcm.mag()/ebeam) > P_mu_cut ) Nmu_cand++;
  }
  // return if Nmu_cand != 2.

  if( Nmu_cand != N_mu_trk_cut ) return;

  //*---
  // get managers for local use.
  //*--
  Mdst_tof_Manager& TofMgr = Mdst_tof_Manager::get_manager();
  Mdst_ecl_Manager& EclMgr = Mdst_ecl_Manager::get_manager();

  //* obsolute!
  //*  Mdst_klm_Manager& KlmMgr = Mdst_klm_Manager::get_manager();
  //* obsolute! 
  //*  Mdst_sim_xref_Manager& SimMgr = Mdst_sim_xref_Manager::get_manager();

  //*---
  // get momentum & ToF info for candidate pair.
  //*---
  Hep3Vector P3;
  P3=Hep3Vector( TrkList_pair_rest[0]->px(), TrkList_pair_rest[0]->py(),
		 TrkList_pair_rest[0]->pz() );
  HepLorentzVector PmumCM( P3, sqrt(P3.mag2()+
				    TrkList_pair_rest[0]->mass()*TrkList_pair_rest[0]->mass()) );
  PmumCM.boost(CMBoost);
  Hep3Vector P3mum = PmumCM.vect();
  P3=Hep3Vector( TrkList_pair_rest[1]->px(), TrkList_pair_rest[1]->py(),
		 TrkList_pair_rest[1]->pz() );
  HepLorentzVector PmupCM( P3, sqrt(P3.mag2()+
				    TrkList_pair_rest[1]->mass()*TrkList_pair_rest[1]->mass()) );
  PmupCM.boost(CMBoost);
  Hep3Vector P3mup = PmupCM.vect();
  
  Hep3Vector P3mumLab( TrkList_pair[0]->px(), TrkList_pair[0]->py(),
		       TrkList_pair[0]->pz() );
  Hep3Vector P3mupLab( TrkList_pair[1]->px(), TrkList_pair[1]->py(),
		       TrkList_pair[1]->pz() );

  float  ToF[2];
  float  ToFLab[2];
  for( int i=0; i<2 ;i++ ){
    Mdst_tof& TOF1 = TrkList_pair_rest[i]->tof();
    Mdst_tof& TOF2 = TrkList_pair[i]->tof();
    if( TOF1 ){
      ToF[i] = TOF1.tof();
    }else{
      ToF[i] = -9999.0;
    }
    if( TOF2 ){
      ToFLab[i] = TOF2.tof();
    }else{
      ToFLab[i] = -9999.0;
    }
  }
  
  //*---
  // get cluster information associated with charged tracks.
  //     using sim_Xref info.
  //*---
  //  float EmumRest;
  //  float EmupRest;
  //  float EmumLab;
  //  float EmupLab;
  //  for(vector<Mdst_sim_xref>::iterator it_sim = SimMgr.begin(); 
  //    it_sim != SimMgr.end(); ++it_sim){
  //    Mdst_sim_xref& Sim = *it_sim;
  //    Mdst_charged& SimChg = Sim.charged();
  //    Mdst_ecl& SimEcl = Sim.ecl();

    // dout(Debugout::DDEBUG,"EventClassify") << " it_sim->ecl "<< it_sim->ecl() <<std::endl;

    // for rest frame.
  //    if ((it_sim -> get_ID()) == (TrkList_pair_rest[0] -> get_ID())){
  //      if (it_sim -> ecl()){
  //	if(SimEcl.quality() == 0){
  //            float ecl_px = SimEcl.energy()*sin(SimEcl.theta())*
  //              cos(SimEcl.phi());
  //            float ecl_py = SimEcl.energy()*sin(SimEcl.theta())*
  //              sin(SimEcl.phi());
  //            float ecl_pz = SimEcl.energy()*cos(SimEcl.theta());
  //            HepLorentzVector Erest ( ecl_px, ecl_py, ecl_pz,
  //				     SimEcl.energy()); 
  //	    Erest.boost(CMBoost);
  //	    EmumRest = Erest.e();

	    
	    //	    dout(Debugout::INFO,"EventClassify") << " SimECL energy: " <<SimEcl.energy()<<std::endl; 
	    //	    dout(Debugout::INFO,"EventClassify") << "      EmumRest: " <<EmumRest<<std::endl; 

  //	}
  //      }
  //    }else
  //    if ((it_sim -> get_ID()) == (TrkList_pair_rest[1] -> get_ID())){
  //      if (it_sim -> ecl()){
  //	if(SimEcl.quality() == 0){
  //            float ecl_px = SimEcl.energy()*sin(SimEcl.theta())*
  //              cos(SimEcl.phi());
  //            float ecl_py = SimEcl.energy()*sin(SimEcl.theta())*
  //              sin(SimEcl.phi());
  //            float ecl_pz = SimEcl.energy()*cos(SimEcl.theta());
  //            HepLorentzVector Erest ( ecl_px, ecl_py, ecl_pz,
  //				     SimEcl.energy()); 
  //	    Erest.boost(CMBoost);
  //	    EmupRest = Erest.e();
  //	}
  //      }
  //    }

    // for Lab frame.
  //    if ((it_sim -> get_ID()) == (TrkList_pair[0] -> get_ID())){
  //      if (it_sim -> ecl()){
  //	if(SimEcl.quality() == 0){
  //	  EmumLab = SimEcl.energy();
  //	}
  //      }
  //    }else
  //    if ((it_sim -> get_ID()) == (TrkList_pair[1] -> get_ID())){
  //      if (it_sim -> ecl()){
  //	if(SimEcl.quality() == 0){
  //	  EmupLab = SimEcl.energy();
  //	}
  //      }
  //    }
  //  }

  //*----
  // energy cut was modified. 99-April-27. IA.
  // So far, cluster energy associated with a charged track was
  // examined(Ecls<E_mu_cut(=1.0GeV). Now, we look at total energy
  // sum observed in ECL(Esum<E_mu_sum_cut=2GeV).
  //*---

  //*---
  // # of ECL good clusters.
  //*---
  int Necl = EclListCM.length();

  float Emupair_allsum = 0.0;
  //*---
  // Cluster energy sum.
  //*---  
  // Modified by G. Majumder on 4 Oct 2004
  // Use track matched cluster energy only

  for(int iecl=0;iecl<Necl;iecl++){
    Emupair_allsum += EclListCM[iecl]->e();
  }
  
  float Emupair_sum = 0.0;
  for (int iecl = 0; iecl < EclListTrack.length(); iecl++) {
    Emupair_sum += EclListTrack[iecl]->e();
  }

  //*---
  // return if Esum == 0 (caused by ECL truncation). '00-May-11. IA
  //*---
  if( Emupair_allsum == 0.0 ) return;

  //*----
  // selection here for rest frame.
  //*---
  float acoll = P3mum.angle( -P3mup )*180.0/pi;
  float Time_diff = abs( ToF[0] - ToF[1] );
  int tof_flag = 1;
  if( ToF[0]<-999.0 || ToF[1]<-999.0 ) tof_flag = 0;

  float acopl = P3mumLab.angle( -P3mupLab )*180.0/pi;
  float cos_theta_Lab = P3mumLab.theta()*180.0/pi;

  if( ( (P3mum.mag()/ebeam) > P_mu_cut ) &&
      ( (P3mup.mag()/ebeam) > P_mu_cut ) ){

    //    mu_rest_flag = 1000;

    // now look at track associated cluster energy.

    if ( (trkecl_crit ==0 && Emupair_allsum < E_mu_sum_cut) || 
	 (trkecl_crit ==1 && Emupair_sum < E_mu_trk_sum_cut)){

      // now events pass through cuts.
      // assign rough flag here.
      // this flag will be over-written if accepted by angle cuts.

      //GM on 12/11/2004 : This was not the case for forward-barrel-backward
      //                   Gap region e.g. muminus in barrel region, but
      //                   muplus in FB region, mu_rest_flag was not overwritten

      // It was accepting full polar angle of muminus, but not muplus


      if( acoll < acol_mu_cut ) {
	mu_rest_flag =40;
      } else {
	mu_rest_flag =41;
      }
	
      float theta_mum = P3mum.theta()*180.0/pi;
      float theta_mup = P3mup.theta()*180.0/pi;

      // this is barrel region.
      if( theta_mum >=  cos_barrel_rest11 &&
	  theta_mum <=  cos_barrel_rest12 &&
	  theta_mup >=  cos_barrel_rest11 &&
	  theta_mup <=  cos_barrel_rest12 ){

	  if( acoll < acol_mu_cut ) {
	    if( tof_flag ){
	      if( Time_diff < T_mu_cut ){
		// collinear mu-pair events w/ time diff OK.
		mu_rest_flag =10;
	      }
	    }else{
	      // collinear mu-pair events w/o time data.
	      mu_rest_flag =11;
	    }
	    // loose mu-pair. time data not checked.
	  } else {
	    if( tof_flag ){
	      if( Time_diff < T_mu_cut ){
		// loose mu-pair events w/ time diff OK.
		mu_rest_flag =12;
	      }
	    }else{
	      // loose mu-pair events w/o time data.
	      mu_rest_flag =13;
	    }
	  }

	  //	}
	
	
	// this is foreward region.
      }else if( theta_mum >=  cos_foreward_rest1 &&
		theta_mum <=  cos_foreward_rest2 &&
		theta_mup >=  cos_backward_rest1 &&
		theta_mup <=  cos_backward_rest2 ){
	  // collinear mu-pair events.
	  if( acoll < acol_angle_cut ) {
	    mu_rest_flag =20;
	    // loose mu-pair.
	  } else {
	    mu_rest_flag =21;
	  }
	  //	}

	// this is backward region.
      }else if( theta_mum >=  cos_backward_rest1 &&
		theta_mum <=  cos_backward_rest2 &&
		theta_mup >=  cos_foreward_rest1 &&
		theta_mup <=  cos_foreward_rest2 ){
	  // collinear mu-pair events.
	  if( acoll < acol_angle_cut ) {
	    mu_rest_flag =30;
	    // loose mu-pair.
	  } else {
	    mu_rest_flag =31;
	  }
	  //	}
	  
      }else if ((theta_mum >=  cos_barrel_rest1 &&
		 theta_mum <=  cos_barrel_rest2 &&
		 ((theta_mup >=  cos_foreward_rest1 &&
		   theta_mup <=  cos_foreward_rest2) ||
		  (theta_mup >=  cos_backward_rest1 &&
		   theta_mup <=  cos_backward_rest2))) ||
		(theta_mup >=  cos_barrel_rest1 &&
		 theta_mup <=  cos_barrel_rest2 &&
		 ((theta_mum >=  cos_foreward_rest1 &&
		   theta_mum <=  cos_foreward_rest2) ||
		  (theta_mum >=  cos_backward_rest1 &&
		   theta_mum <=  cos_backward_rest2)))) {
	
	
	if( acoll < acol_angle_cut ) {
	  mu_rest_flag =40;
	  // loose mu-pair.
	} else {
	  mu_rest_flag =41;
	}
	 
      }else {
	// other regions.
	if( acoll < acol_angle_cut ) {
	  mu_rest_flag =40;
	  // loose mu-pair.
	} else {
	  mu_rest_flag =41;
	}
      }
    }
  }

  
  evtcls_flag[3] = mu_rest_flag;
  
  //*---
  // fill panther table.
  //*---
  if( output_level == 0 ){
    if( evtcls_flag[3] == 0 ) return;
  }

  //*---
  // global informations.
  //*---
  Evtcls_mupair_info_Manager & MuMgr 
    = Evtcls_mupair_info_Manager::get_manager();
  Evtcls_mupair_info& muoninfo = MuMgr.add();

  muoninfo.Ntrk(Ntrk);
  muoninfo.Ncls(EclList.length());
  if( P3mup.mag() > P3mup.mag() ){
    muoninfo.P1st(P3mup.mag());
    muoninfo.P2nd(P3mum.mag());
  }else{
    muoninfo.P1st(P3mum.mag());
    muoninfo.P2nd(P3mup.mag());
  }
  if( tof_flag ){
    muoninfo.ToF_diff(Time_diff);
  }else{
    muoninfo.ToF_diff(-9999.0);
  }
  muoninfo.acollinearity(acoll);
  //*---  
  // fill pointers for mu candidate tracks to MDST_Charged.
  //*---  
  Evtcls_mupair_charged_Manager & EvtChgMgr 
    = Evtcls_mupair_charged_Manager::get_manager();

  if( Ntrk ){
    for( int i=0;i<2;i++ ){
      Evtcls_mupair_charged& evtchg = EvtChgMgr.add();
      evtchg.charged( *TrkList_pair_rest[i] );
    }
  }else{
    Evtcls_mupair_charged& evtchg = EvtChgMgr.add();
    evtchg.reset_charged();
  }

  //*---
  // that's it.
  //*---

}

// calorimeter Bhabha event selection.
void EventClassify::CalBhabha_sel ( int* evtcls_flag, int output_level )
{

  //*---
  // reset local variables.
  //*---
  float Esum = 0.0;
  float const pi = acos(-1.0);

  //*---
  // return if no cluster found.
  //*---
  if( EclListAll.length() < N_CalBhabha_cls_cut ) return;

  //*--
  // compute total energy in ECL.
  //*---
  for( int i=0;i<EclListAll.length();i++ ){
    Esum += EclListAll[i]->energy();
  }

  //*---
  // get two of the most energtic clusters.
  //*---
  float E1st = EclListAll[0]->energy();
  float E2nd = EclListAll[1]->energy();

  // dout(Debugout::DDEBUG,"EventClassify") << " Esum : "<< Esum << std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << " E1st : "<<E1st<<"   E2nd : "<<E2nd<<std::endl;

  //*---
  // return if the 3rd max. cluster energy > 0.5GeV. '99-Oct-19. IA.
  //*---
  if( EclListAll.length() >= 3 ){
    float E3rd = EclListAll[2]->energy();  
    if( E3rd >= E_CalBhabha_3rd_cut ) return;
  }

  //*---
  // selection starts here.
  //*---
  if( Esum > E_CalBhabha_tot_cut ){
    if( ( E1st+E2nd ) < E_CalBhabha_sum12_cut ){
      if( E1st > E_CalBhabha_1st_cut && E2nd > E_CalBhabha_2nd_cut ){
	
	evtcls_flag[9] = 10;
	
      }
    }
  }
  

  //*---
  // check if Bhabha trigger is fired.
  //*---
  if( data_type == 1 ){
    if( bhabha_pre_trg == 1 ){
      if( evtcls_flag[9] ){
	evtcls_flag[9] = evtcls_flag[9] + 100;
      }
    }
  }
  
  
  
}



// tau-pair event selection.
void EventClassify::tausel ( int* evtcls_flag, int output_level )
{

  //*---
  // numerical constants.
  //*---
  float const static pi = acos(-1.0);

  //*---
  // reset local variables.
  //*---
  int tau_flag = 0;
  int tau_flag2 = 0;   // add on 2002-Dec-12 from Inami san.
  int tau_flag3 = 0;   // MOD
  float Psum = 0.0;
  float Esum = 0.0;
  float EsumCM = 0.0;
  float Ptmax = 0.0;
  float Egamma_sum = 0.0;
  int   ChargeSum  = 0;
  int Nbarrel = 0;
  float AngleMax = 0.0;

  //*---
  // for missing momentum calculation.
  //*---
  Hep3Vector Pmiss = -(Hep3Vector)cm;
  Hep3Vector Pmiss_cm = Hep3Vector(0,0,0);

  //*---
  // return if # of good tracks is out of range.
  //*---
  int Ntrk = TrkListTau.length();
  if( Ntrk < Nmin_tau_cut || Ntrk > Nmax_tau_cut ) return;

  //*---
  // primary vertex cut.
  //*---
  if( VtxQuality >= 2 ){
      float Vr = sqrt( PrimaryVtx.x()*PrimaryVtx.x() +
		       PrimaryVtx.y()*PrimaryVtx.y() );
      float Vz = PrimaryVtx.z();
      if( abs( Vz ) >= Vz_tau_cut ) return;  // MOD
      if( abs( Vr ) >= 1.0 ) return;  // MOD
      if( abs( Vr ) >= Vr_tau_cut ) tau_flag3=1;  // MOD
  }

  //*---
  // charged track
  //*---
  for( int i=0; i<Ntrk;i++ ){

    Hep3Vector P (TrkListTau[i]->px(),TrkListTau[i]->py(),
		  TrkListTau[i]->pz());
    HepLorentzVector Prest(P, sqrt(P.mag2()+ pow(pi_mass,2.)) );
    Prest.boost(CMBoost);
    Hep3Vector P3cm = Prest.vect();
    Psum += P3cm.mag();

    if( P.perp() > Ptmax ) Ptmax = P.perp();

    float theta = P.theta()*180/pi;
    if(30<theta && theta<130) Nbarrel++;

    //missing momentum.
    Pmiss -= P;
    Pmiss_cm -= P3cm;

    //charge sum.
    ChargeSum += (int) TrkListTau[i]->charge();

    for(int j=i+1;j<Ntrk;j++){
      Hep3Vector P2(TrkListTau[j]->px(),TrkListTau[j]->py(),
		    TrkListTau[j]->pz());
      HepLorentzVector Prest2(P2, sqrt(P2.mag2()+ pow(pi_mass,2.)) );
      Prest2.boost(CMBoost);
      Hep3Vector P3cm2 = Prest2.vect();

      float angle = P3cm2.angle(P3cm);
      if(angle>AngleMax) AngleMax=angle;
    }
  }

  //*---
  // the simple cut for Pt_max.
  //*---
  if( Ptmax <= Ptmax_tau_cut ) return;

  //*---
  // ECL sum.
  //*---
  for( int i=0;i<EclList.length();i++ ){
    Esum += EclList[i]->energy();
  }
  for( int i=0;i<EclListCM.length();i++ ){
    EsumCM += EclListCM[i]->e();
  }

  //*---
  // Egamma sum.
  //*---
  for( int i=0;i<GamList.length();i++ ){
    Hep3Vector P ( GamList[i]->px(), GamList[i]->py(),
                   GamList[i]->pz()  );
    HepLorentzVector Pgamma_cm( P, P.mag() );
    Pgamma_cm.boost(CMBoost);
    Egamma_sum += Pgamma_cm.e();

    //missing momentum.
    Pmiss -= P;
    Pmiss_cm -= Pgamma_cm.vect();
  }

  //*---
  // missing momentum angle.
  //*---
  float ThetaMiss = Pmiss.theta()*180.0/pi;

  //*---
  // Erec, Etot, Etrk
  //*---
  float Erec = Psum + Egamma_sum;
  float Etot = Erec + Pmiss_cm.mag();
  float Etrk = EsumCM - Egamma_sum;

  //*---
  // selection here.
  //*---
  if( Psum >= Psum_tau_cut ) tau_flag3=1;  // MOD
  if( Esum >= Esum_tau_cut ) tau_flag3=1;  // MOD
  // cuts for 2 track case
  if( Ntrk == 2 ){
    if( Psum >= Psum_tau_2track_cut ) tau_flag3=1;  // MOD
    if( Esum >= 11.0 ) return;  // MOD
    if( Esum >= Esum_tau_2track_cut ) tau_flag3=1;  // MOD
    if(ThetaMiss <= ThetaMiss_tau_min_cut ||
       ThetaMiss >= ThetaMiss_tau_max_cut ) return;
  }
  // 2D cut (Erec,Ptmax)
  if( Erec <= Erec_tau_2Dcut && Ptmax <= Ptmax_tau_2Dcut ) return;
  // cut for 2-4 track case
  if (Ntrk >= 2 && Ntrk <= 4 ){

    // cut shown below updated on '00-Sept-11. By Inami_kun.
    //    if( Etot >= Etot_tau_2_4track_cut ) return;

    if( Etot >= Etot_tau_2_4track_cut && AngleMax*180/pi >= 175 &&  // MOD
	( Esum <= 2 || 10 <= Esum ) ) return; // MOD
    if( Etot >= Etot_tau_2_4track_cut && AngleMax*180/pi >= 175 )  // MOD
      tau_flag3=1; // MOD
    if( Nbarrel<2 && Etrk >= 5.3 ) return;
  }
  if( abs(ChargeSum) > ChargeSum_tau_cut ) return;

  //*---
  // skip : modified : 2002/12/09 K.Inami   ... added
  //*---
  if( AngleMax*180/pi <= 20 ){
    tau_flag2=1;
  }

  //*---
  // after cuts
  //*---
  tau_flag = Ntrk*10 + tau_flag2*1000 + tau_flag3*10000;  // MOD

  //*---
  // store flag into evtcls_flag.
  //*---
  evtcls_flag[4] = tau_flag;



}

// two-photon event selection.
void EventClassify::twophoton_sel ( int* evtcls_flag, int output_level )
{

  // Modified by S.Uehara 7-Mar-2000
  // To set a tighter cut to reduce number of events selected
  //  Added selection criteria: (common for all the event classes)
  //     U1-1) There are at least one positive and one negative tracks
  //                             among level-1 tracks
  //     U1-2) Invariant mass constructed by level-3 tracks <= 4.5 GeV
  //     U1-3) missing mass squared calculated by level-3 tracks
  //            >= 2.0 GeV^2
  //

  //*---
  // reset local variables.
  //*---
  float esum = 0.0;
  float psum = 0.0;

  //U2-2) esum above threshold
  float egamat = 0.0;
  float egamthre = 0.04;

  int   ntr0, ntr1, ntr2, ntr3;
  ntr0 = ntr1 = ntr2 = ntr3 = 0;
  Hep3Vector P_sum(0., 0., 0.);
  float C_sum = 0.0;
  float C_sum_L = 0.0;
  
  // local variable for the added criteria
  float ntr1_p = 0;      //U1-1
  float ntr1_n = 0;      //U1-1
  HepLorentzVector P4_sum_tr3(0., 0., 0., 0.);   //U1-2, U1-3
  float tp_wmax = 4.5;
  float tp_mm2min = 2.0;
  
  //*---
  // # of tracks assigned into each level.
  //*---
  ntr0 = TrkList_TwoPho.length();
  ntr1 = TrkList_TwoPhoLevel1.length();
  ntr2 = TrkList_TwoPhoLevel2.length();
  ntr3 = TrkList_TwoPhoLevel3.length();
  
  //*---
  // return here if no candidates tracks.
  //*---
  if( ! ntr0 ) return;
  
  //*---
  // accumulate ECL deposit energy.
  //*---
  for( int i=0;i<EclList.length(); i++ ){
    esum += EclList[i]->energy();
  }

  //U2-2  
  for( int i=0;i<GamList.length(); i++ ){
    float egam = sqrt( pow(GamList[i]->px(),2) +
       pow(GamList[i]->py(),2) + pow(GamList[i]->pz(),2) );
    if( egam > egamthre ) egamat += egam;  
  }


 
  //*---
  // accumulate charged track momentum.
  //*---
  for( int i=0;i<ntr0; i++ ){
    
    Hep3Vector P3( TrkList_TwoPho[i]->px(), TrkList_TwoPho[i]->py(),
                   TrkList_TwoPho[i]->pz() );
    psum += P3.mag();
  }
  
  //*---
  // level-1 track.
  //*---
  for( int i=0;i<ntr1;i++ ){
    
    Hep3Vector P3( TrkList_TwoPhoLevel1[i]->px(),
                   TrkList_TwoPhoLevel1[i]->py(),
                   TrkList_TwoPhoLevel1[i]->pz() );
    HepLorentzVector P4rest( P3, P3.mag() );
    P4rest.boost(CMBoost);
    Hep3Vector P3_cm = P4rest.vect();
    P_sum += P3_cm;
    
    C_sum_L += TrkList_TwoPhoLevel1[i]->charge();
    
    if( TrkList_TwoPhoLevel1[i]->charge() > 0.0 ) ntr1_p++;
    else ntr1_n++;           //U1-1
    
  }

  //*---
  // level-3 track.
  //*---
  for( int i=0;i<ntr3;i++ ){
    C_sum += TrkList_TwoPhoLevel3[i]->charge();
    
    Hep3Vector P3( TrkList_TwoPhoLevel3[i]->px(),
                   TrkList_TwoPhoLevel3[i]->py(),
                   TrkList_TwoPhoLevel3[i]->pz() );
    HepLorentzVector P4rest( P3, P3.mag() );
    P4rest.boost(CMBoost);
    P4_sum_tr3 += P4rest;             //U1-2, U1-3
    
  }
  
  //Added selection criteria by S.Uehara 7-Mar-2000
  if( !( ntr1_p >0 && ntr1_n >0 ) ) return;  // U1-1
  if( P4_sum_tr3.mag() > tp_wmax ) return;   // U1-2
  HepLorentzVector P4_miss = HepLorentzVector(0.,0.,0.,cm.mag())
    - P4_sum_tr3;
  if( P4_miss.mag2() < tp_mm2min ) return;    // U1-3

  //*---
  // event classification starts here
  //   using global variables.
  //     ntr0, ntr1, ntr2, ntr3, esum, psum.
  //*---
  int evcflag = 0;
  
  float evis = esum + psum;
  
  if( psum < 6.0 && esum < 6.0 && ntr3 >=2 && ntr3<=4  ){
    // low-mul events
    if( ntr3 == 2 && abs(C_sum)<0.01 ){
      if( ntr1 == 2 ) {
        if( P_sum.perp() < 0.2 ) {
          // two-track exclusive (11)
          evcflag = 11;
        }
      }
    }
    
    if( ntr1 == 4 && abs(C_sum_L)<0.01 ) {
      // exclusive 4-track events with tight vertex (110)
      evcflag = 110;
    }
    
    if( evcflag == 0 ) evcflag = 10;   // other low-mul events (10)
  }
  
  if(  evis < 4.0 && ntr3 >=2 ) {
    if( evcflag == 0 ) evcflag = 20;
  }
  
  if( psum < 6.0 && esum < 6.0 && ntr1 == 4
      && abs(C_sum_L)<0.01 ) {
    // 4-track exclusive with loose vertex (100)
    if( evcflag == 0 ) evcflag = 100;
  }


  //U2-1) Require either of Charge sums = 0 for evcflg=10, 20
 if( evcflag == 10 || evcflag == 20 )
    if( abs(C_sum_L) > 0.01 && abs(C_sum) > 0.01 ) return;

 //U2-2) Require minimum egam (0.1GeV) for clusters above the 
 // prefixed threshold (40MeV) for 2-track events of evcflg=10
    if( evcflag == 10 && ntr1 == 2 )
      if(egamat < 0.1 ) return;
         
  //*---
  // now classification completed. try to fill panther tables.
  //*---
  if( evcflag == 0 && output_level == 0 ) return;
  
  evtcls_flag[5] = evcflag;
  


  //++++
  //(1) fill event global informations.
  //++++
  Evtcls_twophoton_info_Manager & EvtInfoMgr
    = Evtcls_twophoton_info_Manager::get_manager();
  Evtcls_twophoton_info& evtinfo = EvtInfoMgr.add();
  
  evtinfo.Ntrk_level1(ntr1);
  evtinfo.Ntrk_level2(ntr2);
  evtinfo.Ntrk_level3(ntr3);
  evtinfo.Psum(psum);
  evtinfo.Esum(esum);
  evtinfo.Qsum(C_sum_L);
  evtinfo.Ptsum(P_sum.perp());
  
  //++++
  //(2) fill pointers for each level track to MDST_Charged.
  //++++
  Evtcls_twophoton_level1_Manager & EvtChLv1Mgr
    = Evtcls_twophoton_level1_Manager::get_manager();
  if( ntr1 ){
    for( int i=0; i<ntr1; ++i ){
      Evtcls_twophoton_level1& evtchg = EvtChLv1Mgr.add();
      evtchg.charged( *TrkList_TwoPhoLevel1[i] );
    }
  }else{
    Evtcls_twophoton_level1& evtchg = EvtChLv1Mgr.add();
    evtchg.reset_charged();
  }
  
  Evtcls_twophoton_level2_Manager & EvtChLv2Mgr
    = Evtcls_twophoton_level2_Manager::get_manager();
  if( ntr2 ){
    for( int i=0; i<ntr2; ++i ){
      Evtcls_twophoton_level2& evtchg = EvtChLv2Mgr.add();
      evtchg.charged( *TrkList_TwoPhoLevel2[i] );
    }
  }else{
    Evtcls_twophoton_level2& evtchg = EvtChLv2Mgr.add();
    evtchg.reset_charged();
  }
  
  Evtcls_twophoton_level3_Manager & EvtChLv3Mgr
    = Evtcls_twophoton_level3_Manager::get_manager();
  if( ntr3 ){
    for( int i=0; i<ntr3; ++i ){
      Evtcls_twophoton_level3& evtchg = EvtChLv3Mgr.add();
      evtchg.charged( *TrkList_TwoPhoLevel3[i] );
    }
  }else{
    Evtcls_twophoton_level3& evtchg = EvtChLv3Mgr.add();
    evtchg.reset_charged();
  }
  
}

// cosmic event selection.
void EventClassify::sel_cosmic ( int* evtcls_flag, int output_level )
{

  //*---
  // reset local variables.
  //*---
  float Psave[2][5];
  float Tsave[5];
  int Isave[2][5];
  int track_flag = 0;
  int tof_flag   = 0;
  for( int j=0;j<5;j++ ){
    Tsave[j]=0.0;
    for( int i=0;i<2;i++ ){
      Psave[i][j]=0.0;
      Isave[i][j]=0;
    }
  }

  //*---
  // Get managers for MDST_Trk, MDST_TOF & MDST_ECL.
  //*---
  Mdst_trk_Manager& TrkMgr = Mdst_trk_Manager::get_manager();
  Mdst_trk_fit_Manager& TrkfitMgr = Mdst_trk_fit_Manager::get_manager();
  Mdst_tof_Manager& ToFMgr = Mdst_tof_Manager::get_manager();

  //*---
  // # of charged tracks coming from off-vertex tracks.
  //   off-vertex:R greater than radius of SVD outer layer
  //              Z is out of region of Zmin <Z< Zmax.
  //*---
  int Ncos = TrkList_offvtx.length();

  //*---
  // return if Ncos <= 1.
  //*---
  //  if( Ncos <= 1 ) return;
  if( Ncos != 2 ) return;

  //*---
  // look at difference between closest approaches
  //   & TOF(if available)
  //   in case of Ncos>=2.
  //*--- 
  for( int i=0; i<Ncos; i++ ){
    Mdst_trk& Trk1 = TrkList_offvtx[i]->trk();
    Mdst_trk_fit& Trkfit1  = Trk1.mhyp(1);
    Mdst_tof& Tof1 = TrkList_offvtx[i]->tof();
    float dr1 = Trkfit1.helix(0);
    float dz1 = Trkfit1.helix(3);
    Hep3Vector P1(TrkList_offvtx[i]->px(),TrkList_offvtx[i]->py(),
		  TrkList_offvtx[i]->pz());

    for( int j=i+1; j<Ncos; j++ ){

      Mdst_trk& Trk2 = TrkList_offvtx[j]->trk();
      Mdst_trk_fit& Trkfit2  = Trk2.mhyp(1);
      Mdst_tof& Tof2 = TrkList_offvtx[j]->tof();
      float dr2 = Trkfit2.helix(0);
      float dz2 = Trkfit2.helix(3);
      Hep3Vector P2(TrkList_offvtx[j]->px(),TrkList_offvtx[j]->py(),
		    TrkList_offvtx[j]->pz());

      // dout(Debugout::DDEBUG,"EventClassify") << " diff dr : " <<  abs((abs(dr1)-abs(dr2))) <<std::endl;
      // dout(Debugout::DDEBUG,"EventClassify") << " diff dz : " <<  abs(dz1-dz2) <<std::endl;
      // dout(Debugout::DDEBUG,"EventClassify") << "  p1-z   : " << P1.z() << "  p2-z : "<< P2.z()<<std::endl;
      // dout(Debugout::DDEBUG,"EventClassify") << " ratio   : " << ( P1.z()+P2.z() )/P1.mag()  << std::endl;

      //look at difference 
      //  between closest approaches of two tracks.
      if( abs((abs(dr1)-abs(dr2)) )<= Rimpact_cut ){
	if( abs(dz1-dz2) <= Zimpact_cut ){

	  //check momnetum balance in Z-component.
	  //it should be balanced if cosmic.
	  if( abs( (P1.z()+P2.z())/P1.mag() ) <= Pz_cosmic_cut ) {
	    track_flag++;
	    Isave[0][track_flag-1]=i;
	    Isave[1][track_flag-1]=j;
	    Psave[0][track_flag-1]=P1.mag();
	    Psave[1][track_flag-1]=P2.mag();
	  }

	  //look at difference 
	  //  between ToF countes' hits.
	  if( Tof1 ){
	    if( Tof2 ){
	      float time_diff = abs( Tof1.tof()-Tof2.tof() );
	      if( time_diff <= tof_time_diff_cut ){
		tof_flag++;
		Tsave[tof_flag-1]=time_diff;
	      }
	    }
	  }
	}
      }
    }
  }

  //*---
  // # of ECL good clusters.
  //*---
  int Necl = EclListCM.length();

  float Esum = 0.0;
  //*---
  // Cluster energy sum.
  //*---  
  for(int iecl=0;iecl<Necl;iecl++){
    Esum += EclListCM[iecl]->e();
  }

  //*---
  // cosmic event classification using global variables obtained above.
  //*---
  if( (Ncos>=2) && (Esum<=Esum_MIP_cut) ){
    if( track_flag == 1){
      if( tof_flag == 1 ){
	evtcls_flag[6] = 10;
      }else if( tof_flag == 0 ){
	evtcls_flag[6] = 11;
      }
    }
  }

  // dout(Debugout::DDEBUG,"EventClassify") << " Esum : " << Esum <<std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << " track_flag : " << track_flag <<std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << " tof_flag : " << tof_flag <<std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << " cos_flag : " << evtcls_flag[6] <<std::endl;
    
  //*---
  // Fill panther table if classified as cosmic event.
  //*---
  if( output_level == 0 ){
    if( evtcls_flag[6] == 0 ) return;
  }

  //output for debugging.
  // dout(Debugout::DDEBUG,"EventClassify") << " cos_flag : " << cos_flag  <<std::endl;
  //++++
  //(1) fill event_flag.
  //++++
  //  Evtcls_flag_Manager& EvtFlgMgr = Evtcls_flag_Manager::get_manager();
  //  Evtcls_flag& evtflag = EvtFlgMgr.add();
  //  evtflag.flag(6,cos_flag);

  //++++  
  //(2) fill cosmic event global informations.
  //++++  
  Evtcls_cosmic_info_Manager & EvtInfoMgr 
    = Evtcls_cosmic_info_Manager::get_manager();
  Evtcls_cosmic_info& evtinfo = EvtInfoMgr.add();
  
  evtinfo.Ntrk(Ncos);
  evtinfo.Ncls(Necl);
  evtinfo.Esum(Esum);
  if( Psave[0][0] > Psave[1][0] ){
    evtinfo.P1st(Psave[0][0]);
    evtinfo.P2nd(Psave[1][0]);
  }else{
    evtinfo.P1st(Psave[1][0]);
    evtinfo.P2nd(Psave[0][0]);
  }
  if( tof_flag ){
    evtinfo.ToF_diff(Tsave[0]);
  }else{
    evtinfo.ToF_diff(-99.0);
  }

  // dout(Debugout::DDEBUG,"EventClassify") << " Ncos : "<< Ncos <<std::endl;
  // dout(Debugout::DDEBUG,"EventClassify") << " Isave[0][0] : "<<Isave[0][0]<<" Isave[1][0]: "<<Isave[1][0]<<std::endl;

  //++++  
  //(3) fill pointers for cosmic tracks to MDST_Charged.
  //++++  
  Evtcls_cosmic_charged_Manager & EvtCosMgr 
    = Evtcls_cosmic_charged_Manager::get_manager();
  if( Isave[0][0] && Isave[1][0] ){
    for( int i = 0; i<2; ++i ){
      Evtcls_cosmic_charged& evtchg = EvtCosMgr.add();
      evtchg.charged( *TrkList_offvtx[Isave[i][0]] );
    }
  }else{
    Evtcls_cosmic_charged& evtchg = EvtCosMgr.add();
    evtchg.reset_charged();
  }

}

// radiative Bhabha event selection.
void EventClassify::rad_bhabha ( int* evtcls_flag, int output_level )
{

  //*---
  // reset local flag.
  //*---
  int rad_bhabha = 0;

  //*---
  // return if # of track == 0.
  //*---
  if( TrkList.length() < Ntrk_eeg_cut ) return;

  //*---
  // # of track should be even. add on '99-Dec-13. IA.
  //*---
  if( (TrkList.length()%2 ) != 0 ) return;

  //*---
  // return if # of ECL clusters < 3.
  //*---
  if( EclList.length() < Ncls_eeg_cut ) return;

  //*---
  // smallest energy cluster > 100MeV.  add on '99-Jul-27.
  //*---
  if( EclList[2]->energy() < E_eeg_min_cut ) return;

  //*---
  // compute the most energitic 3 clusters.
  //*---
  float E3sum = EclList[0]->energy() + EclList[1]->energy()
    + EclList[2]->energy();

  //*---
  // check total charge of tracks. add on '99-Dec-13. IA.
  //*---
  float Qsum = 0.0;
  for( int i=0;i<TrkList.length();i++ ){
    Qsum += TrkList[i]->charge();
  }

  //*---
  // check total energy deposit in ECL. add on '99-Dec-13. IA.
  //*---
  float Esum = 0.0;
  for( int i=0;i<EclListAll.length();i++ ){
    Esum += EclListAll[i]->energy();
  }
    
  //*---
  // selection here.
  //*---
  if( E3sum >= E3sum_eeg_cut ){
    if( Qsum == Q_eeg_cut ){
      if( Esum >= E_eeg_total_cut ){
	rad_bhabha = 10;
      }else{
	rad_bhabha = 20;
      }	
    }
  }

  //*---
  // return if event not pass thru cuts.
  //*---
  if( output_level == 0 ){
    if( rad_bhabha == 0 ) return;
  }

  evtcls_flag[7] = rad_bhabha;

  //*---
  // put event info. into panther table.
  //*---
  Evtcls_radbha_info_Manager & EEGinfoMgr 
    = Evtcls_radbha_info_Manager::get_manager();
  Evtcls_radbha_info& evtinfo = EEGinfoMgr.add();
  evtinfo.Ncls( (int) EclList.length() );
  evtinfo.Ntrk( (int) TrkList.length() );
  for( int i=0;i<3;i++ ) 
    evtinfo.Ecls(i,EclList[i]->energy());

  //*---
  // put pointers into panther table.
  //*---
  Evtcls_radbha_cluster_Manager & EvtEEGMgr 
    = Evtcls_radbha_cluster_Manager::get_manager();
  if( EclList.length() ){
    for( int i = 0; i<3; ++i ){
      Evtcls_radbha_cluster& evtchg = EvtEEGMgr.add();
      evtchg.ecl( *EclList[i] );
    }
  }else{
    Evtcls_radbha_cluster& evtchg = EvtEEGMgr.add();
    evtchg.reset_ecl();
  }
  

}

// radiative mu-pair event selection
void EventClassify::radiative_mu ( int* evtcls_flag, int output_level )
{
  //*---
  //* evtcls_flag[11] assigned.
  //*---
  //*---
  // return if one track is postive and the other negative.
  //*---
  if( TrkList_pos_rest.length() != N_radmu_trk_cut ) return;
  if( TrkList_neg_rest.length() != N_radmu_trk_cut ) return;

  //*---
  // look the two track's momentum at the rest frame.
  //*---

  Hep3Vector P3( TrkList_pos_rest[0]->px(), TrkList_pos_rest[0]->py(),
                 TrkList_pos_rest[0]->pz() );
  HepLorentzVector P1CM( P3, sqrt(P3.mag2()+ 0.1056*0.1056 ) );
  HepLorentzVector P1LAB(P1CM);
  P1CM.boost(CMBoost);
  Hep3Vector P3Rest1 = P1CM.vect();

  P3=Hep3Vector( TrkList_neg_rest[0]->px(), TrkList_neg_rest[0]->py(),
                 TrkList_neg_rest[0]->pz() );
  HepLorentzVector P2CM( P3, sqrt(P3.mag2()+ 0.1056*0.1056 ) );
  HepLorentzVector P2LAB(P2CM);
  P2CM.boost(CMBoost);
  Hep3Vector P3Rest2 = P2CM.vect();

  //*---
  // return if both momentum > 4.8GeV.
  //*---
  //GM  if( P3Rest1.mag() >= P_radmu_cm_cut || P3Rest2.mag() >= P_radmu_cm_cut ) return;

  if( P3Rest1.mag()+P3Rest2.mag() >= 2*P_radmu_cm_cut ) return;

  //*---
  // TOF time difference if exits.
  //*---
  float trk_tof1, trk_tof2;
  Mdst_tof& tof1 = TrkList_pos_rest[0]->tof();
  Mdst_tof& tof2 = TrkList_neg_rest[0]->tof();
  if( tof1 ){
    trk_tof1 = tof1.tof();
  }else{
    trk_tof1 = -99.0;
  }
  if( tof2 ){
    trk_tof2 = tof2.tof();
  }else{
    trk_tof2 = -99.0;
  }
  if( trk_tof1 != -99.0 && trk_tof2 != -99.0 ){
    if(  abs( trk_tof1 - trk_tof2 ) >= tof_radmu_time_diff_cut ) return;
  }

  //*---
  // # of ECL good clusters & total energy.
  //*---
  //* 06/09/2003 GM Replace EclListCM by EclListCMSel
  float ECLsum = 0.0;
  for( int icls=0; icls<EclListCMSel.length(); icls++ ){
    ECLsum += EclListCMSel[icls]->e();
  }

  //*---
  // # of gammas & gamma energy sum to split two cases.
  //*---
  int Ngamma = GamList_rest.length();
  float Egamma = 0.0;
  float EGamMX = 0.0;
  for( int ig=0; ig<Ngamma; ig++ ){
    Hep3Vector P3( GamList_rest[ig]->px(), GamList_rest[ig]->py(),
		   GamList_rest[ig]->pz() );
    if (P3.mag() > EGamMX) EGamMX = P3.mag();
    HepLorentzVector PgamCM( P3, P3.mag());
    PgamCM.boost(CMBoost);
    Egamma += PgamCM.e();
  }

  // G. Majumder 04 Oct 2004
  // Selection criteria without any requirement of observed photon
  // which is based on the criteria of calibatration of high energy
  // photon and its selection efficiency.
  //

  if (EclListCMSel.length() <= N_radmu_ecl_case1 &&
      ECLsum < ECL_radmu_max_case1 && 
      Egamma < Gamma_radmu_max_case1){

    HepLorentzVector cmE(EHER*sin(theta), 0., 
			 EHER*cos(theta)-ELER,
			 EHER+ELER );

    HepLorentzVector exp4v = cmE - P1LAB - P2LAB;
    if ( abs(P2LAB.theta()-0.555) >0.03 &&
	 abs(P1LAB.theta()-2.26) >0.03  &&
	 exp4v.rho() > 0.50 && 
	 abs(exp4v.e()-exp4v.rho()) < 0.7  && 
	 exp4v.theta() > 0.25 && 
	 exp4v.theta() < 2.65 ) {
      
      evtcls_flag[11] = 11;
    }
  }

  //*---
  // case-a -- Gamma's exist.
  //*---
  if( Egamma > Gamma_radmu_min_case1 ){

    if( EclListCMSel.length() <= N_radmu_ecl_case1 ){

      if( ECLsum > ECL_radmu_min_case1 &&
	  ECLsum < ECL_radmu_max_case1 && 
	  EGamMX > Gamma_radmu_min_case1 ){
	
	if( Egamma > Gamma_radmu_min_case1 && 
	    Egamma < Gamma_radmu_max_case1 ){

	  float Ecm_total = P1CM.e()+P2CM.e()+EclListCMSel[0]->e();

	  //Angle criteria add on 12/11/2004 by GM to remove
	  //Bhabha events,where one track passed through the ECL gap
	  //remove only electron in forward side and positron in 
	  //backward side 

	  if( Ecm_total > E_radmu_total_cut_case1 &&
	      abs(P2LAB.theta()-0.555) >0.03 &&
	      abs(P1LAB.theta()-2.26) >0.03 ){
	    
	    evtcls_flag[11] = 10;
	    
	  }
	}
      }
    }
  }else{
    //*---
    // case-b -- No gamma case. 
    //*---  
    // G. Majumder 12/11/2004
    //obsulate : no event passed through 
    //       EclListCM.length() <= N_radmu_ecl_case2 criteria
    //  in general two muon associated ecllist > N_radmu_ecl_case2(1)

    if( EclListCM.length() <= N_radmu_ecl_case2 ){
      if( ECLsum < ECL_radmu_max_case2 ){
    
	//*---
	// look at two track info. in detail.
	//*---
	float acol = P3Rest1.angle( -P3Rest2 )*180.0/acos(-1.0);
	Hep3Vector P2Lab1( TrkList_pos_rest[0]->px(), TrkList_pos_rest[0]->py(),
			   0.0 );
	Hep3Vector P2Lab2( TrkList_neg_rest[0]->px(), TrkList_neg_rest[0]->py(),
			   0.0 );
	float acop = P2Lab1.angle( P2Lab2 )*180.0/acos(-1.0);
	Hep3Vector PsumLab = P2Lab1 + P2Lab2;
	
	if( acol > acol_radmu_cut ){
	  if( acop > acop_radmu_cut ){
	    if( PsumLab.mag() < P_radmu_lab_cut ){

	      evtcls_flag[11] = 12;
	      
	    }
	  }
	}
      }
    }
  }

  //*---
  // that's it.
  //*---

}

// 4e final state events selection.
void EventClassify::eeee_sel ( int* evtcls_flag, int output_level )
{

  //*---
  // local flag reset.
  //*---
  int eeee_flag = 0;

  int N_good_trk_tmp = TrkList_eeee.length();

  //*---
  // return if N_good_trk_tmp <= 1.
  //*---
  if( N_good_trk_tmp <= 1 ) return;

  //*---
  // return if # ECL cluster <= 1.
  //*---
  if ( EclListAll.length() <= 1 ) return;

  //*---
  // initialize.
  //*---
  int N_good_pos_trk = 0;
  int N_good_neg_trk = 0;
  int N_good_trk = 0;
  int N_ext_trk = 0;
  int N_good_trk_2 = 0;
  int N_good_trk_3 = 0;
  int N_e_cand = 0;
  float  c_tot = 0;

  HepLorentzVector p_tot(0.,0.,0.,0.);
  float vz_center = 0.;
  float *charge = new float[N_good_trk_tmp];
  float *vertex_z = new float[N_good_trk_tmp];
  float *vertex_pos_z = new float[N_good_trk_tmp];
  float *vertex_neg_z = new float[N_good_trk_tmp];
  for ( int i=0; i< N_good_trk_tmp; i++){
    vertex_z[i] = TrkList_eeee[i]->trk().mhyp(2).helix(3);
    charge[i] = TrkList_eeee[i]->charge();
    c_tot += charge[i];

    //*---
    // loop if level == 1.
    //*---
    if ( Vz_4e_cut_level==1 ){
      vz_center += vertex_z[i];
    }else if ( Vz_4e_cut_level==2 ){
      if ( charge[i] > 0){
        vertex_pos_z[N_good_pos_trk] = vertex_z[i];
        N_good_pos_trk++;
      }else{
        vertex_neg_z[N_good_neg_trk] = vertex_z[i];
        N_good_neg_trk++;
      }
    }
  }

  if ( c_tot > 1.) return;
  if ( N_good_trk_tmp ) vz_center /= (float)N_good_trk_tmp;

  //*---
  // loop if level == 2.
  //*---
  if ( Vz_4e_cut_level==2 ){
    float vz_min_diff = 9999.;
    for ( int i=0; i< N_good_pos_trk; i++){
      for ( int j=0; j< N_good_neg_trk; j++){
        float vz_diff = abs(vertex_pos_z[i] - vertex_neg_z[j]);
        if ( vz_diff < vz_min_diff ){
          vz_min_diff = vz_diff;
           vz_center = 0.5 * (vertex_pos_z[i] + vertex_neg_z[j]);
        }
      }
    }
  }

  c_tot = 0.;
  atc_pid atcPid(0,-1,0,0,2);
  for ( int i=0; i< N_good_trk_tmp; i++){
    if ( Vz_4e_cut_level!=0 ){
      if ( abs(vz_center) > (Vz_4e_cut-Vz_4e_cut2) ) continue;
      if ( abs(vertex_z[i]-vz_center) > Vz_4e_cut2 ) continue;
    }

    Hep3Vector p3_lab(TrkList_eeee[i]->px(),
		      TrkList_eeee[i]->py(),TrkList_eeee[i]->pz());
    HepLorentzVector p_lab( p3_lab, sqrt(p3_lab.mag2()+e_mass*e_mass));
    N_good_trk++;
    p_tot += p_lab;
    c_tot += charge[i];
    if ( p3_lab.perp() < Pt_4e_cut ) continue;

    int mu_flag = 
      (TrkList_eeee[i]->klm()) ? TrkList_eeee[i]->klm().muon() : 0;
    if (mu_flag) continue;
    N_ext_trk++;
    float p3_mag = p3_lab.mag();
    if (p3_mag > 0.5) N_good_trk_2++;
    if (p3_mag > 0.95) N_good_trk_3++;
    if (!atcPid.best_id(TrkList_eeee[i])) N_e_cand++;
  }

  delete [] charge;
  delete [] vertex_z;
  delete [] vertex_pos_z;
  delete [] vertex_neg_z;

  //*---
  // final selection here.
  //*----
  if ( N_good_trk != 2 ) return;
  if ( c_tot != 0. ) return;
  if ( !N_e_cand ) return;
  p_tot.boost( CMBoost );
  if ( S_Pt_4e_cut >0 && p_tot.perp() > S_Pt_4e_cut ) return;
  if ( W_4e_cut > 0 && p_tot.mag() > W_4e_cut ) return;
  if ( S_Pz_4e_cut > 0 && abs(p_tot.z()) > S_Pz_4e_cut ) return;
  if ( N_ext_trk != 2 ) return;
  if ( !N_good_trk_2 ) return;

  float E_tot = 0.;
  for ( int i=0; i<EclListAll.length(); i++ ){
    E_tot += EclListAll[i]->energy();
  }

  if (E_tot_4e_max > 0 && E_tot > E_tot_4e_max) return;
  if (E_tot_4e_min > 0 && E_tot < E_tot_4e_min) return;

  if ( N_good_trk_3 ) eeee_flag = 20;
  eeee_flag = 10;

  //*---
  // store flag.
  //*---
  evtcls_flag[8] = eeee_flag;

  return;

}

// Gamma+Phi event selection.
void EventClassify::GammaPhi_sel ( int* evtcls_flag, int output_level )
{

  //*---
  // reset.
  //*---
  int sel_loose = 0;
  int sel_tight = 0;
  int sel_kid   = 0;

  //*---
  // return if # cluster < 1.
  //*---
  if( GamList_gphi.length() < N_gphi_cls_cut ) return;

  //*---
  // return if # track !=2.
  //*---
  if( TrkList_gphi.length() != N_gphi_trk_cut ) return;

  //*---
  // charge should be opposite.
  //*---
  float Q1 = TrkList_gphi[0] -> charge();
  float Q2 = TrkList_gphi[1] -> charge();
  if( Q1*Q2 >= 0.0 ) return;

  //*---
  // there should be one cluster with
  // energy at the rest frame : 4.5 < E < 5.5 GeV.
  //*---
  Hep3Vector v0( GamList_gphi[0]->px(), GamList_gphi[0]->py(), 
		 GamList_gphi[0]->pz() );
  HepLorentzVector v0rest4( v0, v0.mag() );
  v0rest4.boost(CMBoost);
  Hep3Vector GamRest0 = v0rest4.vect();

  if( GamRest0.mag() < E_gphi_min_cut || E_gphi_max_cut < GamRest0.mag() )
    return;

  if( GamList_gphi.length() > N_gphi_cls_cut ){

    Hep3Vector v1( GamList_gphi[1]->px(), GamList_gphi[1]->py(), 
		   GamList_gphi[1]->pz() );
    HepLorentzVector v1rest4( v1, v1.mag() );
    v1rest4.boost(CMBoost);
    Hep3Vector GamRest1 = v1rest4.vect();

    if( GamRest1.mag() > E_gphi_min_cut )
      return;
  }

  //*---
  // look at Kshort info.
  //*---
  int  ks_exist = 0;
  Mdst_vee_Manager&  vmgr = Mdst_vee_Manager::get_manager();
  Mdst_vee2_Manager&  v2mgr = Mdst_vee2_Manager::get_manager();
  for( std::vector<Mdst_vee>::iterator v = vmgr.begin();
       v != vmgr.end(); v++ )
    if ( v->kind() == 1 ){ // 1=Ks 2=Lambda 3=Lambda-bar 4=converted-gamma
      ks_exist++;
    }
  for( std::vector<Mdst_vee2>::iterator v = v2mgr.begin();
       v != v2mgr.end(); v++ )
    if ( v->kind() == 1 ){ // 1=Ks 2=Lambda 3=Lambda-bar 4=converted-gamma
      ks_exist++;
    }
  
  //*---
  // look at elecron-ID info.
  //*---
  for( int i=0; i<2; i++ ){
    
    eid  elid( *TrkList_gphi[i] );
    float  eprob = elid.prob( 0,-1,0 );
    if ( eprob < eprob_gphi_tight_cut )
      sel_tight++;
  }

  //*---
  // look at kaon-ID info.
  //*---
  for( int i=0; i<2; i++ ){
    atc_pid  Kpi( 0, 1, 0, 3, 2 );// 0:e 1:mu 2:pi 3:K 4:p
    float kprob = Kpi.prob( TrkList_gphi[i] );
    
    if ( kprob > kid_gphi_cut ) 
      sel_kid++;
  }

  //*---
  // select event.
  //*---
  if ( ( ks_exist == 1 ) ||                    // for KsKL
       ( sel_tight >= 2 && sel_kid >= 2 ) ){   // for K+k-
    evtcls_flag[10] = 10;
  }


}

// Beam background selection.
void EventClassify::beamBGsel ( int* evtcls_flag, int output_level )
{

  //*---
  // if mu-pair flag already fired, return.
  //*---
  if( evtcls_flag[3] ) return;

  //*---
  // reset local variable.
  //*---
  float Esum = 0.0;

  //*---
  // check total energy deposit in ECL.
  //*---
  for( int i=0;i<EclListAll.length();i++ ){
    Esum += EclListAll[i]->energy();
  }

  //*---
  // crude selection. should be modified.
  //*---
  if( Esum < Esum_BG_cut ){
    evtcls_flag[18] = 10;
  }


}

// Monitoring sample for beam gas events.
// hadronic event selection.
void EventClassify::BeamGasMon_sel ( int* evtcls_flag, int output_level )
{

  //*----
  //(Note)
  // events are selected using the same cuts as HadronA.
  // (except for primary vertex cuts)
  //*---

  //*---
  // reset global variables here.
  //*---
  //  for charged tracks.
  float Psum     = 0.0;
  float Esum_chg = 0.0;
  float Psum_cms     = 0.0;
  float Esum_chg_cms = 0.0;
  float Pzsum    = 0.0;
  //  for neutral clusters.
  float Esum_cls     = 0.0;
  float Esum_cls_cms = 0.0;
  //  for ECL clusters.
  float ECL_sum = 0.0;


  //*---
  // calculate charged observables.
  //*---
  int Ntrk = TrkListBG.length();

  // loop over # of good tracks.  
  for( int itrk=0; itrk< Ntrk; itrk++ ){

    Hep3Vector P (TrkListBG[itrk]->px(),TrkListBG[itrk]->py(),
		  TrkListBG[itrk]->pz());
    double E =sqrt(P.mag2()+ pow(pi_mass,2.));
    HepLorentzVector P_cms(P,E);
    P_cms.boost(CMBoost);
    Psum += P.mag();
    Esum_chg += E;
    float ppcms = sqrt( pow(P_cms.e(),2.0) - P_cms.mag2() );
    Psum_cms += sqrt( pow(P_cms.e(),2.0) - P_cms.mag2() );
    Esum_chg_cms += P_cms.e();
    Pzsum += P_cms.pz();

  }

  //*---
  // calculate neutral observables.
  //*---
  int Ncls = GamList.length();

  // loop over # of good clusters.
  for( int icls =0; icls<Ncls; icls++ ){

    Hep3Vector P (GamList[icls]->px(),GamList[icls]->py(),
		  GamList[icls]->pz());
    double E = P.mag();
    Esum_cls += E;
    HepLorentzVector P_cms(P,E);
    P_cms.boost(CMBoost);
    Esum_cls_cms += P_cms.e();
    Pzsum += P_cms.pz();
 
  }
  
  //*---
  // compute visible energy. 
  //*---
  float Evis     = Esum_cls + Esum_chg;
  float Evis_cms = Esum_cls_cms + Esum_chg_cms;

  //*---
  // ECL deposit energy.
  //*---
  for( int i=0;i<EclListCMSel.length();i++){
    Hep3Vector Ecls = EclListCMSel[i]->vect();
    ECL_sum += Ecls.mag();
  }

  //*---
  // Beam-gas event classification.
  //*---
  if( Ntrk >= Ntrk_cut ){
    if( (Evis_cms/ebeam) >= Ev_cut ){
      if( abs(Pzsum/ebeam) <= Pz_cut ) {
	if( (ECL_sum/ebeam) >= Esum_min_cut &&
	    (ECL_sum/ebeam) <= Esum_max_cut ) {
	  
	  evtcls_flag[17] = 10;
	  
	}
      }
    }
  }
  
  
}


// put evtcls_flag into panther table.
void EventClassify::put_flag ( int* evtcls_flag, int* hadronic_flag )
{

  //*---
  // unknown/junk events assgined if no flag put so far.
  //*---
  int flag_sum = 0;
  for( int i=0;i<19;i++){
    flag_sum += evtcls_flag[i];
  }
  if( flag_sum == 0 ) evtcls_flag[19] = 10;

  //*---
  // fill event_flag.
  //*---
  Evtcls_flag_Manager& EvtFlgMgr = Evtcls_flag_Manager::get_manager();
  Evtcls_flag& evtflag = EvtFlgMgr.add();
  Evtcls_flag2_Manager& EvtFlg2Mgr = Evtcls_flag2_Manager::get_manager();
  Evtcls_flag2& evtflag2 = EvtFlg2Mgr.add();

  // Loop over # of event types.
  //                 currently 20 types.
  //                =0:hadron /=1:Bhabha /=2:gamma-pair /=3:mu-pair
  //                =4:tau-pair /=5:two photon /=6:cosmic 
  //                =7:radiative Bhabha(add on '99-May-26 IA)
  //                =8:eeee(add on '99-May-26 IA)
  //                =9:calorimter Bhabha(add on '99-Jul-09 IA)
  //                =10:gamma+phi event(add on '99-Aug-20 IA)
  //                =11:radiative mu-pair event('99-Oct-14/not available now IA)
  //                =17:beam gas monitoring event(add on '99-Jul-07 IA)
  //                =18:beam background /=19:unknown or junk
  for( int i=0;i<10;i++){
    evtflag.flag(i,evtcls_flag[i]);
  }
  for( int i=10;i<20;i++){
    evtflag2.flag(i-10,evtcls_flag[i]);
  }

  //*---
  // fill hadronic_flag.
  //*---
  Evtcls_hadronic_flag_Manager& HadMgr = Evtcls_hadronic_flag_Manager::get_manager();
  Evtcls_hadronic_flag& evt_had_flag = HadMgr.add();
  for( int i=0;i<20;i++){
    evt_had_flag.hadronic_flag(i,hadronic_flag[i]);
  }


}
// clear event data.
void EventClassify::clear ( void )
{

  //*---
  //delete all STL's & HepALists.
  //*---
  //  TrkList.erase( TrkList.begin(), TrkList.end() );
  //  TrkList_offvtx.erase( TrkList_offvtx.begin(), TrkList_offvtx.end() );
  //  GamList.erase( GamList.begin(), GamList.end() );
  //  EclList.erase( EclList.begin(), EclList.end() );

  vlist.erase( vlist.begin(), vlist.end() );
  HepAListDeleteAll( v4list );
  HepAListDeleteAll( EclListCM );
  HepAListDeleteAll( EclListCMSel );
  HepAListDeleteAll( EclListCMAll );
  HepAListDeleteAll( EclListLab );
  HepAListDeleteAll( EclListTrack );  

  Mdst_charged * e;
  while ( ( e = TrkList.last() ) ) {
      TrkList.remove(e);
  }
  while ( ( e = TrkList_offvtx.last() ) ) {
      TrkList_offvtx.remove(e);
  }  
  while ( ( e = TrkList_pos_rest.last() ) ) {
      TrkList_pos_rest.remove(e);
  }
  while ( ( e = TrkList_neg_rest.last() ) ) {
      TrkList_neg_rest.remove(e);
  }
  while ( ( e = TrkList_pos.last() ) ) {
      TrkList_pos.remove(e);
  }
  while ( ( e = TrkList_neg.last() ) ) {
      TrkList_neg.remove(e);
  }
  while ( ( e = TrkList_pair.last() ) ) {
      TrkList_pair.remove(e);
  } 
  while ( ( e = TrkList_pair_rest.last() ) ) {
      TrkList_pair_rest.remove(e);
  }
  while ( ( e = TrkList_TwoPho.last() ) ) {
      TrkList_TwoPho.remove(e);
  }
  while ( ( e = TrkList_TwoPhoLevel1.last() ) ) {
      TrkList_TwoPhoLevel1.remove(e);
  }
  while ( ( e = TrkList_TwoPhoLevel2.last() ) ) {
      TrkList_TwoPhoLevel2.remove(e);
  }
  while ( ( e = TrkList_TwoPhoLevel3.last() ) ) {
      TrkList_TwoPhoLevel3.remove(e);
  }
  while ( ( e = TrkList_eeee.last() ) ) {
      TrkList_eeee.remove(e);
  }
  while ( ( e = TrkListBG.last() ) ) {
      TrkListBG.remove(e);
  }
  while ( ( e = TrkListTau.last() ) ) {
      TrkListTau.remove(e);
  }
  while ( ( e = TrkList_gphi.last() ) ) {
      TrkList_gphi.remove(e);
  }

  Mdst_gamma * f;
  while ( ( f = GamList.last() ) ) {
      GamList.remove(f);
  }  
  while ( ( f = GamList_rest.last() ) ) {
      GamList_rest.remove(f);
  }  
  while ( ( f = GamList_gphi.last() ) ) {
      GamList_gphi.remove(f);
  }  

  Mdst_ecl * g;
  while ( ( g = EclList.last() ) ) {
      EclList.remove(g);
  }
  while ( ( g = EclListAll.last() ) ) {
      EclListAll.remove(g);
  }
  
}

// check if do event classification or not
void EventClassify::check_event ( int limit_mdst_ecl_trk,
				  int limit_mdst_klm_cluster_hit,
				  int* evtcls_flag )
{

  Belle_event_Manager& EvHeadMgr = Belle_event_Manager::get_manager();
  //  if( EvHeadMgr[0].ExpNo() < 39) return;
  
  // dout(Debugout::DDEBUG,"EventClassify") << Mdst_ecl_trk_Manager::get_manager().count() <<" "
  //	    << limit_mdst_ecl_trk << std::endl;
  if(Mdst_ecl_trk_Manager::get_manager().count() > limit_mdst_ecl_trk) {
    evtcls_flag[15]=1;
  }

  // dout(Debugout::DDEBUG,"EventClassify") << Mdst_klm_cluster_hit_Manager::get_manager().count() <<" "
  //	    << limit_mdst_klm_cluster_hit << std::endl;
  if(Mdst_klm_cluster_hit_Manager::get_manager().count() 
     > limit_mdst_klm_cluster_hit ) {
    evtcls_flag[15]+=2;
  }

  // dout(Debugout::DDEBUG,"EventClassify") <<"flag[15]=" << evtcls_flag[15] << std::endl;

}




#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
