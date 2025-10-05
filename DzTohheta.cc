//---------------------------------------//
//for e+e- -> ccbar -> D* X...           //
//       D*+ -> D0 pi_s+ +c.c.           //
//        chn1: D0 -> K- pi+ eta         //
//        chn2: D0 -> pi+ pi- eta        //
//        chn3: D0 -> K+ K- eta          //
//---------------------------------------//

#include "belle.h"
#include "DzTohheta.h"
#include "userinfo.h"
#include "dmixutil.h"
#include <stdlib.h>
#include "kfitter/kmakemother.h"
#include "benergy/BeamEnergy.h"


#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

const double CV = 2.99792458;
const double ANGLE = 22.0;//mrad
const double BF = 1.5; //T
const int IELEC = 0;
const int IMUON = 1;
const int IPION = 2;
const int IKAON = 3;
const int IPROT = 4;

using namespace std;

extern "C" Module_descr *mdcl_DzTohheta()
{
  DzTohheta *module = new DzTohheta;
  Module_descr *dscr = new Module_descr("DzTohheta",module);
  IpProfile::define_global(dscr);
  dscr->define_param ( "mc","mc or data", &module->mc2exp );
  dscr->define_param ( "YnsNo","YnS type", &module->YnsNo );
  dscr->define_param ( "SaveTruth","YnS type", &module->SaveTruth );
  dscr->define_param ( "myDEBUG", "debug flag", &module->myDEBUG );
  dscr->define_param ( "MFitFlag", "M fit flag", &module->MFitFlag );
  return dscr;
}

// ptypeName's defination ref. to ~/ref/decay.dec copied from QQ generator
DzTohheta::DzTohheta(void)
  : dstpsmin(2.3), mdifmax(0.175), md0winL(0.165), md0winR(0.155),
    metawinL(0.08), metawinR(0.06), dstpsmin_fit(2.4), mdifmax_fit(0.168), 
    md0win_fitL(0.155), md0win_fitR(0.145),
    petamin(0.28), Egdiff(1.00),
    kidmin(0.6), pidmax(0.6), pidmax2(0.9), 
    YnsNo(0), SaveTruth(0), myDEBUG(0), MFitFlag(1), mc2exp(1),
    m_ptypeKP("K+"), m_ptypePIP("PI+"),
    m_ptypePI0("PI0"), m_ptypeETA("ETA"),
    m_ptypeD0("D0"),  m_ptypeD0B("D0B"),
    m_ptypeDstarP("D*+"), m_ptypeDstarM("D*-"),
    m_kaonP(), m_kaonM(), m_pionP(), m_pionM(),
    m_pionP2(), m_pionM2(), m_pi0(), m_eta()
{
  m_mc = 1;
  dzcut = 3.0;//cm
  drcut = 1.0;//cm
}

void DzTohheta::hist_def(void)
{
  myDEBUG && std::cout << "Histogram defined" << std::endl;
  extern BelleTupleManager *BASF_Histogram;
  const char *tags =  //MC  
    "ndst_ch1 ndst_ch2 ndst_ch3 ynsno "
    "expno runno evtno chnd0 chndst "
    "pi0flag1 pi0flag2 pi0mmin1 pi0mmin2 "
    "pi0pmin1 pi0pmin2 pi0asye1 pi0asye2 "
    "massd0 massdif md00 mdif0 dstps cosdst "
    "mchip0 m0p0 pp0 ediffp0 "
    "nrk nzk nrp nzp drk dzk drp dzp "
    "kaonid pionid kpnid kmuid pmuid keid peid "
    "kp pp p0p kps pps p0ps kpxy ppxy p0pxy cosk cosp cosp0 "
    "cosgmb1 cosgmb2 mpi0 egm1 egm2 e9251 e9252 cosgm "
    "d0px d0py d0pz d0p d0pxy d0ps cosd0 "
    "m13m m23m m12m m13 m23 m12 m123 "
    "cos1 cos2 cos3 cosh13 cosh23 cosh12 " 
    "chisqv chisqb chisqs sumchisq "
    "spp sppxy spps nrsp nzsp drsp dzsp "
    "spkid speid spmuid cossp cosspd cosdedst "
    "leng lxy el elxy "
    "mc mcspml1 mcspmi1 mcspmn1 mcspmn1g "
    "mcspml2 mcspmi2 mcspmn2 mcspmn2g photo "
    "mckml1 mckmi1 mckmn1 mckmn1g "
    "mckml2 mckmi2 mckmn2 mckmn2g "
    "mckml3 mckmi3 mckmn3 mckmn3g "
    "mcpml1 mcpmi1 mcpmn1 mcpmn1g "
    "mcpml2 mcpmi2 mcpmn2 mcpmn2g "
    "mcpml3 mcpmi3 mcpmn3 mcpmn3g "
    "mcp0ml1 mcp0mi1 mcp0mn1 mcp0mn1g "
    "mcp0ml2 mcp0mi2 mcp0mn2 mcp0mn2g "
    "mcp0ml3 mcp0mi3 mcp0mn3 mcp0mn3g "
    "mcmkp mcmkp0 mcmpp0 mcspi mcspl gspi gspml "
    "mckl gkml mcpl gpml mcp0l gp0ml mcd01 "
    "mcdp mclen1 mclen2 "
    "mcg1ml1 mcg1mi1 mcg2ml1 mcg2mi1 " 
    "mcg1ml2 mcg1mi2 mcg2ml2 mcg2mi2 ";

  const char *tags1 =  //exp
    "ndst_ch1 ndst_ch2 ndst_ch3 ynsno "
    "expno runno evtno chnd0 chndst "
    "pi0flag1 pi0flag2 pi0mmin1 pi0mmin2 "
    "pi0pmin1 pi0pmin2 pi0asye1 pi0asye2 "
    "massd0 massdif md00 mdif0 dstps cosdst "
    "mchip0 m0p0 pp0 ediffp0 "    
    "nrk nzk nrp nzp drk dzk drp dzp "
    "kaonid pionid kpnid kmuid pmuid keid peid "
    "kp pp p0p kps pps p0ps kpxy ppxy p0pxy cosk cosp cosp0 "
    "cosgmb1 cosgmb2 mpi0 egm1 egm2 e9251 e9252 cosgm "
    "d0px d0py d0pz d0p d0pxy d0ps cosd0 "   
    "m13m m23m m12m m13 m23 m12 m123 "
    "cos1 cos2 cos3 cosh13 cosh23 cosh12 "    
    "chisqv chisqb chisqs sumchisq "
    "spp sppxy spps nrsp nzsp drsp dzsp "
    "spkid speid spmuid cossp cosspd cosdedst "
    "leng lxy el elxy "
    "mc ";

  //+-------------
  // save a histogram to check the selection cut flow
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  m_event_counter = BASF_Histogram->histogram("Event Counter",25,0,25);
  m_selection_counter = BASF_Histogram->histogram("Selection counter",20,0,20);
  m_bcs_counter1 = BASF_Histogram->histogram("BCS Counter in chn1", 15, 0, 15);  
  m_bcs_counter2 = BASF_Histogram->histogram("BCS Counter in chn2", 15, 0, 15);  
  m_bcs_counter3 = BASF_Histogram->histogram("BCS Counter in chn3", 15, 0, 15);  
  m_bcs_counter4 = BASF_Histogram->histogram("BCS Counter in chn4", 61, -0.5, 60.5 ); 

  if(mc2exp==0) { //exp. data
    DzNtupe = BASF_Histogram->ntuple("D0_RS_Info_after_BCS",tags1);
//  DzNtupe2 = BASF_Histogram->ntuple("D0_RS_Info_before_BCS",tags1);
  }else { //MC data
    DzNtupe = BASF_Histogram->ntuple("D0_RS_Info_after_BCS",tags);
//  DzNtupe2 = BASF_Histogram->ntuple("D0_RS_Info_before_BCS",tags);
    //+------------
    // ntuple defination for MC Truth info.
    m_d0_hep = BASF_Histogram->ntuple("D0_Event_Info_Hep",
       "expno chnd0 idhep p t moihep mops mocos gm13 gm23 gm12 cosh13 cosh23 cosh12 ");
  }

}

void DzTohheta::init(int *status)
{
  Ptype dummy("E-");
  std::cout << " " << "\n"
                         << "BASF init is called"
                   << "\n"
                         << "D* P cut low (before fit) = "                <<  dstpsmin  << "GeV/c"
                   << "\n"
                         << "D* mass difference cut high (before fit) = " <<  mdifmax   << "GeV"
                   << "\n"
                         << "D* p cut low (after fit) = "                 << dstpsmin   << "GeV/c"
                   << "\n"
                         << "D* mass difference cut high (after fit) = "  << mdifmax    << "GeV"
                   << "\n"
                   << " " << std::endl;
  return;
}

void DzTohheta::begin_run(BelleEvent *evptr, int *status)
{
  eid::init_data();

  //+-----------
  // get IP profile data
  //~~~~~~~~~~~~~~~~~~~~
  IpProfile::begin_run();
  IpProfile::usable();
  IpProfile::dump();

  Belle_runhead_Manager &runMgr = Belle_runhead_Manager::get_manager();

  if ( !runMgr.count() ) {
    std::cerr << "!!!!!Can not Access to Belle_RunHead!!!!!" << std::endl;
    std::cerr << "!!!!!Regard the data as MC sample!!!!!" << std::endl;
  }
  else if ( !runMgr[0].ExpMC() == 1) {
    m_mc = 0;
    std::cout << "              "               << "\n"
              << "Exp Analysis"                 << "\n"
              << "      "    <<std::endl;
  }
  else {
    std::cout << "              "               << "\n"
              << "MC Analysis"                  << "\n"
              << "      "    <<std::endl;
  }
  //+-----------
  //Beam Energy
  m_LER = BeamEnergy::E_LER();  // e+ beam Energy
  m_HER = BeamEnergy::E_HER();  // e- beam energy
  m_CROSS = BeamEnergy::Cross_angle()*1000.0;  // cross angle (mrad)
  if( fabs( m_LER )<1 || fabs( m_HER )<1 )
  {      
    m_CROSS = ANGLE;   // default Y(4S)
    m_HER = 7.998213;
    m_LER = 3.499218;
    if( YnsNo == 1 ) { //Y(1S)
      m_HER = 7.1514;
      m_LER = 3.1287;
    } else if( YnsNo == 2 ){ //Y(2S)
      m_HER = 7.577331;
      m_LER = 3.315082;
    }else if( YnsNo == 3 ){ //Y(3S)
      m_HER = 7.82741;
      m_LER = 3.42449;
    }else if( YnsNo == 4 ) {//Y(4S)_continuum
      m_HER = 7.952693;
      m_LER = 3.479303;
    }else if( YnsNo == 5 ) { //Y(5S)
      m_HER = 8.21504;
      m_LER = 3.59408;
    }
  }
  std::cout << "beam energy " << "\n"
            << "HER E = \t" << m_HER << "GeV" << "\n"
            << "LER E = \t" << m_LER << "GeV" << "\n"
            << "   "   << std::endl;

  return;
}

void DzTohheta::event(BelleEvent *evptr, int *status)
{
  myDEBUG && std::cout << "begin event!!!!!" << std::endl;
  *status = -1;
  m_event_counter->accumulate(0.,1.);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  m_expNo = 0; m_evtNo = 0;
  m_runNo = 0; 
  Belle_event_Manager &evtMgr = Belle_event_Manager::get_manager();
  if( evtMgr.count() !=0 ) {
    const int MASK28BIT = 0x0FFFFFFF;
    m_expNo = evtMgr[0].ExpNo();
    m_runNo = evtMgr[0].RunNo();
    m_evtNo = (int)(evtMgr[0].EvtNo() & MASK28BIT);
  }

  //~save MC info before rec.ed~~
  if( m_mc && SaveTruth ) {
    hepevtInfo(m_d0_hep);
  }
  
  //define particle list
  std::vector<Particle> m_eta0, m_pi00;
  std::vector<Particle> D0Cand, D0BCand;
  std::vector<Particle> DstCand0, DstCand_s, DstCand_bcs;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Mdst_charged_Manager &chgMgr = Mdst_charged_Manager::get_manager();
  Mdst_gamma_Manager &gamMgr = Mdst_gamma_Manager::get_manager();
  Evtcls_hadron_info_Manager &clsMgr = Evtcls_hadron_info_Manager::get_manager();

  //~~~~~~h^+/h^-/pi_s/X~~~~~~~~~~~~~~~~~~~~~~~~
  if( chgMgr.count() < 4) {goto END;}
  m_event_counter->accumulate(1., 1.);

  // event's some flag;
  m_r2 = 0;
  if( clsMgr.count() ) {
    m_r2 = clsMgr[0].R2();
    m_ntrk = clsMgr[0].Ntrk();
    m_ncls = clsMgr[0].Ncls();
    m_eivs = clsMgr[0].Evis();
  }

  makeKPi(m_kaonP, m_kaonM, m_pionP, m_pionM, 0);
  m_pionP2 = m_pionP; m_pionM2 = m_pionM;
  myDEBUG && std::cout << "Kaon && Pion size:" << m_kaonP.size() <<", "<< m_kaonM.size() << ", "
                                             << m_pionP.size() <<", "<< m_pionM.size() << std::endl;

  // K: R=L(K)/(L(K)+L(pi))>=pidmin
  withKaonId( m_kaonP, kidmin, 3,1,5,3,2 );
  withKaonId( m_kaonM, kidmin, 3,1,5,3,2 );
  // pi: R=L(K)/(L(K)+L(pi))<pidmax
  withPionId( m_pionP, pidmax, 3,1,5,3,2 );
  withPionId( m_pionM, pidmax, 3,1,5,3,2 );
  m_event_counter->accumulate(2.,1.);

  withSVD2( m_kaonP, 2, 2 );
  withSVD2( m_kaonM, 2, 2 );
  withSVD2( m_pionP, 2, 2 );
  withSVD2( m_pionM, 2, 2 );
  m_event_counter->accumulate(3.,1.);
  myDEBUG && std::cout << "After PID, (K Pi) size: " << m_kaonP.size() <<", "<< m_kaonM.size() << ", "
                                                   << m_pionP.size() <<", "<< m_pionM.size() << std::endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  withPionId( m_pionP2, pidmax2, 3,1,5,3,2 );
  withPionId( m_pionM2, pidmax2, 3,1,5,3,2 );
  m_event_counter->accumulate(4.,1.);

  withDrDzCut( m_pionP2, drcut, dzcut );
  withDrDzCut( m_pionM2, drcut, dzcut ); //no SVD cut for pionP2,pionM2 to avoid eff loss
  m_event_counter->accumulate(5.,1.);
  myDEBUG && std::cout << "slow pi size: " << m_pionP2.size() << "(+), " << m_pionM2.size() <<"(-)"<< std::endl;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(mc2exp) {
    setGenHepInfoF(m_kaonP); setGenHepInfoF(m_kaonM);
    setGenHepInfoF(m_pionP); setGenHepInfoF(m_pionM);
    setGenHepInfoF(m_pionP2); setGenHepInfoF(m_pionM2);
  }
  m_event_counter->accumulate(6.,1.);

  if( m_pionP2.size()+m_pionM2.size()<1 ) {goto END;}
  m_event_counter->accumulate(7.,1.);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~

  makegam(m_gamma);
  if( m_gamma.size()<2 ) {goto END;}
  myDEBUG && std::cout << "gamma list size:" << m_gamma.size() << std::endl;
   
  combination( m_eta0, m_ptypeETA, m_gamma, m_gamma, metawinL, metawinR );
  setUserInfo( m_eta0, 1 );
  if( m_eta0.size()<1 ) {goto END;}
  myDEBUG && std::cout << "eta0 size:\t" << m_eta0.size() << std::endl;
   
  m_event_counter->accumulate(8.,1.);

  //+-------------------------
  //      m_eta0 -> m_eta
  // with cut min momentum and max energy asymmetry of two gammas
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  selectPi0( m_eta, m_eta0, petamin, Egdiff, 100 );
  if( mc2exp ) {
    setGenHepInfoP2(m_eta);
  }
  myDEBUG && std::cout << "eta size:\t" << m_eta.size() << std::endl;
  if( m_eta.size()<1 ) {goto END;}
  m_event_counter->accumulate(9.,1.);

  //~~~~Self-conjugated D0 multi-channels~~~~~~~~~~~~~~~~
  combination(D0Cand, m_ptypeD0, m_kaonM, m_pionP, m_eta, md0winL, md0winR );
  setUserInfo(D0Cand, 1);
  combination(D0BCand, m_ptypeD0B, m_kaonP, m_pionM, m_eta, md0winL, md0winR );
  setUserInfo(D0BCand, -1);
  combination(D0Cand, m_ptypeD0, m_pionP, m_pionM, m_eta, md0winL, md0winR );
  setUserInfo(D0Cand, 2);
  combination(D0BCand, m_ptypeD0B, m_pionM, m_pionP, m_eta, md0winL, md0winR );
  setUserInfo(D0BCand, -2);
  combination(D0Cand, m_ptypeD0, m_kaonP, m_kaonM, m_eta, md0winL, md0winR );
  setUserInfo(D0Cand, 3);
  combination(D0BCand, m_ptypeD0B, m_kaonM, m_kaonP, m_eta, md0winL, md0winR );
  setUserInfo(D0BCand, -3);
  if( D0Cand.size()+D0BCand.size() == 0 ) {goto END;}
  myDEBUG && std::cout << "D0 size:" << D0Cand.size() <<" D0B size: "<< D0BCand.size() << std::endl;

  //~~~D*+ and D*-~~~~~~~~~~~~~~~~~~
  combination( DstCand0, m_ptypeDstarP, D0Cand, m_pionP2 );
  setUserInfo( DstCand0, 1 );
  combination( DstCand0, m_ptypeDstarM, D0BCand, m_pionM2 );
  setUserInfo( DstCand0, -1 );
  myDEBUG && std::cout << "DstCand0 size:" << DstCand0.size() << std::endl;
  if( DstCand0.size() == 0 ) {goto END;}

  withMassDifCut( DstCand0, 0.135, mdifmax, 0 );
  withPSCut(DstCand0, dstpsmin );
  myDEBUG &&  std::cout << "DstCand0 size with M&pScut: " << DstCand0.size()<<endl;
  m_event_counter->accumulate(10.,1.);

  if( DstCand0.size() == 0 ) {goto END;}
  myDEBUG && std::cout<<"fitD0Event for D0 vertex"<<endl;
  fitD0Event( DstCand0, 0);
  myDEBUG && std::cout<<"DstCand0 size after D0 vertexfit: "<<DstCand0.size()<<endl;
  m_event_counter->accumulate(11.,1.);

  //~~~~~~~~~~~~~~~~~~~~~~~~
  if( DstCand0.size() == 0 ) {goto END; }
  myDEBUG && std::cout << "addbeam2fit for D0 produced vertex" << endl;
  beamPosition( DstCand0, 0 );
  myDEBUG && std::cout << "DstCand0 size after addBeam2fit: " << DstCand0.size() << endl;
  m_event_counter->accumulate(12.,1.);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~
  withMassCut( DstCand0, m_ptypeD0.mass()-md0win_fitL, m_ptypeD0.mass()+md0win_fitR, 0 );
  m_event_counter->accumulate(14.,1.);

  myDEBUG && std::cout<<"pi_s refit to improve Q resolution"<<endl; 
  sprfmassDif( DstCand0, 0, 1 );
  m_event_counter->accumulate(13.,1.);


  //~~~~~~~~~~~~~~~~~~~~~~~~
  withPSCut2( DstCand0, dstpsmin_fit, m_HER, m_LER, 22.0 );
  m_event_counter->accumulate(15.,1.);
  withMassDifCut(DstCand0, 0.135, mdifmax_fit, 0 );
  m_event_counter->accumulate(16.,1.);

  //~~~~~~~~~~~~~~~~~~~~~~
  myDEBUG && std::cout << "DstCand0 size after p,M cut:" << DstCand0.size() << std::endl;
  if( DstCand0.size() == 0 ) {goto END;}
  m_event_counter->accumulate(17.,1.);
  if(mc2exp) {
    // D*->child(0):D0 setGenHepInfo
    setGenHepInfoR2(DstCand0, 0, 1);
    setGenHepInfoR(DstCand0);
  }

//Selection_Dst( DstCand_s, DstCand0, 0 );
//myDEBUG && std::cout<<"DstCand size before BCS: "<<DstCand_s.size()<<endl;
//if( DstCand_s.size() ==0 ) {goto END;}

  Selection_Dst( DstCand_s, DstCand0, 1 );
  myDEBUG && std::cout<<"DstCand size after optimal requirements (before BCS): "<<DstCand_s.size()<<endl;
  if( DstCand_s.size() ==0 ) {goto END;}

  // rec pi0(0.135 GeV) in [0.115, 0.155]GeV/c^2 and mchisq<50 
//combination( m_pi00, m_ptypePI0, m_gamma, m_gamma, 0.020,  0.020 );
  combination( m_pi00, m_ptypePI0, m_gamma, m_gamma, 0.015,  0.015 );
  setUserInfo( m_pi00, 1 );
  selectPi0( m_pi0, m_pi00, 0.0, 1.0, 50 );
  
  for( int i=0; i<DstCand_s.size(); i++ ) {
    int tmp_pi0vetoFlag1(0), tmp_pi0vetoFlag2(0); 
    double mpi0min1(0), mpi0min2(0);
    double ppi0min1(0), ppi0min2(0);
    double easypi01(-1.2), easypi02(-1.2);
    for( int j=0; j<m_pi0.size(); j++ ) {
      if( checkSame( DstCand_s[i].child(0).child(2).child(0), m_pi0[j] ) ) {
        tmp_pi0vetoFlag1 = 1; 
        if( fabs(m_pi0[j].mass()-0.135) < fabs(mpi0min1-0.135) ) {
          mpi0min1 = m_pi0[j].mass();  
          ppi0min1 = m_pi0[j].ptot();  
          easypi01 = (m_pi0[j].child(0).p().e()-m_pi0[j].child(1).p().e())/(m_pi0[j].child(0).p().e()+m_pi0[j].child(1).p().e());
        }
      }
      if( checkSame( DstCand_s[i].child(0).child(2).child(1), m_pi0[j] ) ) {
        tmp_pi0vetoFlag2 = 1; 
        if( fabs(m_pi0[j].mass()-0.135) < fabs(mpi0min2-0.135) ) {
          mpi0min2 = m_pi0[j].mass();  
          ppi0min2 = m_pi0[j].ptot();  
          easypi02 = (m_pi0[j].child(0).p().e()-m_pi0[j].child(1).p().e())/(m_pi0[j].child(0).p().e()+m_pi0[j].child(1).p().e());
        }
      }
    }
    dynamic_cast<DstUserInfo&>(DstCand_s[i].userInfo()).pi0veto1Flag( tmp_pi0vetoFlag1 ); 
    dynamic_cast<DstUserInfo&>(DstCand_s[i].userInfo()).pi0veto2Flag( tmp_pi0vetoFlag2 ); 
    dynamic_cast<DstUserInfo&>(DstCand_s[i].userInfo()).pi0veto1Mmin( mpi0min1 ); 
    dynamic_cast<DstUserInfo&>(DstCand_s[i].userInfo()).pi0veto2Mmin( mpi0min2 ); 
    dynamic_cast<DstUserInfo&>(DstCand_s[i].userInfo()).pi0veto1Pmin( ppi0min1 ); 
    dynamic_cast<DstUserInfo&>(DstCand_s[i].userInfo()).pi0veto2Pmin( ppi0min2 ); 
    dynamic_cast<DstUserInfo&>(DstCand_s[i].userInfo()).pi0veto1Easy( easypi01 ); 
    dynamic_cast<DstUserInfo&>(DstCand_s[i].userInfo()).pi0veto2Easy( easypi02 ); 
  }

  
//writeHist(DstCand_s, DzNtupe); 
//m_event_counter->accumulate(18.,1.);

  for( int i=0; i<DstCand_s.size(); i++ ) {
    double pi0flag1 = dynamic_cast<DstUserInfo&>(DstCand_s[i].userInfo()).m_pi0veto1Flag;  
    double pi0flag2 = dynamic_cast<DstUserInfo&>(DstCand_s[i].userInfo()).m_pi0veto2Flag;  
    double pi0mass1 = dynamic_cast<DstUserInfo&>(DstCand_s[i].userInfo()).m_pi0veto1Mmin ;
    double pi0mass2 = dynamic_cast<DstUserInfo&>(DstCand_s[i].userInfo()).m_pi0veto2Mmin ;
    // tight pi0-veto
//    if( (fabs(pi0mass1-0.135)<0.01&&pi0flag1==1)||(fabs(pi0mass2-0.135)<0.01&&pi0flag2==1) ) 
    // loose pi0-veto
    if( (fabs(pi0mass1-0.135)<0.01&&pi0flag1==1)&&(fabs(pi0mass2-0.135)<0.01&&pi0flag2==1) ) 
    { 
      DstCand_s.erase(DstCand_s.begin()+i); --i; 
    }
  }
   
  //to test BCS efficiency
//writeHist(DstCand_s, DzNtupe2 );
  
  bcselection( DstCand_bcs, DstCand_s);
  myDEBUG && std::cout<<"DstCand size after BCS: "<<DstCand_bcs.size()<<endl;
  if( DstCand_bcs.size()== 0 ) {goto END;}
  m_event_counter->accumulate(19.,1.);
  writeHist(DstCand_bcs, DzNtupe );
  
  m_event_counter->accumulate(20.,1.);
  

END:
  endEvent();
  eraseVector(m_eta0);
  eraseVector(m_pi00);
  eraseVector(D0Cand);
  eraseVector(D0BCand);
  eraseVector(DstCand0);
  eraseVector(DstCand_s);
  eraseVector(DstCand_bcs);

}
  
//==================================================//
// belows are sub-functions                         //
// will be called in above main function            //
// or other in userinfo.cc or dmixutil.cc           //
//==================================================//

//+-------------------------------------------------------
// save MC Truth info. 
// only need for signal MC to obtain the efficiency plane
//+-------------------------------------------------------
void DzTohheta::hepevtInfo(BelleTuple *hist)
{
  Gen_hepevt_Manager &genMgr = Gen_hepevt_Manager::get_manager();
  for(std::vector<Gen_hepevt>::iterator j=genMgr.begin();j!=genMgr.end();++j){
    if(abs(j->idhep()) == 421 &&
       j->daFirst() != 0 &&
       j->daLast() != 0){
      std::set<int,std::less<int> > finaldaID;
      fillChildId1eta(j->get_ID(), finaldaID);
      myDEBUG && std::cout << " D0 daughter number "
                         << finaldaID.size() << std::endl;
        Gen_hepevt gd0 = *j;
        HepLorentzVector mcd0p = HepLorentzVector(gd0.PX(),gd0.PY(),gd0.PZ(),gd0.E());
        HepPoint3D mcprdvtx = HepPoint3D(gd0.VX(),gd0.VY(),gd0.VZ());
        Hep3Vector mcdP = mcd0p.vect();
      if(finaldaID.size()==3) {
        std::vector<Gen_hepevt> finalda;
        for(std::set<int,less<int> >::iterator i=finaldaID.begin();
            i!=finaldaID.end();++i) {
          finalda.push_back(genMgr[(*i)-1]);
        }
        //+-----
        // MC truth decay channel 1 ( D0:421 -> pi:321 pi:211 pi0:111 )
        // MC truth decay channel 2 ( D0:421 -> pi:211 pi:211 eta:221 )
        // MC truth decay channel 3 ( D0:421 -> K:321  K:321  eta:221 )
        // the order is kept same as generator decay table
        //+-----
        for(int k=0; k<3; ++k) {
          int check(0);
          int idhepj=j->idhep();
          if( abs(finalda[k%3].idhep())==321 &&
              abs(finalda[(k+1)%3].idhep())==211 &&
              abs(finalda[(k+2)%3].idhep())==221 ) {
            check = 1; // D0 -> K- pi+ eta
          } else if( abs(finalda[k%3].idhep())==211 &&     
              abs(finalda[(k+1)%3].idhep())==211 && 
              abs(finalda[(k+2)%3].idhep())==221 ) {
            check = 2; // D0 -> pi+ pi- eta
          } else if( abs(finalda[k%3].idhep())==321 &&            
              abs(finalda[(k+1)%3].idhep())==321 && 
              abs(finalda[(k+2)%3].idhep())==221 ) {
            check = 3; // D0 -> K+ K- eta
          } else{check = 0;}

          if( check != 0 ) {
            int idhepj=j->idhep();
            hist->column("expno",  m_expNo);  // Exp #
            hist->column("chnd0", check);
            hist->column("idhep", (float) idhepj );  // D0 id_hep
            hist->column("p", sqrt(j->PX()*j->PX()+j->PY()*j->PY()+j->PZ()*j->PZ()));
            if(j->mother()) {
               Gen_hepevt *mother;
               mother =& j->mother();
               HepLorentzVector moP(mother->PX(),mother->PY(),mother->PZ(),
                                    mother->E());
    				HepLorentzVector moP2(mother->PX(),mother->PY(),mother->PZ(),
                                     mother->E());
               hist->column("moihep", (float) mother->idhep());
               hist->column("mops", pStar(moP,m_HER,m_LER,22.0).vect().mag());
    				double cosTheta = pStar(moP2, m_HER, m_LER, 22.0).vect().unit().dot( Hep3Vector(0,0,1) );
               hist->column("mocos", cosTheta ); 

            }
            int idhep1=finalda[(k)%3].idhep();
            int idhep2=finalda[(k+1)%3].idhep();
            HepLorentzVector genK1P(finalda[k%3].PX(), finalda[k%3].PY(),//D0.child(0)
                                    finalda[k%3].PZ(), finalda[k%3].E());
            HepLorentzVector genK2P(finalda[(k+1)%3].PX(),//D0.child(1)
                                    finalda[(k+1)%3].PY(),
                                    finalda[(k+1)%3].PZ(),
                                    finalda[(k+1)%3].E());
            HepLorentzVector genPi0P(finalda[(k+2)%3].PX(),//D0.child(2)
                                     finalda[(k+2)%3].PY(),
                                     finalda[(k+2)%3].PZ(),
                                     finalda[(k+2)%3].E());
            HepLorentzVector genD0P( j->PX(), j->PY(), j->PZ(), j->E() );  

            double  mK1Pi0 = (genK1P+genPi0P).mag();
            double  mK2Pi0 = (genK2P+genPi0P).mag();
            // K+Pi0 Vs K-Pi0 
            double  mKpPi0  = (genK1P+genPi0P).mag();
            double  mKmPi0 = (genK2P+genPi0P).mag();
            double  mKpKm = (genK1P+genK2P).mag();

            hist->column("gm13", (float)mK1Pi0 );
            hist->column("gm23", (float)mK2Pi0 );
            hist->column("gm12", (float)mKpKm );

            HepLorentzVector genKpP = genK1P ; 
            HepLorentzVector genKmP = genK2P ;
            HepLorentzVector genKpPi0P = genKpP + genPi0P ;
            HepLorentzVector genKpPi0P2 = genKpP + genPi0P ;
            HepLorentzVector genKmP2 = genKmP ;
            HepLorentzVector genKmPi0P = genKmP + genPi0P ;
            HepLorentzVector genKmPi0P2 = genKmP + genPi0P ;
            HepLorentzVector genKpP2 = genKpP ;
            HepLorentzVector genKpKmP = genKpP + genKmP ;
            HepLorentzVector genKpKmP2 = genKpP + genKmP ;
            HepLorentzVector genKpP3 = genKpP ;
            HepLorentzVector genPi0P2 = genPi0P ;
            genKpP.boost( -genKpPi0P.boostVector() );
            //add two lines below on Sep 12
            genKmP2.boost( -genD0P.boostVector() );
            genKpPi0P2.boost( -genD0P.boostVector() );
            genKmP2.boost( -genKpPi0P2.boostVector() );
            hist->column("cosh13", genKpP.vect().unit().dot(genKmP2.vect().unit()));
            genKmP.boost( -genKmPi0P.boostVector() );
            //add two lines below on Sep 12
            genKpP2.boost( -genD0P.boostVector() );
            genKmPi0P2.boost( -genD0P.boostVector() );
            genKpP2.boost( -genKmPi0P2.boostVector() );
            hist->column("cosh23", genKmP.vect().unit().dot(genKpP2.vect().unit()));
            genKpP3.boost( -genKpKmP.boostVector() );
            //add two lines below on Sep 12
            genPi0P2.boost( -genD0P.boostVector() );
            genKpKmP2.boost( -genD0P.boostVector() );
            genPi0P2.boost( -genKpKmP2.boostVector() );
            hist->column("cosh12", genKpP3.vect().unit().dot(genPi0P2.vect().unit()));
           
            hist->dumpData();
            break;
          }
        }
      }
    }
  }
} 

//+-------------------------------------------------------
// select gamma candidates
// most we direct use mdst_gamma_manager 
// without this matrix calculated
//+-------------------------------------------------------
void DzTohheta::makegam(std::vector<Particle> & gammalist)
{
  myDEBUG && std::cout << "Begin make gamma "<<std::endl;
  Ptype ptype_gamma("GAMM");

  Mdst_gamma_Manager& GamMgr  = Mdst_gamma_Manager::get_manager();
  Mdst_ecl_aux_Manager & eclaux = Mdst_ecl_aux_Manager::get_manager(); //GM

  for(Mdst_gamma_Manager::iterator g = GamMgr.begin();
      g!=GamMgr.end(); ++g){
    Mdst_gamma& gam = *g;
    Mdst_ecl& ecl = gam.ecl();
    Particle TmpGam(*g);
    double e9ovr = eclaux[ecl.get_ID()-1].e9oe25();
    if( e9ovr < 0.8 ) continue;
    if( TmpGam.e() < 0.03 || TmpGam.e() > 10.0 ) continue;
    // Here we check the unassociated track and their quality
    // match: the cluster matched with charged tracks in CDC
    //        =0: not match; -1:shower-trk; -2: charged-trk;
    // quality: ECL data quality. in gsim, =0:good cluster.
    //if(ecl.match() == 0 || ecl.quality() == 0)
    if(ecl.match() == 0 && ecl.quality() == 0)
    {
      double E = ecl.energy(), Theta = ecl.theta(), Phi=ecl.phi();
      HepLorentzVector Pgam(E*sin(Theta)*cos(Phi),
                            E*sin(Theta)*sin(Phi),
                            E*cos(Theta), E);
     // Symmetric Matrix
      HepSymMatrix errEcl(3,0); //error in ecl matrix on E, theta, phi


      errEcl[0][0] = ecl.error(0); //Energy
      errEcl[1][0] = ecl.error(1);
      errEcl[1][1] = ecl.error(2); //Phi
      errEcl[2][0] = ecl.error(3);
      errEcl[2][1] = ecl.error(4); //Theta
      errEcl[2][2] = ecl.error(5);

      HepSymMatrix errCart(4,0);
      HepMatrix jacobian(4,3,0);

      jacobian[0][0] =         cos(Phi)*sin(Theta);
      jacobian[0][1] =      -E*sin(Phi)*sin(Theta);
      jacobian[0][2] =       E*cos(Phi)*cos(Theta);
      jacobian[1][0] =         sin(Phi)*sin(Theta);
      jacobian[1][1] =       E*cos(Phi)*sin(Theta);
      jacobian[1][2] =       E*sin(Phi)*cos(Theta);
      jacobian[2][0] =                  cos(Theta);
      jacobian[2][1] =                         0.0;
      jacobian[2][2] =               -E*sin(Theta);
      jacobian[3][0] =                         1.0;
      jacobian[3][1] =                         0.0;
      jacobian[3][2] =                         0.0;

      errCart = errEcl.similarity(jacobian);
      HepSymMatrix dX(3,1); Hep3Vector X;
      dX *= 10.0;                               // position and error matrix   
      TmpGam.name("Gamma");
      TmpGam.momentum().momentum(Pgam,errCart);  //set momenta
      TmpGam.momentum().position(X,dX);
      TmpGam.userInfo(UserInfo());
//    dynamic_cast<UserInfo&>(TmpGam.userInfo()).e9ovr(e9ovr);
      gammalist.push_back(TmpGam);
      //+------------------------
      // setGenHepInfo for gamma list 
      // functions see /src/anal/particle/functions/utility.cc
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      setGenHepInfoG(gammalist);
    } // quality
  } //g = GamMgr.begin()
}

//+-------------------------------------------------------
// select eta candidates
// with min momentum and max energy asymmetry
//+-------------------------------------------------------
void DzTohheta::selectPi0(std::vector<Particle> & pi0list_s, std::vector<Particle> & pi0list, double ppi0min, double ediff, double mchisq=100 )
{
  myDEBUG && std::cout<<"Begin select pi0/eta list"<<endl;
  if( pi0list.size() == 0 ) return;
  myDEBUG && std::cout<< pi0list.size() <<endl; 
  for(int i=0; i< pi0list.size(); ++i ) {
    Particle p = pi0list[i] ;
    if( p.nChildren() != 2 ) {
      myDEBUG && std::cout<<"select pi0 size: "<<p.nChildren()<<endl;
      continue;
    }
    if( !p.child(0).mdstGamma() || !p.child(1).mdstGamma() ) {
      myDEBUG && std::cout<<"select pi0 two child no exist"<<endl;
      continue;
    }
    for(unsigned i=0; i<2; i++) setGammaError1(p.child(i));
    kmassfitter kmf;
    kmf.invariantMass(p.pType().mass());
    for(unsigned i=0; i<2; ++i )  {
      addTrack2fit(kmf,p.child(i));
    }
    int err0 = kmf.fit();
    if( err0 != 0 || kmf.chisq()>mchisq || kmf.chisq()<0. ) {
      myDEBUG && std::cout<<"pi0 massfitter failed"<<endl;
      continue;
    }
    double mmpi0 = (kmf.momentum(0)+kmf.momentum(1)).mag();
    dynamic_cast<UserInfo&>(p.userInfo()).mchisq(kmf.chisq());
    dynamic_cast<UserInfo&>(p.userInfo()).mass0(p.mass());
    dynamic_cast<UserInfo&>(p.userInfo()).mmass(mmpi0);

    double pi0p = p.p().rho();
    Particle pi0g1 = p.child(0);
    Particle pi0g2 = p.child(1);
    double pi0g1e = pi0g1.e();
    double pi0g2e = pi0g2.e();
    double pi0ghe = pi0g1e>=pi0g2e ? pi0g1e: pi0g2e;
    double pi0gle = pi0g1e<pi0g2e ? pi0g1e: pi0g2e;
    dynamic_cast<UserInfo&>(p.userInfo()).ghe(pi0ghe);
    dynamic_cast<UserInfo&>(p.userInfo()).gle(pi0gle);
    double gEdiff = abs( pi0g1e - pi0g2e )/(pi0g1e + pi0g2e);
    dynamic_cast<UserInfo&>(p.userInfo()).ediff(gEdiff);
    dynamic_cast<UserInfo&>(p.userInfo()).moment(pi0p);
    if( pi0p > ppi0min && gEdiff< ediff ) {
      pi0list_s.push_back( p );
    }
  } //loop
  myDEBUG && std::cout<<"End of pi0/eta list selection"<<endl;
}

//+-------------------------------------------------------
// (1)seperate the physical process: bbbar/ccbar/others 
//    according to particle idhep defination character
// (2)label the event kind:signal/rnd bg/cmb bg/et al
//+-------------------------------------------------------
int DzTohheta::GetMcType(Particle &Dst){
  Gen_hepevt_Manager &genMgr = Gen_hepevt_Manager::get_manager();
  int type = 0;
  for(std::vector<Gen_hepevt>::iterator it=genMgr.begin(); it!=genMgr.end();
        ++it){
    if (type==0 && (it->idhep()%1000)/100==4 || it->idhep()/1000==4)
      {  type = 1; } //ccbar
    if ((it->idhep()%1000)/100==5 || it->idhep()/1000==5)
      { type = 2; break; } //bbbar
  }
  Particle & D0 = Dst.relation().child(0);

  int mc=0, mc_sp=0, mc_d0=0, mc_dst=0;
  if( Dst.child(0).genHepevt() ) mc_d0=1;
  if( Dst.child(1).genHepevt() ) mc_sp=2;
  if( Dst.genHepevt() ) mc_dst=4;
  // all(D*/D0/pi_s) exist, mc=7
  mc = mc_d0 + mc_sp + mc_dst;
  return 10*type + mc;
}

//+------------------------------------------------
// D* candidates selection with cut
//+------------------------------------------------
void DzTohheta::Selection_Dst(std::vector<Particle> &dst_s,std::vector<Particle> &dst, int FinalCut )
{
  myDEBUG && std::cout<<"begin Dst selection"<<std::endl;
  if(dst.size() == 0) return;
//  cout<<"Dst list ...."<<endl;             
  for(unsigned i=0; i<dst.size(); ++i) {
    Particle &Dst = dst[i];
    Particle &D0 = Dst.child(0);
    Particle &Kaonp = D0.child(0); //
    Particle &Kaonm = D0.child(1); // 
    Particle &Pion0 = D0.child(2); // eta
    Particle &Slowpi = Dst.child(1);
    int mctype;
    if (mc2exp==1) mctype = GetMcType(Dst);
    m_selection_counter->accumulate(0.5,1.);
    int chnd0i = dynamic_cast<D0UserInfo&>(D0.userInfo()).m_channel; 

    //+--------
    // M,Q window
    double massd0 = dynamic_cast<D0UserInfo&>(D0.userInfo()).mass();
    double massdif = dynamic_cast<DstUserInfo&>(Dst.userInfo()).massDif();
    if ( massd0>1.94 || massd0<1.78 ) continue;  //1.78< M <1.94
    if ( massdif>0.15457 || massdif<0.13957 ) continue;  //0<Q<15
    if( FinalCut && !MFitFlag ) {
      if( SaveTruth==1 ) { //efficiency after applying M-calibration 
        if( abs(chnd0i)==1 && ((massd0-1.86483)*1.18+1.86483<1.842 || (massd0-1.86483)*1.18+1.86483>1.882 ) ) continue;  
        if( abs(chnd0i)==2 && ((massd0-1.86483)*1.18+1.86483<1.840 || (massd0-1.86483)*1.18+1.86483>1.884 ) ) continue;  
        if( abs(chnd0i)==3 && ((massd0-1.86483)*1.18+1.86483<1.850 || (massd0-1.86483)*1.18+1.86483>1.878 ) ) continue;  
      } else {
        if( abs(chnd0i)==1 && (massd0<1.842 || massd0>1.882 ) ) continue;  
        if( abs(chnd0i)==2 && (massd0<1.840 || massd0>1.884 ) ) continue;  
        if( abs(chnd0i)==3 && (massd0<1.850 || massd0>1.878 ) ) continue;  
      }
    }
    if( FinalCut && MFitFlag ) { 
      //efficiency determination after applying calibration on Q for MC
      if( mc2exp==1 && SaveTruth ) {
         if( fabs((massdif*1000.0-139.57-5.85)*1.30+5.85-5.86)>0.8 ) continue; 
      } else {
         if( fabs(massdif*1000.0-139.57-5.86)>0.8 ) continue;  //|Q-5.86|<0.8
      }  
    }
    m_selection_counter->accumulate(1.5,1.);
    //+--------
    //Dst cut
    HepLorentzVector tmpDstP4 = Dst.p(); 
    double  dstps=pStar(tmpDstP4, m_HER, m_LER, 22.0).vect().mag() ;
    if ( dstps<2.45 ) continue;
    if( FinalCut ) if ( dstps<2.70 ) continue;
    m_selection_counter->accumulate(2.5,1.);

    //+--------
    myDEBUG && std::cout<<"2.5 selectDst with D0 lifetime cut"<<endl;
    HepPoint3D dvtx = D0.momentum().decayVertex();
    HepPoint3D bvtx = D0.momentum().vertex();
    HepLorentzVector d0p = D0.momentum().p();

    if( FinalCut && abs(chnd0i)==2 ) {
      HepLorentzVector mKs_p4 = Kaonp.p() + Kaonm.p(); 
      double massKs = mKs_p4.mag(); 
      if( fabs(massKs - 0.497611)<0.010 ) continue;
    }

    //+------------------------
    myDEBUG && std::cout<<"3.5 selectDst with PID info."<<endl;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double kidp = kaonId2(Kaonp);
    double kidm = kaonId2(Kaonm);
    double spid = kaonId2(Slowpi);

    double kmuidp = muonId(Kaonp);
    double kmuidm = muonId(Kaonm);
    double spmuid = muonId(Slowpi);

    double keidp = eId(Kaonp);
    double keidm = eId(Kaonm);
    double speid = eId(Slowpi);

    if( abs(chnd0i)==1 ) { // D0 -> K- pi+ eta
      if( kidp<0.6 || kidm>0.6 ) continue;
      m_selection_counter->accumulate(3.5,1.);
    } else if( abs(chnd0i)==2 ) { // D0 -> pi+ pi- eta
      if( kidp>0.6 || kidm>0.6 ) continue;
      m_selection_counter->accumulate(4.5,1.);
    } else if ( abs(chnd0i)==3 ) { // D0 -> K+ K- eta 
      if (kidp<0.6 || kidm<0.6) continue;
      m_selection_counter->accumulate(5.5,1.);
    }
    
    //slow pion
    if( spid>0.6 ) continue;
    m_selection_counter->accumulate(6.5,1.);

    // kaon,pi via mu electron pid cut
    if ( kmuidp>0.95 || kmuidm>0.95 || spmuid>0.95 ) continue;
    if ( keidp>0.95 || keidm>0.95 || speid>0.95 ) continue;
    m_selection_counter->accumulate(7.5,1.);

    //+----------------
    myDEBUG && std::cout<<"7.5 selectDst with gamma cut"<<endl;
    Mdst_gamma gm0 = Pion0.child(0).relation().mdstGamma() ;
    Mdst_gamma gm1 = Pion0.child(1).relation().mdstGamma() ;
    HepLorentzVector gm0p4 = Pion0.child(0).p();
    HepLorentzVector gm1p4 = Pion0.child(1).p();
    Hep3Vector gm0p3(gm0p4.px(), gm0p4.py(), gm0p4.pz() );
    Hep3Vector gm1p3(gm1p4.px(), gm1p4.py(), gm1p4.pz() );
    double cosgmb1=gm0p4.vect().unit().dot( Hep3Vector(0,0,1));
    double cosgmb2=gm1p4.vect().unit().dot( Hep3Vector(0,0,1));

    double egm1= gm0p3.mag() ;
    double egm2= gm1p3.mag() ;
    Mdst_ecl_aux_Manager& eclaux_mgr = Mdst_ecl_aux_Manager::get_manager();
    double e9251= eclaux_mgr( Panther_ID(gm0.ecl().get_ID()) ).e9oe25() ;
    double e9252= eclaux_mgr( Panther_ID(gm1.ecl().get_ID()) ).e9oe25() ;

    // gamma energy cut for barrel/endcap
    if(cosgmb1>-0.61566 && cosgmb1<0.83867)
      if(egm1<0.05) continue;
    if(cosgmb2>-0.61566 && cosgmb2<0.83867)
      if(egm2<0.05) continue;//barrel

    if(cosgmb1<-0.61566 || cosgmb1>0.83867)
      if(egm1<0.10) continue;
    if(cosgmb2<-0.61566 || cosgmb2>0.83867)
      if(egm2<0.10) continue;//endcap
    m_selection_counter->accumulate(8.5,1.);
    // e925 cut
    if (e9251<0.8 || e9252<0.8) continue;
    m_selection_counter->accumulate(9.5,1.);

    double pi0p = Pion0.ptot();
    double mchipi0 = dynamic_cast<UserInfo&>(Pion0.userInfo()).mchisq() ;                                 
    if ( pi0p<0.30 ) continue;
    if( FinalCut ) { 
      if( pi0p<0.70 ) continue;
      if( mchipi0>8 ) continue;
    }
   
    m_selection_counter->accumulate(10.5, 1.);

    HepLorentzVector gm0p4Tmp = Pion0.child(0).p(); 
    gm0p4Tmp.boost( -Pion0.p().boostVector() );
    double cosHelEta = gm0p4Tmp.vect().unit().dot( Pion0.p().vect().unit() );
    if( FinalCut ) {  
       if( fabs( cosHelEta )>0.85 ) continue; 
    }

    //vertex chisq    
    double chisqv=dynamic_cast<D0UserInfo&>(D0.userInfo()).chisqV();
    double chisqb=dynamic_cast<D0UserInfo&>(D0.userInfo()).chisqB() ;
    double chisqs=dynamic_cast<DstUserInfo&>(Dst.userInfo()).chisqS();

    double ndfv  =dynamic_cast<D0UserInfo&>(D0.userInfo()).ndfV();
    double ndfb  =dynamic_cast<D0UserInfo&>(D0.userInfo()).ndfB();
    double ndfs  =dynamic_cast<DstUserInfo&>(Dst.userInfo()).ndfS();

    if ((chisqv+chisqb+chisqs)/(ndfv+ndfb+ndfs)>100) continue;
    if( FinalCut ) 
       if( (chisqv+chisqb+chisqs)>50 ) continue;
    m_selection_counter->accumulate(11.5,1.);

    dst_s.push_back(Dst);
  }
  myDEBUG && std::cout<<"End select D* candidates: "<< dst_s.size() <<endl;
}

//+------------------------------------------------
// best candidate selection (BCS) for multi-candidates
//+------------------------------------------------
void
DzTohheta::bcselection(std::vector<Particle> &dstr_bcs,
                       std::vector<Particle> &dstr_s  )
{
  myDEBUG && std::cout<<"Begin Best Candidate Selection"<<endl;
  if( dstr_s.size()<1 ) return;
  std::vector<Particle> dstch1, dstch2, dstch3, dstch4 ;
  std::vector<Particle> dstch1f, dstch2f, dstch3f, dstch4f;

  for( int i=0; i<dstr_s.size(); ++i ) {
    int chnd0i = dynamic_cast<D0UserInfo&>(dstr_s[i].child(0).userInfo()).m_channel;
    myDEBUG && std::cout<<"D0 channel No. is "<< chnd0i <<std::endl;
    if( abs(chnd0i)==1 ) 	  dstch1.push_back(dstr_s[i]);
    else if( abs(chnd0i)==2 ) dstch2.push_back(dstr_s[i]);
    else if( abs(chnd0i)==3 ) dstch3.push_back(dstr_s[i]);
    else std::cout<<"RS: not our analysis channel"<<endl;
  }
  myDEBUG && std::cout<< "Size in chn1: "<< dstch1f.size() <<"\tchn2: "<< dstch2f.size()
					<<"\tchn3: "<< dstch3f.size() <<"\tchn4: "<< dstch4f.size() <<std::endl;
  myDEBUG && std::cout<<"all chn BCS, go go go----"<<std::endl;

  if( dstch1.size()>0 ) {
    if( dstch1.size()< 20 ) m_bcs_counter4->accumulate( dstch1.size(), 1 );  
    bcselection0(dstch1f, dstch1, m_bcs_counter1);
  }
  if( dstch2.size()>0 ) {
    if( dstch2.size()<20 ) m_bcs_counter4->accumulate( dstch2.size()+20, 1 );  
    bcselection0(dstch2f, dstch2, m_bcs_counter2);
  }
  if( dstch3.size()>0 ) { 
    if( dstch3.size()<20 ) m_bcs_counter4->accumulate( dstch3.size()+40, 1 );  
    bcselection0(dstch3f, dstch3, m_bcs_counter3);
  }

  if( dstch1f.size()>1 || dstch2f.size()>1 || dstch3f.size()>1 ) {
    std::cout<<"Best candidate selection failed"<<endl;
    return ;
  } 
  
  if( dstch1f.size()>0 ) m_bcs_counter4->accumulate( 0, 1 );  
  if( dstch2f.size()>0 ) m_bcs_counter4->accumulate( 20, 1 );  
  if( dstch3f.size()>0 ) m_bcs_counter4->accumulate( 40, 1 );  

  myDEBUG && std::cout<< "Size in chn1: "<< dstch1f.size() <<"\tchn2: "<< dstch2f.size()
					<<"\tchn3: "<< dstch3f.size() <<"\tchn4: "<< dstch4f.size() <<std::endl;
  for(int i=0; i<dstch1f.size(); ++i) dstr_bcs.push_back( dstch1f[i] );
  for(int i=0; i<dstch2f.size(); ++i) dstr_bcs.push_back( dstch2f[i] );
  for(int i=0; i<dstch3f.size(); ++i) dstr_bcs.push_back( dstch3f[i] );
  if( dstr_bcs.size()>3 ) {
    std::cout<<"Best candidate selection failed"<<std::endl;
    return ;
  } 
  myDEBUG && std::cout<<"End BCS: DstCand.size = "<< dstr_bcs.size() <<std::endl;
}

void
DzTohheta::bcselection0(std::vector<Particle> &dstr_bcs,
                        std::vector<Particle> &dstr_s, BelleHistogram *m_bcs_counter)
{
  myDEBUG && std::cout<<"Begin Each channel BCS"<<endl;
  Particle Dstr_best;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(dstr_s.size()> 0) Dstr_best = dstr_s[0];
  else return;

  if(dstr_s.size()> 1){
    m_bcs_counter->accumulate(14.5, 1.); 
    for(unsigned i=1; i<dstr_s.size(); ++i) {
      Particle &Dst = dstr_s[i];
      Particle &D0 = Dst.child(0);
      Particle &Kaon = D0.child(0);
      Particle &Pion = D0.child(1);
      Particle &Pion0 = D0.child(2); // eta 
      Particle &Slowpi = Dst.child(1);
      Mdst_gamma gm0 = Pion0.child(0).relation().mdstGamma() ;
      Mdst_gamma gm1 = Pion0.child(1).relation().mdstGamma() ;

      Particle &D0_b = Dstr_best.child(0);
      Particle &Slowpi_b = Dstr_best.child(1);
      Particle &Kaon_b = D0_b.child(0);
      Particle &Pion_b = D0_b.child(1);
      Particle &Pion0_b = D0_b.child(2);
      Mdst_gamma gm0_b = Pion0_b.child(0).relation().mdstGamma();
      Mdst_gamma gm1_b = Pion0_b.child(1).relation().mdstGamma();

      if (bcscheckSame(D0_b,D0)){  //1. same d0
        m_bcs_counter->accumulate(0.5,1.);
      } else if (bcscheckSame(Slowpi_b,Slowpi)){ //2. same pi_s
        m_bcs_counter->accumulate(1.5,1.);
        if (bcscheckSame(Pion0_b,Pion0)){ //a. same pi0
          m_bcs_counter->accumulate(2.5,1.);
        } else { //b. dif pi0
          m_bcs_counter->accumulate(3.5,1.);
        }
      }else { //3. diff D0 and diff pi_s
        m_bcs_counter->accumulate(4.5,1.);
        if (bcscheckSame(Pion0_b,Pion0)){ //a. same pi0 
          m_bcs_counter->accumulate(5.5,1.);
          if( bcscheckSame(Kaon_b,Kaon) && bcscheckSame(Slowpi_b,Pion) && 
				  bcscheckSame(Pion_b,Slowpi) ) { //same pi0 and same k( swap pi_s and pi)
            m_bcs_counter->accumulate(6.5,1.);
          } else if( bcscheckSame(Pion_b,Pion) && bcscheckSame(Slowpi_b,Pion) && 
				         bcscheckSame(Kaon_b,Slowpi) ) { //same pi0 and same pi( swap pi_s and kaon)
            m_bcs_counter->accumulate(7.5,1.);
          }else { // not (same pi0 and same K and swap(pi_s,pi))
            m_bcs_counter->accumulate(8.5,1.);
          }
        } else { //b. dif pi0 
          m_bcs_counter->accumulate(9.5,1.);
        }
      }
     
      //double chisq = ( dynamic_cast<DstUserInfo&>(Dst.userInfo()).chisqS()+dynamic_cast<D0UserInfo&>(D0.userInfo()).chisqV()+dynamic_cast<D0UserInfo&>(D0.userInfo()).chisqB() ) / ( dynamic_cast<DstUserInfo&>(Dst.userInfo()).ndfS()+dynamic_cast<D0UserInfo&>(D0.userInfo()).ndfV()+dynamic_cast<D0UserInfo&>(D0.userInfo()).ndfB() );
      //double chisq_b = ( dynamic_cast<DstUserInfo&>(Dstr_best.userInfo()).chisqS()+dynamic_cast<D0UserInfo&>(D0_b.userInfo()).chisqV()+dynamic_cast<D0UserInfo&>(D0_b.userInfo()).chisqB() ) / ( dynamic_cast<DstUserInfo&>(Dstr_best.userInfo()).ndfS()+dynamic_cast<D0UserInfo&>(D0_b.userInfo()).ndfV()+dynamic_cast<D0UserInfo&>(D0_b.userInfo()).ndfB() );
      double chisq = dynamic_cast<DstUserInfo&>(Dst.userInfo()).chisqS()+dynamic_cast<D0UserInfo&>(D0.userInfo()).chisqV()+dynamic_cast<D0UserInfo&>(D0.userInfo()).chisqB() ; 
      double chisq_b = dynamic_cast<DstUserInfo&>(Dstr_best.userInfo()).chisqS()+dynamic_cast<D0UserInfo&>(D0_b.userInfo()).chisqV()+dynamic_cast<D0UserInfo&>(D0_b.userInfo()).chisqB();

      double mchipi0 = dynamic_cast<UserInfo&>(Pion0.userInfo()).mchisq() ;                                 
      double mchipi0_b = dynamic_cast<UserInfo&>(Pion0_b.userInfo()).mchisq() ;                                 
 
      if ( chisq_b+mchipi0_b > chisq+mchipi0 )  {
        Dstr_best = Dst;
      }
    }//end dstr_s
  } 
  if (dstr_s.size()>0) dstr_bcs.push_back(Dstr_best);
  myDEBUG && std::cout<<"=======End Best Candidate Selection======"<<endl;
}

//+------------------------------------------------
// write MC hsitogram
//+------------------------------------------------
void DzTohheta::writeMcHist(Particle &Dst, BelleTuple *hist)
{
  myDEBUG && std::cout << "Begin write MC Histogram" << endl;
  Gen_hepevt_Manager &genMgr = Gen_hepevt_Manager::get_manager();

  Particle & D0 = Dst.relation().child(0);
  int chnd0 = dynamic_cast<D0UserInfo&>(D0.userInfo()).m_channel;
  hist->column("mc",GetMcType(Dst));

  //+-----------------
  // slow pion mother(only save 400<|id|<600) info
  // getGenHepInfoF in dmixutil.cc
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  std::vector<moinfo> spMlist;
  getGenHepInfoFeta(Dst.relation().child(1), spMlist);
  int spMsize = spMlist.size();
  myDEBUG && std::cout<<"Save pi_s info.--"<<std::endl;
  hist->column("mcspml1",spMsize>=1 ? spMlist[0].motherLundID : 0);
  hist->column("mcspmi1",spMsize>=1 ? spMlist[0].motherID : 0);
  hist->column("mcspmn1",spMsize>=1 ? spMlist[0].idset.size() : 0);
  hist->column("mcspmn1g",spMsize>=1 ? spMlist[0].idset2.size() : 0);
  hist->column("mcspml2",spMsize>=2 ? spMlist[1].motherLundID : 0);
  hist->column("mcspmi2",spMsize>=2 ? spMlist[1].motherID : 0);
  hist->column("mcspmn2",spMsize>=2 ? spMlist[1].idset.size() : 0);
  hist->column("mcspmn2g",spMsize>=2 ? spMlist[1].idset2.size() : 0);

  myDEBUG && std::cout<<"save chg tracks info."<<std::endl;
  std::vector<moinfo> kMlist, pMlist, p0Mlist;
  // Pi0: getGenHepInfoF in dmixutil.cc
  // Eta: getGenHepInfoFeta in dmixutil.cc
  // 2: with FSR photon
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int photo = -1;
  getGenHepInfoF2eta(D0.relation().child(0), kMlist, &photo);
  getGenHepInfoFeta(D0.relation().child(1), pMlist);
  if(D0.child(2).genHepevt()) getGenHepInfoFeta(D0.relation().child(2), p0Mlist);
  int kMsize = kMlist.size();
  int pMsize = pMlist.size();
  int p0Msize = p0Mlist.size();

  hist->column("photo", photo);
  hist->column("mckml1",kMsize>=1 ? kMlist[0].motherLundID : 0);
  hist->column("mckmi1",kMsize>=1 ? kMlist[0].motherID : 0);
  hist->column("mckmn1",kMsize>=1 ? kMlist[0].idset.size() : 0);
  hist->column("mckmn1g",kMsize>=1 ? kMlist[0].idset2.size() : 0);
  hist->column("mckml2",kMsize>=2 ? kMlist[1].motherLundID : 0);
  hist->column("mckmi2",kMsize>=2 ? kMlist[1].motherID : 0);
  hist->column("mckmn2",kMsize>=2 ? kMlist[1].idset.size() : 0);
  hist->column("mckmn2g",kMsize>=2 ? kMlist[1].idset2.size() : 0);
  hist->column("mckml3",kMsize>=3 ? kMlist[2].motherLundID : 0);
  hist->column("mckmi3",kMsize>=3 ? kMlist[2].motherID : 0);
  hist->column("mckmn3",kMsize>=3 ? kMlist[2].idset.size() : 0);
  hist->column("mckmn3g",kMsize>=3 ? kMlist[2].idset2.size() : 0);
  hist->column("mcpml1",pMsize>=1 ? pMlist[0].motherLundID : 0);
  hist->column("mcpmi1",pMsize>=1 ? pMlist[0].motherID : 0);
  hist->column("mcpmn1",pMsize>=1 ? pMlist[0].idset.size() : 0);
  hist->column("mcpmn1g",pMsize>=1 ? pMlist[0].idset2.size() : 0);
  hist->column("mcpml2",pMsize>=2 ? pMlist[1].motherLundID : 0);
  hist->column("mcpmi2",pMsize>=2 ? pMlist[1].motherID : 0);
  hist->column("mcpmn2",pMsize>=2 ? pMlist[1].idset.size() : 0);
  hist->column("mcpmn2g",pMsize>=2 ? pMlist[1].idset2.size() : 0);
  hist->column("mcpml3",pMsize>=3 ? pMlist[2].motherLundID : 0);
  hist->column("mcpmi3",pMsize>=3 ? pMlist[2].motherID : 0);
  hist->column("mcpmn3",pMsize>=3 ? pMlist[2].idset.size() : 0);
  hist->column("mcpmn3g",pMsize>=3 ? pMlist[2].idset2.size() : 0);
  hist->column("mcp0ml1",p0Msize>=1 ? p0Mlist[0].motherLundID : 0);
  hist->column("mcp0mi1",p0Msize>=1 ? p0Mlist[0].motherID : 0);
  hist->column("mcp0mn1",p0Msize>=1 ? p0Mlist[0].idset.size() : 0);
  hist->column("mcp0mn1g",p0Msize>=1 ? p0Mlist[0].idset2.size() : 0);
  hist->column("mcp0ml2",p0Msize>=2 ? p0Mlist[1].motherLundID : 0);
  hist->column("mcp0mi2",p0Msize>=2 ? p0Mlist[1].motherID : 0);
  hist->column("mcp0mn2",p0Msize>=2 ? p0Mlist[1].idset.size() : 0);
  hist->column("mcp0mn2g",p0Msize>=2 ? p0Mlist[1].idset2.size() : 0);
  hist->column("mcp0ml3",p0Msize>=3 ? p0Mlist[2].motherLundID : 0);
  hist->column("mcp0mi3",p0Msize>=3 ? p0Mlist[2].motherID : 0);
  hist->column("mcp0mn3",p0Msize>=3 ? p0Mlist[2].idset.size() : 0);
  hist->column("mcp0mn3g",p0Msize>=3 ? p0Mlist[2].idset2.size() : 0);

  int mcspi, mcspl, gspi, gspl, gspmi, gspml;
  int mcki, mckl, gki, gkl, gkmi, gkml;
  int mcpi, mcpl, gpi, gpl, gpmi, gpml;
  int mcp0i, mcp0l, gp0i, gp0l, gp0mi, gp0ml;
  HepLorentzVector mcspP, mckP, mcpP, mcp0P;
  HepLorentzVector gspP, gkP, gpP, gp0P;
  HepPoint3D mcspvtx, mckvtx, mcpvtx, mcp0vtx;
  HepPoint3D gspvtx, gkvtx, gpvtx, gp0vtx;
  HepPoint3D mcbvtx, mcdvtx;
  HepLorentzVector mcd0p;
  Hep3Vector mcdP;
  double mcbt = 0, gkp13 = 0.;

  getGenHepInfoS(Dst.child(1),
                 mcspi, mcspl, mcspP, mcspvtx,
                 gspi , gspl, gspP, gspvtx,
                 gspmi, gspml);
  getGenHepInfoS(D0.child(0),
                 mcki, mckl, mckP, mckvtx,
                 gki , gkl, gkP, gkvtx,
                 gkmi, gkml);
  getGenHepInfoS(D0.child(1),
                 mcpi, mcpl, mcpP, mcpvtx,
                 gpi , gpl, gpP, gpvtx,
                 gpmi, gpml);
  getGenHepInfoS(D0.child(2),      
       mcp0i,mcp0l,mcp0P,mcp0vtx,
       gp0i , gp0l, gp0P,gp0vtx,   
       gp0mi,gp0ml); 

  // need to calculate the mass resolution
  if (D0.child(0).genHepevt() && D0.child(1).genHepevt() && D0.child(2).genHepevt()) {
     mcp0P.setX(D0.child(2).genHepevt().PX());
     mcp0P.setY(D0.child(2).genHepevt().PY());
     mcp0P.setZ(D0.child(2).genHepevt().PZ());
     mcp0P.setT(D0.child(2).genHepevt().E());
             
     double mcmkp=(mckP+mcpP).mag();
     double mcmpp0=(mcp0P+mcpP).mag();
     double mcmkp0=(mckP+mcp0P).mag();
             
     hist->column("mcmkp",mcmkp);
     hist->column("mcmkp0",mcmkp0);
     hist->column("mcmpp0",mcmpp0);
  } else {    
     hist->column("mcmkp",0);
     hist->column("mcmkp0",0);
     hist->column("mcmpp0",0);
  }

  hist->column("mcspi", mcspi);
  hist->column("mcspl", mcspl);
  hist->column("gspi",gspi);
  hist->column("gspl",gspl);
  hist->column("gspml", gspml);
  hist->column("mckl", mckl);
  hist->column("gkml", gkml);
  hist->column("mcpl", mcpl);
  hist->column("gpml", gpml);
  hist->column("mcp0l",
               (D0.child(2).genHepevt())?(float) (D0.child(2).genHepevt().idhep()) : 0);
  hist->column("gp0ml", gp0ml ); 

  myDEBUG && std::cout<<"save Decay time "<<std::endl;
  if(D0.genHepevt()) {
    hist->column("mcd01", D0.genHepevt().idhep());
    Gen_hepevt gd0 = gen_level(D0.genHepevt());
    if(gd0) {
      mcd0p = HepLorentzVector(gd0.PX(), gd0.PY(), gd0.PZ(), gd0.E());
      mcdvtx = HepPoint3D(gd0.VX(), gd0.VY(), gd0.VZ());
      mcdP = mcd0p.vect();

      HepPoint3D mcdLen1 = gkvtx*10-mcdvtx; //mm
      HepPoint3D mcdLen2 = gpvtx*10-mcdvtx; //mm

      double mclen1 = mcdLen1*(mcdP.unit()); //mm
      double mclen2 = mcdLen2*(mcdP.unit()); //mm

      hist->column( "mcdp", mcdP.mag() );
      hist->column( "mclen1", mclen1 );
      hist->column( "mclen2", mclen2 );
     }else {
      hist->column( "mcdp", -1. );
      hist->column( "mclen1", 0. );
      hist->column( "mclen2", 0. );
     }
   }
   if(D0.child(2).child(0).genHepevt() && D0.child(2).child(0).genHepevt().mother() ) {
     hist->column("mcg1ml1",
                  (float) (D0.child(2).child(0).genHepevt().mother().idhep()));
     hist->column("mcg1mi1",
                  (int) (D0.child(2).child(0).genHepevt().mother().get_ID()));
  } else {         
     hist->column("mcg1mi1", 0 );
     hist->column("mcglmi1", 0 );
  }
  if(D0.child(2).child(1).genHepevt() && D0.child(2).child(1).genHepevt().mother() ) {
    hist->column("mcg2ml1",
                 (float) (D0.child(2).child(1).genHepevt().mother().idhep()));
    hist->column("mcg2mi1",
                 (int) (D0.child(2).child(1).genHepevt().mother().get_ID()));
  } else {         
    hist->column("mcg2ml1", 0);
    hist->column("mcg2mi1", 0);
  }
   if(D0.child(2).child(0).genHepevt() && D0.child(2).child(0).genHepevt().mother()
      && D0.child(2).child(0).genHepevt().mother().mother()) {
     hist->column("mcg1ml2",
                  (float) (D0.child(2).child(0).genHepevt().mother().mother().idhep()));
     hist->column("mcg1mi2",
                  (int) (D0.child(2).child(0).genHepevt().mother().mother().get_ID()));
  } else {
     hist->column("mcg1mi2", 0 );
     hist->column("mcglmi2", 0 );
  }
  if(D0.child(2).child(1).genHepevt() && D0.child(2).child(1).genHepevt().mother()
     && D0.child(2).child(1).genHepevt().mother().mother()) {
    hist->column("mcg2ml2",
                 (float) (D0.child(2).child(1).genHepevt().mother().mother().idhep()));
    hist->column("mcg2mi2",
                 (int) (D0.child(2).child(1).genHepevt().mother().mother().get_ID()));
  } else {
    hist->column("mcg2ml2", 0);
    hist->column("mcg2mi2", 0);
  }
  myDEBUG && std::cout<<"End writeMcHist--------"<<std::endl;
}

void DzTohheta::writeHist(std::vector<Particle> &plist, BelleTuple *hist)
{

  myDEBUG && std::cout << "Begin write Histogram" << endl;
  if(plist.size() == 0) return;

  int n_ch1(0), n_ch2(0), n_ch3(0), n_ch4(0), n_ch5(0);
  for(unsigned i=0; i<plist.size(); ++i) {
    Particle &Dst = plist[i];
    Particle &D0 = Dst.child(0);
    int chnd0 = dynamic_cast<D0UserInfo&>(D0.userInfo()).m_channel;
    if( abs(chnd0)==1 ) n_ch1++;
    else if( abs(chnd0)==2 ) n_ch2++;
    else if( abs(chnd0)==3 ) n_ch3++;
  }

  for(unsigned i=0; i<plist.size(); ++i) {
    hist->column("ndst_ch1", n_ch1 );
    hist->column("ndst_ch2", n_ch2 );
    hist->column("ndst_ch3", n_ch3 );

    Particle &Dst = plist[i];
    Particle &D0 = Dst.child(0);
    Particle &sp = Dst.child(1);
    Particle &Kaon = D0.child(0);
    Particle &Pion = D0.child(1); 
    Particle &Pion0 = D0.child(2); 

    // event info.
    hist->column("ynsno",  YnsNo  );// YnS #
    hist->column("expno",  m_expNo);// Exp #
    hist->column("runno",  m_runNo);// Run #
    hist->column("evtno",  m_evtNo);// Event #

    int chnd0 = dynamic_cast<D0UserInfo&>(D0.userInfo()).m_channel; 
    hist->column("chnd0", chnd0 );
    hist->column("chndst", dynamic_cast<DstUserInfo&>(Dst.userInfo()).m_channel );
    hist->column("pi0flag1", dynamic_cast<DstUserInfo&>(Dst.userInfo()).m_pi0veto1Flag );
    hist->column("pi0flag2", dynamic_cast<DstUserInfo&>(Dst.userInfo()).m_pi0veto2Flag );
    hist->column("pi0mmin1", dynamic_cast<DstUserInfo&>(Dst.userInfo()).m_pi0veto1Mmin );
    hist->column("pi0mmin2", dynamic_cast<DstUserInfo&>(Dst.userInfo()).m_pi0veto2Mmin );
    hist->column("pi0pmin1", dynamic_cast<DstUserInfo&>(Dst.userInfo()).m_pi0veto1Pmin );
    hist->column("pi0pmin2", dynamic_cast<DstUserInfo&>(Dst.userInfo()).m_pi0veto2Pmin );
    hist->column("pi0asye1", dynamic_cast<DstUserInfo&>(Dst.userInfo()).m_pi0veto1Easy );
    hist->column("pi0asye2", dynamic_cast<DstUserInfo&>(Dst.userInfo()).m_pi0veto2Easy );

    int ikaonf = 3; 
    int ipionf = 2;
    if( abs(chnd0)==2 ) { ikaonf=2; ipionf=2; }
    else if( abs(chnd0)==3 ) { ikaonf=3; ipionf=3; } 
 
    // dmass info.
    hist->column("massd0", dynamic_cast<D0UserInfo&>(D0.userInfo()).mass());
    hist->column("massdif",
                 dynamic_cast<DstUserInfo&>(Dst.userInfo()).massDif() );
    hist->column("md00", (Kaon.p()+Pion.p()+Pion0.p()).mag() );
    hist->column("mdif0", (Kaon.p()+Pion.p()+Pion0.p()+sp.p()).mag()
                 - (Kaon.p()+Pion.p()+Pion0.p()).mag() );
    HepLorentzVector DstP4 = Dst.p(); 
    double dstpsTmp = pStar(DstP4, m_HER, m_LER, 22.0).vect().mag(); 
    hist->column("dstps", dstpsTmp );
    HepLorentzVector DstP42 = Dst.p(); 
    double cosTheta = pStar(DstP42, m_HER, m_LER, 22.0).vect().unit().dot( Hep3Vector(0,0,1) );  
    hist->column("cosdst", cosTheta);

    //for pi0 information
    hist->column("mchip0", dynamic_cast<UserInfo&>(Pion0.userInfo()).mchisq() );
    hist->column("m0p0", dynamic_cast<UserInfo&>(Pion0.userInfo()).mass0() );
    hist->column("pp0", dynamic_cast<UserInfo&>(Pion0.userInfo()).moment() );
    hist->column("ediffp0", dynamic_cast<UserInfo&>(Pion0.userInfo()).ediff() );

    /*=======================Kaon, Pion Info.================================*/
    // SVD hits number: nR and nZ
    // nhits(5): No. of associated hits in CDC+SVD.
    // 0: axial-wire; 1: stereo-wire; 2: cathode; 3: SVD-rphi; 4: SVD-z.
    hist->column("nrk", Kaon.mdstCharged().trk().mhyp(ikaonf).nhits(3));
    hist->column("nzk", Kaon.mdstCharged().trk().mhyp(ikaonf).nhits(4));
    hist->column("nrp", Pion.mdstCharged().trk().mhyp(ipionf).nhits(3));
    hist->column("nzp", Pion.mdstCharged().trk().mhyp(ipionf).nhits(4));

    // Track Fitting Info.
    hist->column("drk", Kaon.mdstCharged().trk().mhyp(ikaonf).helix(0));
    hist->column("dzk", Kaon.mdstCharged().trk().mhyp(ikaonf).helix(3));
    hist->column("drp", Pion.mdstCharged().trk().mhyp(ipionf).helix(0));
    hist->column("dzp", Pion.mdstCharged().trk().mhyp(ipionf).helix(3));

    // K/pi ID
    hist->column("kaonid", kaonId2(Kaon));
    hist->column("pionid", kaonId2(Pion));
    // Proton/K
    hist->column("kpnid", protonId(Kaon));
    // moun Id
    hist->column("kmuid", muonId(Kaon));
    hist->column("pmuid", muonId(Pion));
    // Eid
    hist->column("keid", eId(Kaon));
    hist->column("peid", eId(Pion));

    hist->column("kp", Kaon.ptot());
    hist->column("pp", Pion.ptot());
    hist->column("p0p", Pion0.ptot());
    hist->column("kpxy", Kaon.p().vect().perp());
    hist->column("ppxy", Pion.p().vect().perp());
    hist->column("p0pxy", Pion0.p().vect().perp());
    double cosk_tmp  = Kaon.p().vect().unit().dot( Hep3Vector(0,0,1) );
    double cosp_tmp  = Pion.p().vect().unit().dot( Hep3Vector(0,0,1) );
    double cosp0_tmp = Pion0.p().vect().unit().dot( Hep3Vector(0,0,1) );
    hist->column("cosk", cosk_tmp );
    hist->column("cosp", cosp_tmp );
    hist->column("cosp0", cosp0_tmp );

    HepLorentzVector tmpK_P4 = Kaon.p(); 
    HepLorentzVector tmpPi_P4 = Pion.p(); 
    HepLorentzVector tmpPi0_P4 = Pion0.p(); 
    hist->column("kps", pStar(tmpK_P4, m_HER, m_LER, 22.0).vect().mag());
    hist->column("pps", pStar(tmpPi_P4, m_HER, m_LER, 22.0).vect().mag());
    hist->column("p0ps", pStar(tmpPi0_P4, m_HER, m_LER, 22.0).vect().mag());

    // PI0 and Gamma info
    Mdst_gamma gm0 = Pion0.child(0).relation().mdstGamma() ;
    Mdst_gamma gm1 = Pion0.child(1).relation().mdstGamma() ;
    HepLorentzVector gm0p4 = Pion0.child(0).p();
    HepLorentzVector gm1p4 = Pion0.child(1).p();
    Hep3Vector gm0p3(gm0p4.px(), gm0p4.py(), gm0p4.pz() );
    Hep3Vector gm1p3(gm1p4.px(), gm1p4.py(), gm1p4.pz() );
    hist->column("cosgmb1", gm0p4.vect().unit().dot( Hep3Vector(0,0,1)));
    hist->column("cosgmb2", gm1p4.vect().unit().dot( Hep3Vector(0,0,1)));
    hist->column("mpi0", Pion0.p().mag() );
    hist->column("egm1", gm0p3.mag() );
    hist->column("egm2", gm1p3.mag() );
    Mdst_ecl_aux_Manager& eclaux_mgr = Mdst_ecl_aux_Manager::get_manager();
    hist->column("e9251", eclaux_mgr( Panther_ID(gm0.ecl().get_ID()) ).e9oe25() );
    hist->column("e9252", eclaux_mgr( Panther_ID(gm1.ecl().get_ID()) ).e9oe25() );
    gm0p4.boost(-Pion0.p().boostVector());
    hist->column("cosgm", gm0p4.vect().unit().dot(Pion0.p().vect().unit()));


    /*======================= D0 info ============================*/
    HepPoint3D dvtx = D0.momentum().decayVertex();
    HepPoint3D bvtx = D0.momentum().vertex();
    HepLorentzVector d0p = D0.momentum().p();

    hist->column("d0px", d0p.x());
    hist->column("d0py", d0p.y());
    hist->column("d0pz", d0p.z());
    hist->column("d0p", D0.ptot());
    hist->column("d0pxy", d0p.vect().perp());
    HepLorentzVector D0P4 = D0.momentum().p(); 
    double d0psTmp = pStar(D0P4, m_HER, m_LER, 22.0).vect().mag(); 
    hist->column("d0ps", d0psTmp);
    HepLorentzVector D0P42 = D0.momentum().p(); 
    double cosThetaD0 = pStar(D0P42, m_HER, m_LER, 22.0).vect().unit().dot( Hep3Vector(0,0,1) );  
    hist->column("cosd0", cosThetaD0 );

    double m12(-1), m13(-1), m23(-1);
    Dalitz( D0, m12, m13, m23 );
    hist->column("m13m", m13 );
    hist->column("m23m", m23 );
    hist->column("m12m", m12 );

    Particle &pKp = D0.child(0);
    Particle &pKm = D0.child(1);
    hist->column("m13", (pKp.p()+Pion0.p()).mag() );
    hist->column("m23", (pKm.p()+Pion0.p()).mag() );
    hist->column("m12", (pKp.p()+pKm.p()).mag());
    hist->column("m123", (pKp.p()+pKm.p()+Pion0.p()).mag());

    HepLorentzVector KaonP = D0.child(0).p();
    HepLorentzVector PionP = D0.child(1).p();
    HepLorentzVector Pion0P = D0.child(2).p();;
    KaonP.boost(-D0.p().boostVector());
    PionP.boost(-D0.p().boostVector());
    Pion0P.boost(-D0.p().boostVector());
    // decay angle
    hist->column("cos1", KaonP.vect().unit().dot(D0.p().vect().unit()));
    hist->column("cos2", PionP.vect().unit().dot(D0.p().vect().unit()));
    hist->column("cos3", Pion0P.vect().unit().dot(D0.p().vect().unit()));


    HepLorentzVector KaonP2 = D0.child(0).p();
    HepLorentzVector PionP2 = D0.child(1).p(); 
    HepLorentzVector Pion0P2 = D0.child(2).p(); 

    HepLorentzVector genKpPi0P = KaonP2 + Pion0P2 ;
    HepLorentzVector PionP3 = D0.child(1).p();
    HepLorentzVector genKmPi0P = PionP2 + Pion0P2 ;
    HepLorentzVector KaonP3 = D0.child(0).p();
    HepLorentzVector genKpKmP = KaonP2 + PionP2 ;
    HepLorentzVector Pion0P3 = D0.child(2).p();
    HepLorentzVector KaonP4 = KaonP2 ;
    // helicity angle
    KaonP2.boost( -genKpPi0P.boostVector() );
    //add two lines below on Sep 12 (after published...)
    PionP3.boost( -D0.p().boostVector() );
    genKpPi0P.boost( -D0.p().boostVector() ); 
    PionP3.boost( -genKpPi0P.boostVector() );
    hist->column("cosh13", KaonP2.vect().unit().dot(PionP3.vect().unit()));

    PionP2.boost( -genKmPi0P.boostVector() );
    //add two lines below on Sep 12 (after published...)
    KaonP3.boost( -D0.p().boostVector() );
    genKmPi0P.boost( -D0.p().boostVector() );
    KaonP3.boost( -genKmPi0P.boostVector() );
    hist->column("cosh23", PionP2.vect().unit().dot(KaonP3.vect().unit()));

    KaonP4.boost( -genKpKmP.boostVector() );
    //add two lines below on Sep 12 (after published...)
    Pion0P3.boost( -D0.p().boostVector() );
    genKpKmP.boost( -D0.p().boostVector() );
    Pion0P3.boost( -genKpKmP.boostVector() );
    hist->column("cosh12", KaonP4.vect().unit().dot(Pion0P3.vect().unit()));

    // kinematic fit quality
    //hist->column("clv", dynamic_cast<D0UserInfo&>(D0.userInfo()).clV());
    hist->column("chisqv", dynamic_cast<D0UserInfo&>(D0.userInfo()).chisqV());
    //hist->column("ndfv", dynamic_cast<D0UserInfo&>(D0.userInfo()).ndfV());
    myDEBUG && cout<<"ndfv: "<< dynamic_cast<D0UserInfo&>(D0.userInfo()).ndfV() <<endl; 
    //hist->column("clb", dynamic_cast<D0UserInfo&>(D0.userInfo()).clB());
    hist->column("chisqb", dynamic_cast<D0UserInfo&>(D0.userInfo()).chisqB() );
    //hist->column("ndfb", dynamic_cast<D0UserInfo&>(D0.userInfo()).ndfB());
    myDEBUG && cout<<"ndfb: "<< dynamic_cast<D0UserInfo&>(D0.userInfo()).ndfB() <<endl; 
    //hist->column("cls", dynamic_cast<DstUserInfo&>(Dst.userInfo()).clS());
    hist->column("chisqs", dynamic_cast<DstUserInfo&>(Dst.userInfo()).chisqS());
    //hist->column("ndfs", dynamic_cast<DstUserInfo&>(Dst.userInfo()).ndfS());
    myDEBUG && cout<<"ndfs: "<< dynamic_cast<DstUserInfo&>(Dst.userInfo()).ndfS() <<endl; 
    double sumchisq = dynamic_cast<D0UserInfo&>(D0.userInfo()).chisqV() + dynamic_cast<D0UserInfo&>(D0.userInfo()).chisqB() + dynamic_cast<DstUserInfo&>(Dst.userInfo()).chisqS() ; 
    hist->column("sumchisq", sumchisq );

    // Slow Pion Info.
    HepLorentzVector spP = sp.p();
    hist->column("spp", spP.rho());
    hist->column("sppxy", spP.vect().perp());
    double sppsTmp = pStar(spP,m_HER,m_LER,22.0).vect().mag(); 
    hist->column("spps", sppsTmp );
    HepLorentzVector spP2 = sp.p();
    double cosspTmp = spP2.vect().unit().dot( Hep3Vector(0,0,1) );
    hist->column("cossp", cosspTmp );

    hist->column("nrsp", sp.mdstCharged().trk().mhyp(IPION).nhits(3));
    hist->column("nzsp", sp.mdstCharged().trk().mhyp(IPION).nhits(4));
    //hist->column("chisqsp", sp.mdstCharged().trk().mhyp(IPION).chisq());
    //hist->column("ndfsp", sp.mdstCharged().trk().mhyp(IPION).ndf());
    hist->column("drsp", sp.mdstCharged().trk().mhyp(IPION).helix(0));
    hist->column("dzsp", sp.mdstCharged().trk().mhyp(IPION).helix(3));
    hist->column("spkid", kaonId2(sp));
    hist->column("speid", eId(sp));
    hist->column("spmuid", muonId(sp));

    HepLorentzVector spP3 = sp.p();
    hist->column("cosspd", spP3.vect().unit().dot(D0.p().vect().unit()));
    HepLorentzVector spP4 = sp.p();
    HepLorentzVector DstP43 = Dst.p(); 
    spP4.boost(-DstP43.boostVector());
    hist->column("cosdedst", spP4.vect().unit().dot(DstP43.vect().unit()));   

    //D0 lifetime and decay length(3D)
    HepPoint3D dLen = dvtx - bvtx; //cm
    Hep3Vector dP = d0p.vect();
    double leng = dLen*(dP.unit()); //cm
    hist->column("leng",leng*10000.);//um

    // D0 lifetime and decay length(2D)
    dLen.setZ(0.); dP.setZ(0.);
    double lenxy  = dLen*(dP.unit());//cm
    hist->column("lxy", lenxy*10000.);//um

    // D0 lifetime error (3D and 2D)
    double errT, errTxy, errTx, errTy, errTz;
    // D0 decay length error(3D and 2D)
    double errL, errLxy, errLx, errLy, errLz;
    errDecayLength(errL,errLxy,errLx,errLy,errLz,
                   dvtx,bvtx,d0p,
                   D0.momentum().dVertex(),
                   D0.momentum().dpx(),
                   dynamic_cast<D0UserInfo&>(D0.userInfo()).errBP(),1);
    if(errL>=0.)  hist->column("el", sqrt(errL)*10000.);//um, with correlation
    else hist->column("el", -sqrt(-errL)*10000.);//um
    if(errLxy>=0.)  hist->column("elxy", sqrt(errLxy)*10000.);//um
    else hist->column("elxy", -sqrt(-errLxy)*10000.);//um

    if(mc2exp){
      // MC Data
      if( D0.genHepevt() && !(D0.genHepevt()) )  std::cout << "Error: MC matching(1)" << std::endl;
      if( !(D0.genHepevt()) && D0.genHepevt() )  std::cout << "Error: MC matching(2)" << std::endl;
      writeMcHist(Dst, hist);
    } else {
      // Real Data
      hist->column("mc", -1);
    }
    hist->dumpData();
  }
}

//+-------------------------------------------------------
// erase particle list
//+-------------------------------------------------------
void DzTohheta::endEvent(void)
{
  eraseVector(m_kaonP);
  eraseVector(m_kaonM);
  eraseVector(m_pionP);
  eraseVector(m_pionM);
  eraseVector(m_pionP2);
  eraseVector(m_pionM2);
  eraseVector(m_gamma);
  eraseVector(m_eta);
  eraseVector(m_pi0);
}

//+--------------------------------------------------
// D0 decay length error matrix 
//+-------------------------------------------------
void DzTohheta::errDecayLength(double &errL, double &errLxy,
                    double &errLx, double &errLy, double &errLz,
                    const HepPoint3D &dec, const HepPoint3D &ip, const HepLorentzVector &p,
                    const HepSymMatrix &errIp, const HepSymMatrix &errP, const HepMatrix &covIPP,
                    int flag)
{
  double f   = 1.0/p.vect().mag();
  double fxy = 1.0/p.vect().perp();
  double fx  = 1.0/p.x();
  double fy  = 1.0/p.y();
  double fz  = 1.0/p.z();
  //
  //     |  Dec  : D-IP :  D-M  |
  //     |                      |
  // V = | IP-D  :  IP  : IP-M  |
  //     |                      |
  //     |  M-D  : M-IP :  Mom  |
  //
  HepSymMatrix V(9,0);
  V.sub(1, errP.sub(5,7));
  V.sub(4, errIp);
  V.sub(7, errP.sub(1,3));
  if(flag){
    V[0][6] = errP[4][0];//x-px
    V[0][7] = errP[4][1];//x-py
    V[0][8] = errP[4][2];//x-pz
    V[1][6] = errP[5][0];//y-px
    V[1][7] = errP[5][1];//y-py
    V[1][8] = errP[5][2];//y-pz
    V[2][6] = errP[6][0];//z-px
    V[2][7] = errP[6][1];//z-py
    V[2][8] = errP[6][2];//z-pz

    V[3][6] = covIPP[0][0];//ipx-px
    V[3][7] = covIPP[0][1];//ipx-py
    V[3][8] = covIPP[0][2];//ipx-pz
    V[4][6] = covIPP[1][0];//ipy-px
    V[4][7] = covIPP[1][1];//ipy-py
    V[4][8] = covIPP[1][2];//ipy-pz
    V[5][6] = covIPP[2][0];//ipz-px
    V[5][7] = covIPP[2][1];//ipz-py
    V[5][8] = covIPP[2][2];//ipz-pz

    V[3][0] = covIPP[0][4];//ipx-x
    V[3][1] = covIPP[0][5];//ipx-y
    V[3][2] = covIPP[0][6];//ipx-z
    V[4][0] = covIPP[1][4];//ipy-x
    V[4][1] = covIPP[1][5];//ipy-y
    V[4][2] = covIPP[1][6];//ipy-z
    V[5][0] = covIPP[2][4];//ipz-x
    V[5][1] = covIPP[2][5];//ipz-y
    V[5][2] = covIPP[2][6];//ipz-z
  }

  HepMatrix dtdV(9,1,0);
  dtdV[0][0] =  f*p.x();
  dtdV[1][0] =  f*p.y();
  dtdV[2][0] =  f*p.z();
  dtdV[3][0] = -f*p.x();
  dtdV[4][0] = -f*p.y();
  dtdV[5][0] = -f*p.z();
  dtdV[6][0] = f*(dec-ip).x()-f*((dec-ip)*(p.vect()))*p.x()/p.vect().mag2();
  dtdV[7][0] = f*(dec-ip).y()-f*((dec-ip)*(p.vect()))*p.y()/p.vect().mag2();
  dtdV[8][0] = f*(dec-ip).z()-f*((dec-ip)*(p.vect()))*p.z()/p.vect().mag2();
  errL = ((dtdV.T())*V*dtdV)[0][0];

  Hep3Vector momP = p.vect();
  momP.setZ(0.);
  HepMatrix dtxydV(9,1,0);
  dtxydV[0][0] =  fxy*p.x();
  dtxydV[1][0] =  fxy*p.y();
  dtxydV[2][0] =  0.; // fxy*p.z();
  dtxydV[3][0] = -fxy*p.x();
  dtxydV[4][0] = -fxy*p.y();
  dtxydV[5][0] =  0.; // -fxy*p.z();
  dtxydV[6][0] =  fxy*(dec-ip).x()-fxy*((dec-ip)*momP)*p.x()/p.vect().perp2();
  dtxydV[7][0] =  fxy*(dec-ip).y()-fxy*((dec-ip)*momP)*p.y()/p.vect().perp2();
  dtxydV[8][0] =  0.; // fxy*(dec-ip).z()-fxy*((dec-ip)*momP)*p.z()/p.vect().perp2();
  errLxy = ((dtxydV.T())*V*dtxydV)[0][0];

  // x
  HepMatrix dtxdV(9,1,0);
  dtxdV[0][0] =  1.;
  dtxdV[1][0] =  0.; // fx*p.y();
  dtxdV[2][0] =  0.; // fx*p.z();
  dtxdV[3][0] = -1.;
  dtxdV[4][0] =  0.; // -fx*p.y();
  dtxdV[5][0] =  0.; // -fx*p.z();
  dtxdV[6][0] =  0.;
  dtxdV[7][0] =  0.; //fx*(dec-ip).y()-fx*((dec-ip)*(p.vect()))*p.y()/(p.x()*p.x());
  dtxdV[8][0] =  0.; //fx*(dec-ip).z()-fx*((dec-ip)*(p.vect()))*p.z()/(p.x()*p.x());
  errLx = ((dtxdV.T())*V*dtxdV)[0][0];

  // y
  HepMatrix dtydV(9,1,0);
  dtydV[0][0] =  0.; // fy*p.x();
  dtydV[1][0] =  1.;
  dtydV[2][0] =  0.; // fy*p.z();
  dtydV[3][0] =  0.; // -fy*p.x();
  dtydV[4][0] =  -1.;
  dtydV[5][0] =  0.; // -fy*p.z();
  dtydV[6][0] =  0.; //fy*(dec-ip).x()-fy*((dec-ip)*(p.vect()))*p.x()/(p.y()*p.y());
  dtydV[7][0] =  0.;
  dtydV[8][0] =  0.; //fy*(dec-ip).z()-fy*((dec-ip)*(p.vect()))*p.z()/(p.y()*p.y());
  errLy = ((dtydV.T())*V*dtydV)[0][0];

  // z
  HepMatrix dtzdV(9,1,0);
  dtzdV[0][0] =  0.; // fz*p.x();
  dtzdV[1][0] =  0.; // fz*p.y();
  dtzdV[2][0] =  1.;
  dtzdV[3][0] =  0.; // -fz*p.x();
  dtzdV[4][0] =  0.; // -fz*p.y();
  dtzdV[5][0] = -1.;
  dtzdV[6][0] =  0.; //fz*(dec-ip).x()-fz*((dec-ip)*(p.vect()))*p.x()/(p.z()*p.z());
  dtzdV[7][0] =  0.; //fz*(dec-ip).y()-fz*((dec-ip)*(p.vect()))*p.y()/(p.z()*p.z());
  dtzdV[8][0] =  0.;
  errLz = ((dtzdV.T())*V*dtzdV)[0][0];
  myDEBUG && std::cout<<"errz? 2_5 = " <<V[2][5] << "  "<<V[5][2]<<std::endl;

  if(fabs(errL)   > 1.e+30)errL   = 1.e+30;
  if(fabs(errLxy) > 1.e+30)errLxy = 1.e+30;
  if(fabs(errLx)  > 1.e+30)errLx  = 1.e+30;
  if(fabs(errLy)  > 1.e+30)errLy  = 1.e+30;
  if(fabs(errLz)  > 1.e+30)errLz  = 1.e+30;
}

//+--------------------------------------------------
// D0 mass constraint fit for Dalitz plot variables
// to improve the mass resolution of DP 
//+-------------------------------------------------
int DzTohheta::Dalitz(Particle &p, double &m12, double &m13,
                 double &m23)
{
  kmassfitter km;
  km.invariantMass(p.pType().mass());
  for(unsigned i=0; i<p.nChildren(); ++i)
    addTrack2fit(km,p.child(i));
  if(!km.fit()) {
      m13 = (km.momentum(0)+km.momentum(2)).mag();
      m23 = (km.momentum(1)+km.momentum(2)).mag();
      m12 = (km.momentum(0)+km.momentum(1)).mag();
    return 1;
  }
  return 0;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
