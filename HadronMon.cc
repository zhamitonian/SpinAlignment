// -*- C++ -*-
//
// Package:     <package>
// Module:      HadronMon.cc
// 
// Description: <one line class summary>
//
// Usage:
//    <usage>
//
// Author:      Adachi Ichiro
// Created:     Mon Dec 06 11:31:44 JST 1999
// $Id: HadronMon.cc 10375 2008-01-22 23:54:23Z katayama $
//
// Revision history
//
//

#include "belle.h"
#include <fstream>
#include <list>

#include "belleCLHEP/Alist/AList.h"
#include "belleCLHEP/Vector/LorentzVector.h"

#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"

#include "basf/module.h"
#include "basf/module_descr.h"

#include "panther/panther.h"

#include TRG_H
#include MDST_H
#include EVTCLS_H
#include EVTVTX_H
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


//*---
// static constants.
//*---
static const float ELER(3.5);
static const float EHER(7.9965); // didn't count different C.M. energy , not resonable ?
static const float theta(0.022);
static const float mass_pi(0.13956995);
static const float pi(3.141592653589);    
static const int   hbase(16300);

// Class definition
class HadronMon : public Module {
 public:
  HadronMon ( void );
  ~HadronMon ( void ){};
  void init ( int *status );
  void term ( void );
  void hist_def ( void );
  void disp_stat ( const char* ) {} ;
  void event ( BelleEvent*, int* );
  void begin_run ( BelleEvent*, int *status ){ *status = 0; };
  void end_run ( BelleEvent*, int *status ){ *status = 0; };
  void other ( int*, BelleEvent*, int* ){};

 public:
  int MonitorLevel;

 private:
  HepLorentzVector CMVector;
  Hep3Vector CMboost;

 private:
  BelleHistogram* H[100];

};


// BASF Interface Function 
extern "C" Module_descr *mdcl_HadronMon ()
{
/* main */

  HadronMon *module = new HadronMon;
  Module_descr *dscr = new Module_descr ( "HadronMon", module );

  //*---
  // monitoring level.
  //*---
  dscr->define_param ( "MonitorLevel",
                       "monitoring level(=3:both/1:HadronC/2:HadronA)",
                       &module->MonitorLevel );

  return dscr;

}

// Member Functions

// Constructor
HadronMon::HadronMon ( void ) 
{

  //*---
  // monitor level.
  //*---
  MonitorLevel = 3;

  //*---
  // Boost back vector.
  //*---
  CMVector=HepLorentzVector(-EHER*sin(theta), 0., -EHER*cos(theta)+ELER,
			    EHER+ELER );
  CMboost=CMVector.boostVector();


}

// init function
void HadronMon::init ( int *status ) 
{
  *status = 0; 
  dout(Debugout::INFO,"HadronMon") << "------------------------------" << std::endl;
  dout(Debugout::INFO,"HadronMon") << "%-HadronMon initialization.      "<<std::endl;
  dout(Debugout::INFO,"HadronMon") << "------------------------------" << std::endl;
  dout(Debugout::INFO,"HadronMon") << "    MonitorLevel : "<<MonitorLevel<<std::endl;
  dout(Debugout::INFO,"HadronMon") << "------------------------------" << std::endl;

}

void HadronMon::term (void)
{

}

// hist_def function
void HadronMon::hist_def ( void )
{
  extern BelleTupleManager* BASF_Histogram;
  BelleTupleManager* tm = BASF_Histogram;

  //*---
  // HadronC events.
  //*---
  H[10] = tm->histogram ( "HadronC:good track:dr",   100, -2.0,  2.0, hbase+10 );
  H[11] = tm->histogram ( "HadronC:good track:dz",   160, -4.0,  4.0 );
  H[12] = tm->histogram ( "HadronC:P",               100, 0.0, 10.0 );
  H[13] = tm->histogram ( "HadronC:PT",              100, 0.0, 10.0 );
  H[14] = tm->histogram ( "HadronC:cos(track)",      100,-1.0, 1.0 );
  H[15] = tm->histogram ( "HadronC:phi(track)",      180, 0.0,360.0 );
  H[16] = tm->histogram ( "HadronC:E",               100, 0.0, 10.0 );

  H[20] = tm->histogram ( "HadronC:trigger bits",     51,-1.0, 50.0, hbase+20 );
  H[21] = tm->histogram ( "HadronC:good tracks",      20, 0.0, 20.0 );
  H[22] = tm->histogram ( "HadronC:good clusters",    20, 0.0, 20.0 );
  H[23] = tm->histogram ( "HadronC:Momentum Sum(CM)",120, 0.0, 12.0 );
  H[24] = tm->histogram ( "HadronC:Energy Sum(CM)",  120, 0.0, 12.0 );
  H[25] = tm->histogram ( "HadronC:Gamma Sum(CM)",   120, 0.0, 12.0 );
  H[26] = tm->histogram ( "HadronC:Visible E(CM)",   150, 0.0, 15.0 );
  H[27] = tm->histogram ( "HadronC:Pz(CM)",          120,-3.0,  3.0 );
  H[28] = tm->histogram ( "HadronC:R2",               50, 0.0,  1.0 );
  H[29] = tm->histogram ( "HadronC:R2_cut",           50, 0.0,  1.0 );
  H[30] = tm->histogram ( "HadronC:Primary Vertex:r",100, 0.0,  2.0 );
  H[31] = tm->histogram ( "HadronC:Primary Vertex:z",100,-4.0,  4.0 );

  //*---
  // HadronA events.
  //*---
  H[40] = tm->histogram ( "HadronA:good track:dr",   100, -2.0,  2.0, hbase+40 );
  H[41] = tm->histogram ( "HadronA:good track:dz",   160, -4.0,  4.0 );
  H[42] = tm->histogram ( "HadronA:P",               100, 0.0, 10.0 );
  H[43] = tm->histogram ( "HadronA:PT",              100, 0.0, 10.0 );
  H[44] = tm->histogram ( "HadronA:cos(track)",      100,-1.0, 1.0 );
  H[45] = tm->histogram ( "HadronA:phi(track)",      180, 0.0,360.0 );
  H[46] = tm->histogram ( "HadronA:E",               100, 0.0, 10.0 );

  H[50] = tm->histogram ( "HadronA:trigger bits",     51,-1.0, 50.0, hbase+50 );
  H[51] = tm->histogram ( "HadronA:good tracks",      20, 0.0, 20.0 );
  H[52] = tm->histogram ( "HadronA:good clusters",    20, 0.0, 20.0 );
  H[53] = tm->histogram ( "HadronA:Momentum Sum(CM)",120, 0.0, 12.0 );
  H[54] = tm->histogram ( "HadronA:Energy Sum(CM)",  120, 0.0, 12.0 );
  H[55] = tm->histogram ( "HadronA:Gamma Sum(CM)",   120, 0.0, 12.0 );
  H[56] = tm->histogram ( "HadronA:visible E(CM)",   150, 0.0, 15.0 );
  H[57] = tm->histogram ( "HadronA:Pz(CM)",          120,-3.0,  3.0 );
  H[58] = tm->histogram ( "HadronA:R2",               50, 0.0,  1.0 );
  H[59] = tm->histogram ( "HadronA:R2_cut",           50, 0.0,  1.0 );
  H[60] = tm->histogram ( "HadronA:Primary Vertex:r",100, 0.0,  2.0 );
  H[61] = tm->histogram ( "HadronA:Primary Vertex:z",100,-4.0,  4.0 );

  //*---
  // HadronB events.
  //*---
  H[70] = tm->histogram ( "HadronB:good track:dr",   100, -2.0,  2.0, hbase+70 );
  H[71] = tm->histogram ( "HadronB:good track:dz",   160, -4.0,  4.0 );
  H[72] = tm->histogram ( "HadronB:P",               100, 0.0, 10.0 );
  H[73] = tm->histogram ( "HadronB:PT",              100, 0.0, 10.0 );
  H[74] = tm->histogram ( "HadronB:cos(track)",      100,-1.0, 1.0 );
  H[75] = tm->histogram ( "HadronB:phi(track)",      180, 0.0,360.0 );
  H[76] = tm->histogram ( "HadronB:E",               100, 0.0, 10.0 );

  H[80] = tm->histogram ( "HadronB:trigger bits",     51,-1.0, 50.0, hbase+80 );
  H[81] = tm->histogram ( "HadronB:good tracks",      20, 0.0, 20.0 );
  H[82] = tm->histogram ( "HadronB:good clusters",    20, 0.0, 20.0 );
  H[83] = tm->histogram ( "HadronB:Momentum Sum(CM)",120, 0.0, 12.0 );
  H[84] = tm->histogram ( "HadronB:Energy Sum(CM)",  120, 0.0, 12.0 );
  H[85] = tm->histogram ( "HadronB:Gamma Sum(CM)",   120, 0.0, 12.0 );
  H[86] = tm->histogram ( "HadronB:visible E(CM)",   150, 0.0, 15.0 );
  H[87] = tm->histogram ( "HadronB:Pz(CM)",          120,-3.0,  3.0 );
  H[88] = tm->histogram ( "HadronB:R2",               50, 0.0,  1.0 );
  H[89] = tm->histogram ( "HadronB:R2_cut",           50, 0.0,  1.0 );
  H[90] = tm->histogram ( "HadronB:Primary Vertex:r",100, 0.0,  2.0 );
  H[91] = tm->histogram ( "HadronB:Primary Vertex:z",100,-4.0,  4.0 );

}

// event function
void HadronMon::event ( BelleEvent* evptr, int* status )
{

  //*---
  // return if not evtcls tables.
  //*---
  if( ! BsCouTab(EVTCLS_FLAG) ) return;

  //*---
  // Check trigger table for MC events.
  //*---
  int Ntrig  = BsCouTab(RECTRG_SUMMARY);
  int Ntrig3 = BsCouTab(RECTRG_SUMMARY3);

  //*---
  // Get managers.
  //*---
  Evtcls_flag_Manager& EvtflagMgr = Evtcls_flag_Manager::get_manager();
  Evtcls_flag_Manager::iterator it_evtcls_flag = EvtflagMgr.begin();

  Evtcls_flag2_Manager& Evtflag2Mgr = Evtcls_flag2_Manager::get_manager();
  Evtcls_flag2_Manager::iterator it_evtcls_flag2 = Evtflag2Mgr.begin();

  Evtcls_hadronic_flag_Manager& HadTagMgr 
    = Evtcls_hadronic_flag_Manager::get_manager();
  Evtcls_hadronic_flag_Manager::iterator it_hadronic_flag = HadTagMgr.begin();

  //hadron
  Evtcls_hadron_info_Manager& EvtHadInfoMgr
    = Evtcls_hadron_info_Manager::get_manager();
  Evtcls_hadron_info_Manager::iterator it_hadron_info = EvtHadInfoMgr.begin();
  Evtcls_hadron_charged_Manager& EvtHadChgMgr
    = Evtcls_hadron_charged_Manager::get_manager();
  Evtcls_hadron_neutral_Manager& EvtHadClsMgr
    = Evtcls_hadron_neutral_Manager::get_manager();

  //mdst_ecl
  Mdst_ecl_Manager& MdstEclMgr = Mdst_ecl_Manager::get_manager();

  //triggger
  Rectrg_summary_Manager& RecTrgMgr = Rectrg_summary_Manager::get_manager();
  Rectrg_summary_Manager::iterator it_rectrg = RecTrgMgr.begin();
  Rectrg_summary3_Manager& RecTrg3Mgr = Rectrg_summary3_Manager::get_manager();
  Rectrg_summary3_Manager::iterator it_rectrg3 = RecTrg3Mgr.begin();

  //primary vertex
  Evtvtx_primary_vertex_Manager& EvtVtxMgr 
    = Evtvtx_primary_vertex_Manager::get_manager();
  Evtvtx_primary_vertex_Manager::iterator it_evtvtx = EvtVtxMgr.begin();

  //*---
  // protection.
  //*---
  if( it_evtcls_flag == EvtflagMgr.end() || !(*it_evtcls_flag) ) return;
  if( (*it_evtcls_flag).flag(0)<10 )   return;
  int Nvtx =  BsCouTab(EVTVTX_PRIMARY_VERTEX);
  if( Nvtx ){
    if( it_evtvtx == EvtVtxMgr.end() || !(*it_evtvtx) ) return;
  }

  //*---
  // check trigger bits.
  //*---
  int bit[48];
  if(  Ntrig && it_rectrg != RecTrgMgr.end() && (*it_rectrg)){
    for( int i=1; i<9; i++ ){
      bit[4*(i-1)  ] = ( ( (*it_rectrg).final(0)>> 4*(i-1) ) & 0x01 );
      bit[4*(i-1)+1] = ( ( (*it_rectrg).final(0)>> 4*(i-1) ) & 0x02 );
      bit[4*(i-1)+2] = ( ( (*it_rectrg).final(0)>> 4*(i-1) ) & 0x04 );
      bit[4*(i-1)+3] = ( ( (*it_rectrg).final(0)>> 4*(i-1) ) & 0x08 );
    }
    for( int i=1; i<5; i++ ){
      bit[4*(i-1)+32] = ( ( (*it_rectrg).final(1)>> 4*(i-1) ) & 0x01 );
      bit[4*(i-1)+33] = ( ( (*it_rectrg).final(1)>> 4*(i-1) ) & 0x02 );
      bit[4*(i-1)+34] = ( ( (*it_rectrg).final(1)>> 4*(i-1) ) & 0x04 );
      bit[4*(i-1)+35] = ( ( (*it_rectrg).final(1)>> 4*(i-1) ) & 0x08 );
    }
  }else if(  Ntrig3 && it_rectrg3 != RecTrg3Mgr.end() && (*it_rectrg3)){
    for( int i=1; i<9; i++ ){
      bit[4*(i-1)  ] = ( ( (*it_rectrg3).final(0)>> 4*(i-1) ) & 0x01 );
      bit[4*(i-1)+1] = ( ( (*it_rectrg3).final(0)>> 4*(i-1) ) & 0x02 );
      bit[4*(i-1)+2] = ( ( (*it_rectrg3).final(0)>> 4*(i-1) ) & 0x04 );
      bit[4*(i-1)+3] = ( ( (*it_rectrg3).final(0)>> 4*(i-1) ) & 0x08 );
    }
    for( int i=1; i<5; i++ ){
      bit[4*(i-1)+32] = ( ( (*it_rectrg3).final(1)>> 4*(i-1) ) & 0x01 );
      bit[4*(i-1)+33] = ( ( (*it_rectrg3).final(1)>> 4*(i-1) ) & 0x02 );
      bit[4*(i-1)+34] = ( ( (*it_rectrg3).final(1)>> 4*(i-1) ) & 0x04 );
      bit[4*(i-1)+35] = ( ( (*it_rectrg3).final(1)>> 4*(i-1) ) & 0x08 );
    }
  }else{
    for( int i=0;i<48;i++ ) bit[i] = 0;
  }

  //*---
  // HadronA or C case.
  //*---
  if( (*it_evtcls_flag).flag(0) >= 10 ){

    //*---
    // charged track.
    //*---
    float Psum  = 0.0;
    float Pzsum = 0.0;
    float Qsum  = 0.0;
    for(Evtcls_hadron_charged_Manager::iterator it = EvtHadChgMgr.begin();
	it != EvtHadChgMgr.end(); ++it ){
      
      Mdst_charged& Charged = (*it).charged();
      Mdst_trk& Trk = Charged.trk();
      Mdst_trk_fit& Trkfit  = Trk.mhyp(2);
      Hep3Vector P3( Charged.px(), Charged.py(), Charged.pz() );
      HepLorentzVector P4( P3, sqrt(P3.mag2()+mass_pi*mass_pi ) );
      HepLorentzVector P4Rest = P4;
      P4Rest.boost(CMboost);
      Hep3Vector P3Rest = P4Rest.vect();
      Psum  += P3Rest.mag();
      Qsum  += Charged.charge();
      Pzsum += P3Rest.z();
	  
      if( MonitorLevel == 1 || MonitorLevel == 3 ){
	if( (*it_evtcls_flag).flag(0) >= 30 ){
	  H[10] -> accumulate( (float)Trkfit.helix(0) );
	  H[11] -> accumulate( (float)Trkfit.helix(3) );
	  H[12] -> accumulate( (float)P3.mag() );
	  H[13] -> accumulate( (float)P3.perp() );
	  H[14] -> accumulate( (float)P3.cosTheta() );
	  H[15] -> accumulate( (float)(P3.phi()*180.0/pi+180.) );
	}
      }
      if( MonitorLevel == 2 || MonitorLevel == 3 ){
	//	if( (*it_evtcls_flag).flag(0) >= 10 ){
	// new HadronA.
	if( (*it_hadronic_flag).hadronic_flag(1) > 0 ){
	  H[40] -> accumulate( (float)Trkfit.helix(0) );
	  H[41] -> accumulate( (float)Trkfit.helix(3) );
	  H[42] -> accumulate( (float)P3.mag() );
	  H[43] -> accumulate( (float)P3.perp() );
	  H[44] -> accumulate( (float)P3.cosTheta() );
	  H[45] -> accumulate( (float)(P3.phi()*180.0/pi+180.) );
	}
	//	if( (*it_evtcls_flag2).flag(3) >= 10 ){
	// HadronB.
	if( (*it_hadronic_flag).hadronic_flag(2) > 0 ){
	  H[70] -> accumulate( (float)Trkfit.helix(0) );
	  H[71] -> accumulate( (float)Trkfit.helix(3) );
	  H[72] -> accumulate( (float)P3.mag() );
	  H[73] -> accumulate( (float)P3.perp() );
	  H[74] -> accumulate( (float)P3.cosTheta() );
	  H[75] -> accumulate( (float)(P3.phi()*180.0/pi+180.) );
	}
      }

    }
    
    //*---
    // neutral clusters.
    //*---
    float Esum  = 0.0;
    float Ezsum = 0.0;
    for(Evtcls_hadron_neutral_Manager::iterator it = EvtHadClsMgr.begin();
	  it!=EvtHadClsMgr.end(); ++it ){

      Mdst_gamma& cls = (*it).gamma();
      if( cls ){
	Hep3Vector P3( cls.px(), cls.py(), cls.pz() );
	HepLorentzVector P4( P3, P3.mag() );
	HepLorentzVector P4Rest = P4;
	P4Rest.boost(CMboost);
	Hep3Vector P3Rest = P4Rest.vect();
	Esum  += P3Rest.mag();
	Ezsum += P3Rest.z();
	
	if( MonitorLevel == 1 || MonitorLevel == 3 ){
	  if( (*it_evtcls_flag).flag(0) >= 30 ){
	    H[16] -> accumulate( (float)P3.mag() );
	  }
	}
	if( MonitorLevel == 2 || MonitorLevel == 3 ){
	  //	  if( (*it_evtcls_flag).flag(0) >= 10 ){
	  // new HadronA.
	  if( (*it_hadronic_flag).hadronic_flag(1) > 0 ){
	    H[46] -> accumulate( (float)P3.mag() );
	  }
	  //	  if( (*it_evtcls_flag2).flag(3) >= 10 ){
	  // HadronB.
	  if( (*it_hadronic_flag).hadronic_flag(2) > 0 ){
	    H[76] -> accumulate( (float)P3.mag() );
	  }
	}
      }
    }
    //*---
    // ECL sum.
    //*---
    float ECLsum = 0.0;
    for(Mdst_ecl_Manager::iterator it = MdstEclMgr.begin();
	it != MdstEclMgr.end(); ++it ){
      
      Mdst_ecl& cls = *it;
      if( cls.quality() != 0 ) continue; // not associated to track
      if( cls.energy() < 0.1 ) continue; // low energy cut
      if( cls.theta()*180.0/pi < 17.0 ) continue; // ECL barrel cut
      if( cls.theta()*180.0/pi >150.0 ) continue;
      Hep3Vector P3( cls.energy()*sin( cls.theta() )*cos( cls.phi() ),
		     cls.energy()*sin( cls.theta() )*sin( cls.phi() ),
		     cls.energy()*cos( cls.theta() ) );
      HepLorentzVector P4( P3, P3.mag() );
      HepLorentzVector P4Rest = P4;
      P4Rest.boost(CMboost);
      Hep3Vector P3Rest = P4Rest.vect();
      ECLsum += P3Rest.mag();
    }  

    //*---
    // primary event vertex.
    //*---
    float PVr, PVz;
    bool PVok(false);
    if( Nvtx ){
      if( (*it_evtvtx).quality() >= 2 ){
	PVr = sqrt( (*it_evtvtx).PV_x()*(*it_evtvtx).PV_x()+
		    (*it_evtvtx).PV_y()*(*it_evtvtx).PV_y() );
	PVz = (*it_evtvtx).PV_z();
	PVok=true;
      }
    }

    //*---
    // global quantities.
    //*---
    float Evis = Psum + Esum;
    float Pz   = Pzsum + Ezsum;

    // HadronC
    if( MonitorLevel == 1 || MonitorLevel == 3 ){
      
      if( (*it_evtcls_flag).flag(0) >= 30 ){

	//trigger bits.      
	H[20] -> accumulate( (float) -1.0 );
	for( int i=0;i<48;i++ ){
	  if( bit[i] ) H[20] -> accumulate( (float)i );
	}      
	
	H[21] -> accumulate( (float)(*it_hadron_info).Ntrk() );
	H[22] -> accumulate( (float)(*it_hadron_info).Ncls() );
	H[23] -> accumulate( (float)Psum );
	H[24] -> accumulate( (float)ECLsum );
	H[25] -> accumulate( (float)Esum );
	H[26] -> accumulate( (float)Evis );
	H[27] -> accumulate( (float)Pz );
	H[28] -> accumulate( (float)(*it_hadron_info).R2() );
	if( (*it_hadron_info).R2()<0.96 ) 
	  H[29] -> accumulate( (float)(*it_hadron_info).R2() );
	if(PVok) {
	  H[30] -> accumulate( (float)PVr );
	  H[31] -> accumulate( (float)PVz );
	}
      }
    }

    //HadronA
    if( MonitorLevel == 2 || MonitorLevel == 3 ){

      //      if( (*it_evtcls_flag).flag(0) >= 10 ){
      // new HadronA.
      if( (*it_hadronic_flag).hadronic_flag(1) > 0 ){

	//trigger bits.      
	H[50] -> accumulate( (float) -1.0 );
	for( int i=0;i<48;i++ ){
	  if( bit[i] ) H[50] -> accumulate( (float)i );
	}      
	
	H[51] -> accumulate( (float)(*it_hadron_info).Ntrk() );
	H[52] -> accumulate( (float)(*it_hadron_info).Ncls() );
	H[53] -> accumulate( (float)Psum );
	H[54] -> accumulate( (float)ECLsum );
	H[55] -> accumulate( (float)Esum );
	H[56] -> accumulate( (float)Evis );
	H[57] -> accumulate( (float)Pz );
	H[58] -> accumulate( (float)(*it_hadron_info).R2() );
	if( (*it_hadron_info).R2()<0.96 ) 
	  H[59] -> accumulate( (float)(*it_hadron_info).R2() );
	if (PVok) {
	  H[60] -> accumulate( (float)PVr );
	  H[61] -> accumulate( (float)PVz );
	}
      }

      //      if( (*it_evtcls_flag2).flag(3) >= 10 ){
      // HadronB.
      if( (*it_hadronic_flag).hadronic_flag(2) > 0 ){

	//trigger bits.      
	H[80] -> accumulate( (float) -1.0 );
	for( int i=0;i<48;i++ ){
	  if( bit[i] ) H[80] -> accumulate( (float)i );
	}      
	
	H[81] -> accumulate( (float)(*it_hadron_info).Ntrk() );
	H[82] -> accumulate( (float)(*it_hadron_info).Ncls() );
	H[83] -> accumulate( (float)Psum );
	H[84] -> accumulate( (float)ECLsum );
	H[85] -> accumulate( (float)Esum );
	H[86] -> accumulate( (float)Evis );
	H[87] -> accumulate( (float)Pz );
	H[88] -> accumulate( (float)(*it_hadron_info).R2() );
	if( (*it_hadron_info).R2()<0.96 ) 
	  H[89] -> accumulate( (float)(*it_hadron_info).R2() );
	if(PVok) {
	  H[90] -> accumulate( (float)PVr );
	  H[91] -> accumulate( (float)PVz );
	}
      }

    }

  }

      


}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
