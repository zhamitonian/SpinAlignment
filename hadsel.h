//
// $Id: hadsel.h 10002 2007-02-26 06:56:17Z katayama $
//
// $Log$
// Revision 1.2  2002/01/03 11:04:32  katayama
// Point3D and other header files are cleaned
//
// Revision 1.1  2001/06/20 12:39:31  katayama
// New from Garmash san
//
//
// This is my hadronic selection function

// The imputs are the panther table managers 
// for the mdst_charged table and the mdst_ecl table 
// and the CMS energy

//the output is a flag = 0 for fail; 1 for pass the cuts
// The next 4 are the 4 variables I use, Esum, Evis, nECL and HJM.


// If you didn't already include these
#include "belle.h"
#include <list>
#include "belleCLHEP/Vector/ThreeVector.h"
#include "belleCLHEP/Vector/LorentzVector.h"
#include "toolbox/Thrust.h"
#include "toolbox/FuncPtr.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif



struct HadB{

  int type;     // = 0 for pass, >0 for fail
  float Esum;   // cluster energy sum
  float Evis;   // sum of track and gamma energy
  int nECL;     // cluster multiplicity
  float HJM;    // invariant mass of one hemisphere (jet)

};

// here, the variables used to apply the hadronic cuts has difference with var in evtcls_hadron_info 
// event select difference is about 1/10000
HadB good_hadronB(Mdst_charged_Manager& ChgMgr, Mdst_ecl_Manager& EclMgr, float Ecm){

  float Esum = 0;  
  float Evis = 0;  
  int nECL = 0;    
  float HJM = 0;   

  std::vector<Hep3Vector> GoodParticles3;
  std::vector<HepLorentzVector> GoodParticles4;

  for(std::vector<Mdst_charged>::iterator it = ChgMgr.begin();
      
      it != ChgMgr.end(); ++it){

    Mdst_charged& Charged = *it;  

    float dz = 1000;
    float dr = 1000;

    Mdst_trk& Trk = Charged.trk();
    Mdst_trk_fit& Fit = Trk.mhyp(2);  // assume pion mass 

    if(Trk){
      if(Fit){

	dz = Fit.helix(3);  //get the impact parameter
	dr = Fit.helix(0);  // from the nominal IP

      }
    }

    Hep3Vector p3lab ( Charged.px(), Charged.py(), Charged.pz() ); 

    if(p3lab.mag()>0 && p3lab.mag()<20){  //remove crazy numbers
            
      float EPi =  sqrt(p3lab.mag2() + 0.0195);  //assume pion mass
      
      HepLorentzVector p4cm (p3lab.x(), p3lab.y(), p3lab.z(), EPi);
      
      p4cm.boost(-0.0153, 0.0, -0.3909);  //belle boost

      Hep3Vector p3cm (p4cm.x(), p4cm.y(), p4cm.z());

      if(p3lab.mag()>0.1){  // 'good track' cuts
	if(abs(dr)<2){
	  if(abs(dz)<4){
    
	    Evis = Evis + p4cm.t();  // Sum the track energies
	    
	    GoodParticles3.push_back(p3cm);  //add tracks to 
	    GoodParticles4.push_back(p4cm);  //the particle list
	  }
	}
      }
    }
  }

  for(std::vector<Mdst_ecl>::iterator icl = EclMgr.begin();
      
      icl != EclMgr.end(); ++icl){
    
    Mdst_ecl& Clus = *icl;

    int match = Clus.match();
    int quality = Clus.quality();

    float energy = Clus.energy();
    float theta =  Clus.theta();
    float phi =  Clus.phi();

    float Ex = energy*sin(theta)*cos(phi);
    float Ey = energy*sin(theta)*sin(phi);
    float Ez = energy*cos(theta);
    
    HepLorentzVector Clus4cm (Ex, Ey, Ez, energy);   

    Clus4cm.boost(-0.0153, 0.0, -0.3909);

    Hep3Vector Clus3cm (Clus4cm.x(), Clus4cm.y(), Clus4cm.z());
    
    if(energy>0.1){  // 'good cluster' cuts
      if(quality==0){
	
	Esum = Esum + Clus4cm.t();  // sum the cluster energies

	nECL++;   // sum up the clusters

	if(match==0){  //'good unmatched cluster'
	
	  Evis = Evis + Clus4cm.t();  //add gamma energies 
	                              // with track energies

	  GoodParticles3.push_back(Clus3cm);  //add gammas 
	                                      // to the particle list
	  GoodParticles4.push_back(Clus4cm);  //add gammas to the 
	                                      //particle list

	}
      }
    }
  }

  Hep3Vector thr = thrust( GoodParticles3.begin(), //calculate the
			GoodParticles3.end(),   //thrust vector using
			SelfFunc(Hep3Vector()));   // Nakao-san's function

  Hep3Vector tn  = thr.unit();  // thrust axis 

  HepLorentzVector Vpos (0.0, 0.0, 0.0, 0.0);
  HepLorentzVector Vneg (0.0, 0.0, 0.0, 0.0);

  for(std::vector<HepLorentzVector>::iterator ip4 = GoodParticles4.begin();
      
      ip4 != GoodParticles4.end(); ++ip4){
    
    HepLorentzVector& p4 = *ip4;

    Hep3Vector p3 (p4.x(), p4.y(), p4.z());

    double scalar = p3*tn;
 
    if( scalar >= 0.0 ){  // cos theta > 0 = one hemisphere

      Vpos = Vpos + p4;

    }
    else{  // cos theta < 0 = other hemisphere
      
      Vneg = Vneg + p4;
       
    }
  }
    
  if( Vpos.mag() >= Vneg.mag() ){

    HJM = Vpos.mag();
  
  }
  else{
  
    HJM = Vneg.mag();
    
  }

  // now for the cuts

  int type = 3;  // passes no cuts;

  // cut1: 
  //  Tight Esum cut to get rid of the low energy 
  // two photon, tau's, and beam gas events.
  // To recover some of the continuum efficiency loss
  // make an 'or' with the heavy jet mass > 2 GeV

  if(Esum/Ecm > 0.18 || HJM > 1.8){

    type = 2; //passes the first cut.

    //cut2:
    // heavy jet mass cut
    // For hadronic events, the heavy jet mass is about 1/2 the 
    // visible energy.
    //  For tau's, the heavy jet mass is always below 2 GeV
    // For QED, the heavy jet mass is almost always below 2 GeV
    //  This cut also gets rid of a few beam gas and 2 photon events
    // To recover some of the continuum efficiency loss
    // make an 'or' with the heavy jet mass > 1.8 GeV

    if(HJM > 0.25*Evis || HJM > 1.8){

      type = 1;  //passes second cut.

      //cut3:
      // average cluster energy
      // For hadronic events, the average cluster energy is 
      // usually always below 1 GeV
      // For the QED events still in the skim, 
      // its always above 1 GeV

      if(nECL>0){
	if(Esum/nECL < 1){

	  type = 0;

	}
      }
    }
  }

  HadB out;
  
  out.type = type;
  out.Esum = Esum;
  out.Evis = Evis;
  out.nECL = nECL;
  out.HJM = HJM;

  return out;

}

struct NewHadA{

  int type;     // = 0 for pass, >0 for fail
  float Esum;   // cluster energy sum
  int nECLfvc;   // cluster multiplicity in -0.7 - 0.9 cos theta

};

NewHadA new_hadronA(Mdst_ecl_Manager& EclMgr, float Ecm){

  int type = 2;
  float Esum = 0;
  int nECLfvc = 0;

  float thetamin1 = 17*3.1416/180;
  float thetamax1 = 150*3.1416/180;

  float costhetamin2 = -0.7;
  float costhetamax2 = 0.9;

  for(std::vector<Mdst_ecl>::iterator icl = EclMgr.begin();
      
      icl != EclMgr.end(); ++icl){
    
    Mdst_ecl& Clus = *icl;

    int match = Clus.match();
    int quality = Clus.quality();
    
    float energy = Clus.energy();
    float theta =  Clus.theta();
    float phi =  Clus.phi();

    float Ex = energy*sin(theta)*cos(phi);
    float Ey = energy*sin(theta)*sin(phi);
    float Ez = energy*cos(theta);
    
    HepLorentzVector Clus4cm (Ex, Ey, Ez, energy);   

    Clus4cm.boost(-0.0153, 0.0, -0.3909);

    Hep3Vector Clus3cm (Clus4cm.x(), Clus4cm.y(), Clus4cm.z());
    
    if(energy>0.1){  // 'good cluster' cuts
      if(quality==0){

	if(theta>thetamin1 && theta<thetamax1){
		
	  Esum = Esum + Clus4cm.t();  // sum the cluster energies

	}

	if(cos(theta)>costhetamin2 && 
	   cos(theta)<costhetamax2){

	  nECLfvc++;   // sum up the clusters

	}
      }
    }
  }
  
  float nEsum = Esum/Ecm;

  //cut1:  Esum cut:

  if(nEsum > 0.1 && nEsum < 0.8){

    type = 1;

    //cut2:  cluster multiplicity 
    //       in a feducial volume

    if(nECLfvc > 1){

      type = 0;

    }
  }
  
  NewHadA out;

  out.type = type;
  out.Esum = Esum;
  out.nECLfvc = nECLfvc;

  return out;

}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
