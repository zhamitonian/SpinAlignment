#include "belle.h"
#include "panther/panther.h"
#include "particle/utility.h"
#include "particle/combination.h"
#include "ip/IpProfile.h"
#include "tables/mdst.h"
#include "userinfo.h"
#include "myutils.h"
#include "user_def.h"
#include "kfit.h"

#include "kmassfitterExt.h"
#include "kmassfitterExt.cc"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

/////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////           
// Copy from Li Jin     
// Calculates the error on LifeTime calculated from Decay - Production Vertex   
// Output: errT   - 3D           lifetime error  
//         errTxy - 2D(xy-plane) lifetime error 
//         errTx  -    
//         errTy  - 1D           lifetime errors 
//         errTz  -  
// Input:  dec    - decayVertex  
//         ip     - productionVertex   
//         p      - momentum 
//         errIP  - error matrix of Production vertex    
//         errP   - error matrix of momentum p    
//         covIPP - correltation matrix between IP & P   
void errLifetime(double &errT, double &errTxy, double &errTx, double &errTy, double &errTz,
		 const HepPoint3D &dec, const HepPoint3D &ip, const HepLorentzVector &p,
		 const HepSymMatrix &errIp, const HepSymMatrix &errP, const HepMatrix &covIPP,
		 int flag) {
  double f   = 10.*Ptype("D0").mass()/CV*10000.0/(p.vect().mag2());
  double fxy = 10.*Ptype("D0").mass()/CV*10000.0/(p.vect().perp2());
  double fx  = 10.*Ptype("D0").mass()/CV*10000.0/(p.x()*p.x());
  double fy  = 10.*Ptype("D0").mass()/CV*10000.0/(p.y()*p.y());
  double fz  = 10.*Ptype("D0").mass()/CV*10000.0/(p.z()*p.z());
  // 
  //     |  Dec  : D-IP :  D-M  |                                                                                            
  //     |                      |                                                                                 
  // V = | IP-D  :  IP  : IP-M  |  
  //     |                      | 
  //     |  M-D  : M-IP :  Mom  | 
  //
  HepSymMatrix V(9,0);
  V.sub(1, errP.sub(5,7)); // DEC                                                                                                             
  V.sub(4, errIp);         // IP                                                                                                              
  V.sub(7, errP.sub(1,3)); // Mom                                                                                                             
  if(flag){
    // set correlation                                                                                                                       
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

  // xyz                                                                                                                                     
  HepMatrix dtdV(9,1,0);
  dtdV[0][0] =  f*p.x();
  dtdV[1][0] =  f*p.y();
  dtdV[2][0] =  f*p.z();
  dtdV[3][0] = -f*p.x();
  dtdV[4][0] = -f*p.y();
  dtdV[5][0] = -f*p.z();
  dtdV[6][0] = f*(dec-ip).x()-2.*f*((dec-ip)*(p.vect()))*p.x()/p.vect().mag2();
  dtdV[7][0] = f*(dec-ip).y()-2.*f*((dec-ip)*(p.vect()))*p.y()/p.vect().mag2();
  dtdV[8][0] = f*(dec-ip).z()-2.*f*((dec-ip)*(p.vect()))*p.z()/p.vect().mag2();
  errT = ((dtdV.T())*V*dtdV)[0][0];

  // xy                                                                                                                                      
  Hep3Vector momP = p.vect();
  momP.setZ(0.);
  HepMatrix dtxydV(9,1,0);
  dtxydV[0][0] =  fxy*p.x();
  dtxydV[1][0] =  fxy*p.y();
  dtxydV[2][0] =  0.; // fxy*p.z();                                                                                                          
  dtxydV[3][0] = -fxy*p.x();
  dtxydV[4][0] = -fxy*p.y();
  dtxydV[5][0] =  0.; // -fxy*p.z();                                                                                                         
  dtxydV[6][0] =  fxy*(dec-ip).x()-2.*fxy*((dec-ip)*momP)*p.x()/p.vect().perp2();
  dtxydV[7][0] =  fxy*(dec-ip).y()-2.*fxy*((dec-ip)*momP)*p.y()/p.vect().perp2();
  dtxydV[8][0] =  0.; // fxy*(dec-ip).z()-2.*fxy*((dec-ip)*momP)*p.z()/p.vect().perp2();                                                     
  errTxy = ((dtxydV.T())*V*dtxydV)[0][0];

  // x                                                                                                                                       
  HepMatrix dtxdV(9,1,0);
  dtxdV[0][0] =  fx*p.x();
  dtxdV[1][0] =  0.; // fx*p.y();                                                                                                            
  dtxdV[2][0] =  0.; // fx*p.z();                                                                                                            
  dtxdV[3][0] = -fx*p.x();
  dtxdV[4][0] =  0.; // -fx*p.y();                                                                                                           
  dtxdV[5][0] =  0.; // -fx*p.z();                                                                                                           
  dtxdV[6][0] = -fx*(dec-ip).x(); // -2.*fx*((dec-ip)*(p.vect()))/p.x();                                                                     
  dtxdV[7][0] =  0.; //fx*(dec-ip).y()-2.*fx*((dec-ip)*(p.vect()))*p.y()/(p.x()*p.x());                                                      
  dtxdV[8][0] =  0.; //fx*(dec-ip).z()-2.*fx*((dec-ip)*(p.vect()))*p.z()/(p.x()*p.x());                                                      
  errTx = ((dtxdV.T())*V*dtxdV)[0][0];

  // y                                                                                                                                       
  HepMatrix dtydV(9,1,0);
  dtydV[0][0] =  0.; // fy*p.x();                                                                                                            
  dtydV[1][0] =  fy*p.y();
  dtydV[2][0] =  0.; // fy*p.z();                                                                                                            
  dtydV[3][0] =  0.; // -fy*p.x();                                                                                                           
  dtydV[4][0] = -fy*p.y();
  dtydV[5][0] =  0.; // -fy*p.z();                                                                                                           
  dtydV[6][0] =  0.; // fy*(dec-ip).x()-2.*fy*((dec-ip)*(p.vect()))*p.x()/(p.y()*p.y());                                                     
  dtydV[7][0] = -fy*(dec-ip).y(); //-2.*fy*((dec-ip)*(p.vect()))/p.y();                                                                      
  dtydV[8][0] =  0.; //fy*(dec-ip).z()-2.*fy*((dec-ip)*(p.vect()))*p.z()/(p.y()*p.y());                                                      
  errTy = ((dtydV.T())*V*dtydV)[0][0];

  // z                                                                                                                                       
  HepMatrix dtzdV(9,1,0);
  dtzdV[0][0] =  0.; // fz*p.x();                                                                                                            
  dtzdV[1][0] =  0.; // fz*p.y();                                                                                                            
  dtzdV[2][0] =  fz*p.z();
  dtzdV[3][0] =  0.; // -fz*p.x();                                                                                                           
  dtzdV[4][0] =  0.; // -fz*p.y();                                                                                                           
  dtzdV[5][0] = -fz*p.z();
  dtzdV[6][0] =  0.; //fz*(dec-ip).x()-2.*fz*((dec-ip)*(p.vect()))*p.x()/(p.z()*p.z());                                                      
  dtzdV[7][0] =  0.; //fz*(dec-ip).y()-2.*fz*((dec-ip)*(p.vect()))*p.y()/(p.z()*p.z());                                                      
  dtzdV[8][0] = -fz*(dec-ip).z(); //-2.*fz*((dec-ip)*(p.vect()))/p.z();                                                                      
  errTz = ((dtzdV.T())*V*dtzdV)[0][0];

  if(fabs(errT)   > 1.e+30)errT   = 1.e+30;
  if(fabs(errTxy) > 1.e+30)errTxy = 1.e+30;
  if(fabs(errTx)  > 1.e+30)errTx  = 1.e+30;
  if(fabs(errTy)  > 1.e+30)errTy  = 1.e+30;
  if(fabs(errTz)  > 1.e+30)errTz  = 1.e+30;
}


  // Mass constrained fit
  
  int doKmFit(Particle &p) {
    if(! &p.userInfo() ) setUserInfo(p);

    dynamic_cast<UserInfo&>(p.userInfo()).setMassBeforeVertexFit(p.mass());

    // Do Mass fit
    kmassfitter kmf;
    int nfit=0;
    for(unsigned i=0; i<p.nChildren(); i++) {
      if(p.child(i).nChildren()<2){
      addTrack2fit(kmf,p.child(i));
      nfit++;
      }
       else {
        for(unsigned j=0; j<p.child(i).nChildren(); j++) {
                addTrack2fit(kmf,p.child(i).child(j));
                nfit++;
          }
      }
    }
    kmf.invariantMass(p.pType().mass());
//    kmf.notDecayPoint();
    unsigned err = kmf.fit();
    p.usable(!err);

//    if(err == 0) {
      // save fitted Momentum and Dalitz Variables
//      kmakemother kmm;
//      int ok =  makeMother(kmm,kmf,p,0);
//      p.usable(ok);
//      return  ok;
//    }
//   else return 0;
    return 0;
  }
  
  void doKmFit(std::vector<Particle> &plist) {
    for(unsigned i=0;i<plist.size();++i) 
      doKmFit(plist[i]);
    
    rm_unusable( plist );
  }
  

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// Vertex FIT followed by MASS fit
// momentum of mother particle (In Particle object) is updated after vertex fit
// momentum of mother particle after mass fit is saved to UserInfo class
// chi2_cut -> cut on chi2/ndf

void doKvmFit(Particle &p, float chi2_cut) {
  if(! &p.userInfo() ) setUserInfo(p);

  dynamic_cast<UserInfo&>(p.userInfo()).setMassBeforeVertexFit(p.mass());

  kvertexfitter kv;

  int nfit=0;
  for(unsigned i=0; i<p.nChildren(); i++) {
    addTrack2fit(kv,p.child(i));
    nfit++;
  }
  
  // Do vertex fit
  int err=kv.fit();
  
  // set usable flag; if fit is OK
  p.usable(!err);
  if(!err) {
    if(kv.chisq()<0) {
      p.usable(0);
      return;
    }
    
    dynamic_cast<UserInfo&>(p.userInfo()).setType(p.lund());
    dynamic_cast<UserInfo&>(p.userInfo()).setVChiSq(kv.chisq());
    dynamic_cast<UserInfo&>(p.userInfo()).setNDF(kv.dgf());
    p.momentum().decayVertex(kv.vertex(), kv.errVertex());
    
    float chi = kv.chisq();
    int ndf = kv.dgf();
    float prob = prob_chi2( chi, ndf );
    p.usable(chi/ndf < chi2_cut);
    if(p.usable()) {
      dynamic_cast<UserInfo&>(p.userInfo()).setVProb(prob);
      int ok=makeMother(kv, p);
      p.usable(ok);
      
      // Do Mass fit
      kmassfitter kmf;
      nfit=0;
      for(unsigned i=0; i<p.nChildren(); i++) {
	addTrack2fit(kmf,p.child(i));
	nfit++;
      }
      kmf.invariantMass(p.pType().mass());
      kmf.vertex(kv.vertex());
      kmf.errVertex(kv.errVertex());
      err = kmf.fit();
      p.usable(!err);
      // save fitted Momentum and Dalitz Variables
      ok = saveFittedMomentumAfterMassFit(p, kmf);
//      saveDalitzPlotVariables(p, kmf);
      dynamic_cast<UserInfo&>(p.userInfo()).setvtex(kmf.vertex());   
      p.usable(ok);

    }
  } else {
    //std::cout << "Error doKvFit: " << err << std::endl;
  }
}

void doKvmFit(std::vector<Particle> &plist, float chi2_cut) {
  for(unsigned i=0;i<plist.size();++i) 
    doKvmFit(plist[i], chi2_cut);

  rm_unusable( plist );
}

void doKvmFit(Particle &p, float chi2_cut, double MASS) {
	if(! &p.userInfo() ) setUserInfo(p);

	dynamic_cast<UserInfo&>(p.userInfo()).setMassBeforeVertexFit(p.mass());

	kvertexfitter kv;

	int nfit=0;
	for(unsigned i=0; i<p.nChildren(); i++) {
		addTrack2fit(kv,p.child(i));
		nfit++;
	}

	// Do vertex fit
	int err=kv.fit();

	// set usable flag; if fit is OK
	p.usable(!err);
	if(!err) {
		if(kv.chisq()<0) {
			p.usable(0);
			return;
		}

		dynamic_cast<UserInfo&>(p.userInfo()).setType(p.lund());
		dynamic_cast<UserInfo&>(p.userInfo()).setVChiSq(kv.chisq());
		dynamic_cast<UserInfo&>(p.userInfo()).setNDF(kv.dgf());
		p.momentum().decayVertex(kv.vertex(), kv.errVertex());

		float chi = kv.chisq();
		int ndf = kv.dgf();
		float prob = prob_chi2( chi, ndf );
		p.usable(chi/ndf < chi2_cut);
		if(p.usable()) {
			dynamic_cast<UserInfo&>(p.userInfo()).setVProb(prob);
			int ok=makeMother(kv, p);
			p.usable(ok);

			// Do Mass fit
			kmassfitter kmf;
			nfit=0;
			for(unsigned i=0; i<p.nChildren(); i++) {
				addTrack2fit(kmf,p.child(i));
				nfit++;
			}
			kmf.invariantMass(MASS);
			kmf.vertex(kv.vertex());
			kmf.errVertex(kv.errVertex());
			err = kmf.fit();
			p.usable(!err);
			// save fitted Momentum and Dalitz Variables
			ok = saveFittedMomentumAfterMassFit(p, kmf);
			saveDalitzPlotVariables(p, kmf);

			p.usable(ok);

		}
	} else {
		//std::cout << "Error doKvFit: " << err << std::endl;
	}
}

void doKvmFit(std::vector<Particle> &plist, float chi2_cut, double MASS) {
	for(unsigned i=0;i<plist.size();++i) doKvmFit(plist[i], chi2_cut, MASS);
	rm_unusable( plist );
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// Vertex FIT followed by MASS fit for candidates with Pi0 child
// momentum of mother particle (In Particle object) is updated after vertex fit
// momentum of mother particle after mass fit is saved to UserInfo class
// chi2_cut -> cut on chi2/ndf

void doKvmFitWPi0(Particle &p, float chi2_cut) {
	if(! &p.userInfo() ) setUserInfo(p);

	dynamic_cast<UserInfo&>(p.userInfo()).setMassBeforeVertexFit(p.mass());

	kvertexfitter kv;

	int nfit=0;
	int pi0 = -1;
	for(unsigned i=0; i<p.nChildren(); i++) {
		if(p.child(i).lund()!=111) {
			addTrack2fit(kv,p.child(i));
			nfit++;
		} else {
			pi0 = i;
		}
	}
	// Do vertex fit (with only charged tracks)
	int err=kv.fit();

	// set usable flag; if fit is OK
	p.usable(!err);
	if(!err) {
		// mass fit of PI0
		Particle g0 = p.child(pi0).child(0);
		Particle g1 = p.child(pi0).child(1);
		setGammaError(g0, kv.vertex(), kv.errVertex());
		setGammaError(g1, kv.vertex(), kv.errVertex());

		kmassfitter kmf_pi0;
		addTrack2fit(kmf_pi0, g0);
		addTrack2fit(kmf_pi0, g1);
		kmf_pi0.invariantMass(p.child(pi0).pType().mass());
		kmf_pi0.vertex(kv.vertex());
		kmf_pi0.errVertex(kv.errVertex());
		err = kmf_pi0.fit();

		p.usable(!err);
		if(!err) {
			Particle pi0_part = p.child(pi0);
			kmakemother kmm_pi0;
			makeMother(kmm_pi0,kmf_pi0,pi0_part,0);

			// redo the vertex fit (now include pi0)
			kvertexfitter kv2;

			nfit=0;
			for(unsigned i=0; i<p.nChildren(); i++) {
				if(p.child(i).lund()!=111) {
					addTrack2fit(kv2,p.child(i));
					nfit++;
				} else {
					addTrack2fit(kv2,pi0_part);
					nfit++;
				}
			}

			err = kv2.fit();
			p.usable(!err);
			if(!err) {
				if(kv2.chisq()<0) {
					p.usable(0);
					return;
				}
				dynamic_cast<UserInfo&>(p.userInfo()).setType(p.lund());
				dynamic_cast<UserInfo&>(p.userInfo()).setVChiSq(kv2.chisq());
				dynamic_cast<UserInfo&>(p.userInfo()).setNDF(kv2.dgf());
				p.momentum().decayVertex(kv2.vertex(), kv2.errVertex());

				float chi = kv2.chisq();
				int ndf = kv2.dgf();
				float prob = prob_chi2( chi, ndf );
				p.usable(chi/ndf < chi2_cut);

				if(p.usable()) {
					dynamic_cast<UserInfo&>(p.userInfo()).setVProb(prob);
					int ok=makeMother(kv2, p);
					p.usable(ok);

					// Do Mass fit
					kmassfitter kmf;
					nfit=0;
					for(unsigned i=0; i<p.nChildren(); i++) {
						if(p.child(i).lund()!=111) {
							addTrack2fit(kmf,p.child(i));
							nfit++;
						} else {
							addTrack2fit(kmf,pi0_part);
							nfit++;
						}
					}
					kmf.invariantMass(p.pType().mass());
					kmf.vertex(kv2.vertex());
					kmf.errVertex(kv2.errVertex());
					err = kmf.fit();
					p.usable(!err);
					// save fitted Momentum and Dalitz Variables
					ok = saveFittedMomentumAfterMassFit(p, kmf);
					saveDalitzPlotVariables(p, kmf);

					p.usable(ok);
				}

			}
		}
	} else {
		//std::cout << "Error doKvFit: " << err << std::endl;
	}
}
void doKvmFitWPi0(std::vector<Particle> &plist, float chi2_cut) {
	for(unsigned i=0;i<plist.size();++i) doKvmFitWPi0(plist[i], chi2_cut);
	rm_unusable( plist );
}

void doKvmFitWPi0(Particle &p, float chi2_cut, double MASS) {
	if(! &p.userInfo() ) setUserInfo(p);

	dynamic_cast<UserInfo&>(p.userInfo()).setMassBeforeVertexFit(p.mass());

	kvertexfitter kv;

	int nfit=0;
	int pi0 = -1;
	for(unsigned i=0; i<p.nChildren(); i++) {
		if(p.child(i).lund()!=111) {
			addTrack2fit(kv,p.child(i));
			nfit++;
		} else {
			pi0 = i;
		}
	}
	// Do vertex fit (with only charged tracks)
	int err=kv.fit();

	// set usable flag; if fit is OK
	p.usable(!err);
	if(!err) {
		// mass fit of PI0
		Particle g0 = p.child(pi0).child(0);
		Particle g1 = p.child(pi0).child(1);
		setGammaError(g0, kv.vertex(), kv.errVertex());
		setGammaError(g1, kv.vertex(), kv.errVertex());

		kmassfitter kmf_pi0;
		addTrack2fit(kmf_pi0, g0);
		addTrack2fit(kmf_pi0, g1);
		kmf_pi0.invariantMass(p.child(pi0).pType().mass());
		kmf_pi0.vertex(kv.vertex());
		kmf_pi0.errVertex(kv.errVertex());
		err = kmf_pi0.fit();

		p.usable(!err);
		if(!err) {
			Particle pi0_part = p.child(pi0);
			kmakemother kmm_pi0;
			makeMother(kmm_pi0,kmf_pi0,pi0_part,0);

			// redo the vertex fit (now include pi0)
			kvertexfitter kv2;

			nfit=0;
			for(unsigned i=0; i<p.nChildren(); i++) {
				if(p.child(i).lund()!=111) {
					addTrack2fit(kv2,p.child(i));
					nfit++;
				} else {
					addTrack2fit(kv2,pi0_part);
					nfit++;
				}
			}

			err = kv2.fit();
			p.usable(!err);
			if(!err) {
				if(kv2.chisq()<0) {
					p.usable(0);
					return;
				}
				dynamic_cast<UserInfo&>(p.userInfo()).setType(p.lund());
				dynamic_cast<UserInfo&>(p.userInfo()).setVChiSq(kv2.chisq());
				dynamic_cast<UserInfo&>(p.userInfo()).setNDF(kv2.dgf());
				p.momentum().decayVertex(kv2.vertex(), kv2.errVertex());

				float chi = kv2.chisq();
				int ndf = kv2.dgf();
				float prob = prob_chi2( chi, ndf );
				p.usable(chi/ndf < chi2_cut);

				if(p.usable()) {
					dynamic_cast<UserInfo&>(p.userInfo()).setVProb(prob);
					int ok=makeMother(kv2, p);
					p.usable(ok);

					// Do Mass fit
					kmassfitter kmf;
					nfit=0;
					for(unsigned i=0; i<p.nChildren(); i++) {
						if(p.child(i).lund()!=111) {
							addTrack2fit(kmf,p.child(i));
							nfit++;
						} else {
							addTrack2fit(kmf,pi0_part);
							nfit++;
						}
					}
					kmf.invariantMass(MASS);
					kmf.vertex(kv2.vertex());
					kmf.errVertex(kv2.errVertex());
					err = kmf.fit();
					p.usable(!err);
					// save fitted Momentum and Dalitz Variables
					ok = saveFittedMomentumAfterMassFit(p, kmf);
					saveDalitzPlotVariables(p, kmf);

					p.usable(ok);
				}

			}
		}
	} else {
		//std::cout << "Error doKvFit: " << err << std::endl;
	}
}
void doKvmFitWPi0(std::vector<Particle> &plist, float chi2_cut, double MASS) {
	for(unsigned i=0;i<plist.size();++i) doKvmFitWPi0(plist[i], chi2_cut, MASS);
	rm_unusable( plist );
}

unsigned saveFittedMomentumAfterMassFit(Particle &mother, kmassfitter &km) {

	kmakemother kmm;
	unsigned n = km.tracks();

	for(unsigned i=0;i<n;++i){
		kmm.addTrack(km.momentum(i),
				km.position(i),
				km.error(i),
				mother.relation().child(i).pType().charge());

		if(km.fitWithVertex())
			kmm.errVertexTrack(km.errVertexTrack(i));

		for(unsigned j=i+1;j<n;++j){
			kmm.correlation(km.correlation(i,j));
		}
	}

	kmm.vertex(km.vertex());

	if(km.fitWithVertex()){
		kmm.errVertex(km.errVertex());
	}

	unsigned err = kmm.make();
	if(err != 0)
		return 0;

	// save Momentum (HepLorentzVector)
	//dynamic_cast<UserInfo&>(mother.userInfo()).setFittedMomentum(kmm.momentum());
	dynamic_cast<UserInfo&>(mother.userInfo()).getFittedMomentum().momentumPosition(kmm.momentum(),kmm.position(),kmm.error());
	dynamic_cast<UserInfo&>(mother.userInfo()).getFittedMomentum().decayVertex(km.vertex(), km.errVertex());

	return 1;
}

void saveDalitzPlotVariables(Particle &mother, kmassfitter &km) {

	unsigned n = km.tracks();

	if(n==2)
		return;

	if(n==3) {
		dynamic_cast<UserInfo&>(mother.userInfo()).setM12sq((km.momentum(0)+km.momentum(1)).mag2());
		dynamic_cast<UserInfo&>(mother.userInfo()).setM23sq((km.momentum(1)+km.momentum(2)).mag2());
		dynamic_cast<UserInfo&>(mother.userInfo()).setM13sq((km.momentum(0)+km.momentum(2)).mag2());
	}

	if(n==4) {
		dynamic_cast<UserInfo&>(mother.userInfo()).setM12sq((km.momentum(0)+km.momentum(1)).mag2());
		dynamic_cast<UserInfo&>(mother.userInfo()).setM23sq((km.momentum(1)+km.momentum(2)).mag2());
		dynamic_cast<UserInfo&>(mother.userInfo()).setM13sq((km.momentum(0)+km.momentum(2)).mag2());
		dynamic_cast<UserInfo&>(mother.userInfo()).setM14sq((km.momentum(0)+km.momentum(3)).mag2());
		dynamic_cast<UserInfo&>(mother.userInfo()).setM24sq((km.momentum(1)+km.momentum(3)).mag2());
		dynamic_cast<UserInfo&>(mother.userInfo()).setM34sq((km.momentum(2)+km.momentum(3)).mag2());

		dynamic_cast<UserInfo&>(mother.userInfo()).setM123sq((km.momentum(0)+km.momentum(1)+km.momentum(2)).mag2());
		dynamic_cast<UserInfo&>(mother.userInfo()).setM124sq((km.momentum(0)+km.momentum(1)+km.momentum(3)).mag2());
		dynamic_cast<UserInfo&>(mother.userInfo()).setM134sq((km.momentum(0)+km.momentum(2)+km.momentum(3)).mag2());
		dynamic_cast<UserInfo&>(mother.userInfo()).setM234sq((km.momentum(1)+km.momentum(2)+km.momentum(3)).mag2());
	}
}

void refitDKXSystem(Particle &p, HepLorentzVector BEAM, int withXFrag, int withProton) {
	if(! &p.userInfo() ) setUserInfo(p);

	dynamic_cast<UserInfo&>(p.userInfo()).setMassBeforeVertexFit(p.mass());
	std::vector<Particle> neutralPi0s;
	std::vector<Particle> fitOne;
	std::vector<Particle> fitAll;

	//std::cout << " ****** refitDKXSystem ********* " << std::endl;
	//std::cout << " Dtag  = " << dynamic_cast<UserInfo&>(p.child(0).userInfo()).decayMode() << std::endl;
	//std::cout << " Kaon  = " << p.child(1).lund() << std::endl;


	kvertexfitter kv;

	// add IP profile to vertex fit
	HepPoint3D bp = IpProfile::position(EVENT_DEPEND_IP);
	HepSymMatrix berr = IpProfile::position_err(EVENT_DEPEND_IP);
	kv.initialVertex(bp);
	kv.vertexProfile(berr);

	// add Dtag (mass vertex fitted momentum)
	int ddm = dynamic_cast<UserInfo&>(p.child(0).userInfo()).decayMode();
	int gamma_from_DTagStar = 0;
	if(ddm<100) { // D0 or D+
		Particle Dtag = p.child(0);
		kv.addTrack(dynamic_cast<UserInfo&>(Dtag.userInfo()).getFittedMomentum().p(),
				dynamic_cast<UserInfo&>(Dtag.userInfo()).getFittedMomentum().x(),
				dynamic_cast<UserInfo&>(Dtag.userInfo()).getFittedMomentum().dpx(),
				Dtag.pType().charge(),
				Dtag.pType().mass());

		fitOne.push_back(Dtag);
	} else { // D*
		Particle Dtag = p.child(0).child(0);
		kv.addTrack(dynamic_cast<UserInfo&>(Dtag.userInfo()).getFittedMomentum().p(),
				dynamic_cast<UserInfo&>(Dtag.userInfo()).getFittedMomentum().x(),
				dynamic_cast<UserInfo&>(Dtag.userInfo()).getFittedMomentum().dpx(),
				Dtag.pType().charge(),
				Dtag.pType().mass());

		fitOne.push_back(Dtag);

		if(ddm<430) {
			Particle Pis = p.child(0).child(1);
			if(abs(Pis.lund())==211) {
				addTrack2fit(kv,Pis);
				fitOne.push_back(Pis);
			} else if(abs(Pis.lund())==111) {
				neutralPi0s.push_back(Pis);
			} else if(abs(Pis.lund())==22) {
				gamma_from_DTagStar = 1;
			}
		} else {
			Particle Pis1 = p.child(0).child(1);
			addTrack2fit(kv,Pis1);
			fitOne.push_back(Pis1);

			Particle Pis2 = p.child(0).child(2);
			addTrack2fit(kv,Pis2);
			fitOne.push_back(Pis2);
		}
	}

	// add Kaon
	addTrack2fit(kv,p.child(1));
	fitOne.push_back(p.child(1));

	// add proton
	if(withProton) {
		addTrack2fit(kv,p.child(2));
		fitOne.push_back(p.child(2));
	}

	// add X
	//int isX = 0;
	//if(p.nChildren()>2)
	//	isX = 1;

	if(withXFrag) {
		Particle xFrag;
		if(withProton)
			xFrag = p.child(3);
		else
			xFrag = p.child(2);

		//std::cout << " Xfrag = " << dynamic_cast<UserInfo&>(xFrag.userInfo()).decayMode() << std::endl;
		std::vector<Particle> fragChildren;
		getFSParticles2(xFrag, fragChildren);

		const int xfrag_nptcle =fragChildren.size();
		for(int  l = 0; l < xfrag_nptcle; l++) {
			Particle xpi = fragChildren[l];
			if(abs(xpi.lund())==211) {
				addTrack2fit(kv,xpi);
				fitOne.push_back(xpi);
			} else if (abs(xpi.lund())==111) {
				neutralPi0s.push_back(xpi);
			}
		}
	}
	//std::cout << " fitOne      = " << fitOne.size() << std::endl;
	//std::cout << " neutralPi0s = " << neutralPi0s.size() << std::endl;

	// Do vertex fit (Dtag K and charged pions from frag. system)
	int err=kv.fit();

	// set usable flag; if fit is OK
	p.usable(!err);

	// update mother
	if(!err) {
		//std::cout << " -> FitOne success!" << std::endl;

		// mass fits of PI0s
		int refitPi0 = 1;
		//std::cout << "-------- Refiting pi0s ---------" << std::endl;
		for(unsigned  i = 0; i < neutralPi0s.size(); i++) {
			Particle g0 = neutralPi0s[i].child(0);
			Particle g1 = neutralPi0s[i].child(1);
			setGammaError(g0, kv.vertex(), kv.errVertex());
			setGammaError(g1, kv.vertex(), kv.errVertex());

			kmassfitter kmf_pi0;
			addTrack2fit(kmf_pi0, g0);
			addTrack2fit(kmf_pi0, g1);
			kmf_pi0.invariantMass(neutralPi0s[i].pType().mass());
			kmf_pi0.vertex(kv.vertex());
			kmf_pi0.errVertex(kv.errVertex());
			err = kmf_pi0.fit();

			p.usable(!err);

			if(!err) {
				kmakemother kmm_pi0;
				makeMother(kmm_pi0,kmf_pi0,neutralPi0s[i],0);
				//std::cout << " - " << i << " Pi0 refit fit success" << std::endl;
			} else {
				refitPi0 = 0;
				p.usable(0);
				//std::cout << " - " << i << " Pi0 refit fit failed" << std::endl;
				break;
			}
		}
		// redo the vertex fit (now include pi0)
		if(refitPi0) {
			//std::cout << " /////// Second fit /////////// " << std::endl;

			kvertexfitter kv2;

			// add all particles to vertex fitter
			// Dtag
			kv2.addTrack(dynamic_cast<UserInfo&>(fitOne[0].userInfo()).getFittedMomentum().p(),
					dynamic_cast<UserInfo&>(fitOne[0].userInfo()).getFittedMomentum().x(),
					dynamic_cast<UserInfo&>(fitOne[0].userInfo()).getFittedMomentum().dpx(),
					fitOne[0].pType().charge(),
					fitOne[0].pType().mass());
			fitAll.push_back(fitOne[0]);
			//std::cout << " -> Adding Dtag = " << fitOne[0].lund() << std::endl;

			for(unsigned i=1; i<fitOne.size(); i++) {
				addTrack2fit(kv2,fitOne[i]);
				fitAll.push_back(fitOne[i]);
				//std::cout << " -> Adding part " << i << " = " << fitOne[i].lund() << std::endl;
			}

			for(unsigned i=0; i<neutralPi0s.size(); i++) {
				addTrack2fit(kv2,neutralPi0s[i]);
				fitAll.push_back(neutralPi0s[i]);
				//std::cout << " -> Adding part " << i << " = " << neutralPi0s[i].lund() << std::endl;
			}

			err = kv2.fit();
			p.usable(!err);
			if(!err) {
				if(kv2.chisq()<0) {
					//std::cout << " *** Second fit failed! " << std::endl;
					p.usable(0);
					return;
				}
				//std::cout << " *** Second fit succeess! " << std::endl;
				dynamic_cast<UserInfo&>(p.userInfo()).setType(p.lund());
				dynamic_cast<UserInfo&>(p.userInfo()).setVChiSq(kv2.chisq());
				dynamic_cast<UserInfo&>(p.userInfo()).setNDF(kv2.dgf());
				p.momentum().decayVertex(kv2.vertex(), kv2.errVertex());

				float chi = kv2.chisq();
				int ndf = kv2.dgf();
				float prob = prob_chi2( chi, ndf );
				//p.usable(chi/ndf < chi2_cut);

				if(p.usable()) {
					//std::cout << " *** Make mother! " << std::endl;
					kmakemother kmm;
					for(unsigned i=0;i<fitAll.size();++i){
						kmm.addTrack(kv2.momentum(i),
								kv2.position(i),
								kv2.error(i),
								fitAll[i].pType().charge());

						kmm.errVertexTrack(kv2.errVertexTrack(i));
						for(unsigned j=i+1;j<fitAll.size();++j){
							kmm.correlation(kv2.correlation(i,j));
						}
					}
					kmm.vertex(kv2.vertex());
					kmm.errVertex(kv2.errVertex());
					unsigned err = kmm.make();

					if(err != 0) {
						//std::cout << " *** Make mother failed! " << std::endl;
						p.usable(0);
						return;
					}

					//std::cout << " *** Make mother success! " << std::endl;
					//std::cout << " Old momentum " << p.momentum().p() << std::endl;
					if(gamma_from_DTagStar==0) {
						p.momentum().momentumPosition(BEAM-kmm.momentum(),
								kmm.position(),
								kmm.error());
					} else {
						setGammaError(p.child(0).child(1), kv2.vertex(), kv2.errVertex());
						p.momentum().momentumPosition(BEAM-kmm.momentum()-p.child(0).child(1).p(),
								kmm.position(),
								kmm.error());
					}
					p.momentum().decayVertex(kv2.vertex(), kv2.errVertex());
					//std::cout << " New momentum " << p.momentum().p() << std::endl;
					dynamic_cast<UserInfo&>(p.userInfo()).setVProb(prob);
				}

			}
		}
	} else {
		//std::cout << "Error doKvFit: " << err << std::endl;
	}


}

void refitDXSystem(Particle &p, HepLorentzVector BEAM) {
	if(! &p.userInfo() ) setUserInfo(p);

	dynamic_cast<UserInfo&>(p.userInfo()).setMassBeforeVertexFit(p.mass());
	std::vector<Particle> neutralPi0s;
	std::vector<Particle> fitOne;
	std::vector<Particle> fitAll;

	kvertexfitter kv;

	// add IP profile to vertex fit
	HepPoint3D bp = IpProfile::position(EVENT_DEPEND_IP);
	HepSymMatrix berr = IpProfile::position_err(EVENT_DEPEND_IP);
	kv.initialVertex(bp);
	kv.vertexProfile(berr);

	// add Dtag (mass vertex fitted momentum)
	int ddm = dynamic_cast<UserInfo&>(p.child(0).userInfo()).decayMode();
	int gamma_from_DTagStar = 0;
	if(ddm<100) { // D0 or D+
		Particle Dtag = p.child(0);
		kv.addTrack(dynamic_cast<UserInfo&>(Dtag.userInfo()).getFittedMomentum().p(),
				dynamic_cast<UserInfo&>(Dtag.userInfo()).getFittedMomentum().x(),
				dynamic_cast<UserInfo&>(Dtag.userInfo()).getFittedMomentum().dpx(),
				Dtag.pType().charge(),
				Dtag.pType().mass());

		fitOne.push_back(Dtag);
	} else { // D*
		Particle Dtag = p.child(0).child(0);
		kv.addTrack(dynamic_cast<UserInfo&>(Dtag.userInfo()).getFittedMomentum().p(),
				dynamic_cast<UserInfo&>(Dtag.userInfo()).getFittedMomentum().x(),
				dynamic_cast<UserInfo&>(Dtag.userInfo()).getFittedMomentum().dpx(),
				Dtag.pType().charge(),
				Dtag.pType().mass());

		fitOne.push_back(Dtag);
		Particle Pis = p.child(0).child(1);
		if(abs(Pis.lund())==211) {
			addTrack2fit(kv,Pis);

			fitOne.push_back(Pis);
		} else if(abs(Pis.lund())==111) {
			neutralPi0s.push_back(Pis);
		} else if(abs(Pis.lund())==22) {
			gamma_from_DTagStar = 1;
		}
	}

	// add X
	int isX = 0;
	if(p.nChildren()>1)
		isX = 1;

	if(isX) {
		Particle xFrag = p.child(1);
		std::vector<Particle> fragChildren;
		getFSParticles2(xFrag, fragChildren);

		const int xfrag_nptcle =fragChildren.size();
		for(int  l = 0; l < xfrag_nptcle; l++) {
			Particle xpi = fragChildren[l];
			if(abs(xpi.lund())==211) {
				addTrack2fit(kv,xpi);
				fitOne.push_back(xpi);
			} else if (abs(xpi.lund())==111) {
				neutralPi0s.push_back(xpi);
			}
		}
	}

	// Do vertex fit (Dtag and charged pions from frag. system)
	int err=kv.fit();

	// set usable flag; if fit is OK
	p.usable(!err);

	// update mother
	if(!err) {
		//std::cout << " -> FitOne success!" << std::endl;

		// mass fits of PI0s
		int refitPi0 = 1;
		//std::cout << "-------- Refiting pi0s ---------" << std::endl;
		for(unsigned  i = 0; i < neutralPi0s.size(); i++) {
			Particle g0 = neutralPi0s[i].child(0);
			Particle g1 = neutralPi0s[i].child(1);
			setGammaError(g0, kv.vertex(), kv.errVertex());
			setGammaError(g1, kv.vertex(), kv.errVertex());

			kmassfitter kmf_pi0;
			addTrack2fit(kmf_pi0, g0);
			addTrack2fit(kmf_pi0, g1);
			kmf_pi0.invariantMass(neutralPi0s[i].pType().mass());
			kmf_pi0.vertex(kv.vertex());
			kmf_pi0.errVertex(kv.errVertex());
			err = kmf_pi0.fit();

			p.usable(!err);

			if(!err) {
				kmakemother kmm_pi0;
				makeMother(kmm_pi0,kmf_pi0,neutralPi0s[i],0);
				//std::cout << " - " << i << " Pi0 refit fit success" << std::endl;
			} else {
				refitPi0 = 0;
				p.usable(0);
				//std::cout << " - " << i << " Pi0 refit fit failed" << std::endl;
				break;
			}
		}
		// redo the vertex fit (now include pi0)
		if(refitPi0) {
			//std::cout << " /////// Second fit /////////// " << std::endl;

			kvertexfitter kv2;

			// add all particles to vertex fitter
			// Dtag
			kv2.addTrack(dynamic_cast<UserInfo&>(fitOne[0].userInfo()).getFittedMomentum().p(),
					dynamic_cast<UserInfo&>(fitOne[0].userInfo()).getFittedMomentum().x(),
					dynamic_cast<UserInfo&>(fitOne[0].userInfo()).getFittedMomentum().dpx(),
					fitOne[0].pType().charge(),
					fitOne[0].pType().mass());
			fitAll.push_back(fitOne[0]);
			//std::cout << " -> Adding Dtag = " << fitOne[0].lund() << std::endl;

			for(unsigned i=1; i<fitOne.size(); i++) {
				addTrack2fit(kv2,fitOne[i]);
				fitAll.push_back(fitOne[i]);
				//std::cout << " -> Adding part " << i << " = " << fitOne[i].lund() << std::endl;
			}

			for(unsigned i=0; i<neutralPi0s.size(); i++) {
				addTrack2fit(kv2,neutralPi0s[i]);
				fitAll.push_back(neutralPi0s[i]);
				//std::cout << " -> Adding part " << i << " = " << neutralPi0s[i].lund() << std::endl;
			}

			err = kv2.fit();
			p.usable(!err);
			if(!err) {
				if(kv2.chisq()<0) {
					//std::cout << " *** Second fit failed! " << std::endl;
					p.usable(0);
					return;
				}
				//std::cout << " *** Second fit succeess! " << std::endl;
				dynamic_cast<UserInfo&>(p.userInfo()).setType(p.lund());
				dynamic_cast<UserInfo&>(p.userInfo()).setVChiSq(kv2.chisq());
				dynamic_cast<UserInfo&>(p.userInfo()).setNDF(kv2.dgf());
				p.momentum().decayVertex(kv2.vertex(), kv2.errVertex());

				float chi = kv2.chisq();
				int ndf = kv2.dgf();
				float prob = prob_chi2( chi, ndf );
				//p.usable(chi/ndf < chi2_cut);

				if(p.usable()) {
					//std::cout << " *** Make mother! " << std::endl;
					kmakemother kmm;
					for(unsigned i=0;i<fitAll.size();++i){
						kmm.addTrack(kv2.momentum(i),
								kv2.position(i),
								kv2.error(i),
								fitAll[i].pType().charge());

						kmm.errVertexTrack(kv2.errVertexTrack(i));
						for(unsigned j=i+1;j<fitAll.size();++j){
							kmm.correlation(kv2.correlation(i,j));
						}
					}
					kmm.vertex(kv2.vertex());
					kmm.errVertex(kv2.errVertex());
					unsigned err = kmm.make();

					if(err != 0) {
						//std::cout << " *** Make mother failed! " << std::endl;
						p.usable(0);
						return;
					}

					//std::cout << " *** Make mother success! " << std::endl;
					//std::cout << " Old momentum " << p.momentum().p() << std::endl;
					if(gamma_from_DTagStar==0) {
						p.momentum().momentumPosition(BEAM-kmm.momentum(),
								kmm.position(),
								kmm.error());
					} else {
						setGammaError(p.child(0).child(1), kv2.vertex(), kv2.errVertex());
						p.momentum().momentumPosition(BEAM-kmm.momentum()-p.child(0).child(1).p(),
								kmm.position(),
								kmm.error());
					}
					p.momentum().decayVertex(kv2.vertex(), kv2.errVertex());
					//std::cout << " New momentum " << p.momentum().p() << std::endl;
					dynamic_cast<UserInfo&>(p.userInfo()).setVProb(prob);
				}

			}
		}
	} else {
		//std::cout << "Error doKvFit: " << err << std::endl;
	}


}

void massFitDXSystem(std::vector<Particle> &plist, HepLorentzVector BEAM) {
	for(unsigned i=0;i<plist.size();++i) massFitDXSystem(plist[i], BEAM);
	rm_unusable( plist );
}

void massFitDXSystem(Particle &p, HepLorentzVector BEAM) {
	if(! &p.userInfo() ) setUserInfo(p);

	dynamic_cast<UserInfo&>(p.userInfo()).setMassBeforeVertexFit(p.mass());
	std::vector<Particle> neutralPi0s;
	std::vector<Particle> fitOne;
	std::vector<Particle> fitAll;

	Particle ipParticle(-BEAM, Ptype(22));

	// add Dtag (mass vertex fitted momentum)
	int ddm = dynamic_cast<UserInfo&>(p.child(0).userInfo()).decayMode();
	int gamma_from_DTagStar = 0;
	if(ddm<100) { // D0 or D+
		Particle Dtag = p.child(0);
		fitOne.push_back(Dtag);
	} else { // D*
		Particle Dtag = p.child(0).child(0);
		fitOne.push_back(Dtag);

		Particle Pis = p.child(0).child(1);
		if(abs(Pis.lund())==211) {
			fitOne.push_back(Pis);
		} else if(abs(Pis.lund())==111) {
			neutralPi0s.push_back(Pis);
		} else if(abs(Pis.lund())==22) {
			gamma_from_DTagStar = 1;
		}
	}

	// add X
	int isX = 0;
	if(p.nChildren()>1)
		isX = 1;

	if(isX) {
		Particle xFrag = p.child(1);
		std::vector<Particle> fragChildren;
		getFSParticles2(xFrag, fragChildren);

		const int xfrag_nptcle =fragChildren.size();
		for(int  l = 0; l < xfrag_nptcle; l++) {
			Particle xpi = fragChildren[l];
			if(abs(xpi.lund())==211) {
				fitOne.push_back(xpi);
			} else if (abs(xpi.lund())==111) {
				neutralPi0s.push_back(xpi);
			}
		}
	}

	int refitPi0 = 1;

	for(unsigned  i = 0; i < neutralPi0s.size(); i++) {
		Particle g0 = neutralPi0s[i].child(0);
		Particle g1 = neutralPi0s[i].child(1);
		setGammaError(g0, p.momentum().decayVertex(), p.momentum().dDecayVertex());
		setGammaError(g1, p.momentum().decayVertex(), p.momentum().dDecayVertex());

		kmassfitter kmf_pi0;
		addTrack2fit(kmf_pi0, g0);
		addTrack2fit(kmf_pi0, g1);
		kmf_pi0.invariantMass(neutralPi0s[i].pType().mass());
		kmf_pi0.vertex(p.momentum().decayVertex());
		kmf_pi0.errVertex(p.momentum().dDecayVertex());
		unsigned err = kmf_pi0.fit();

		p.usable(!err);

		if(!err) {
			kmakemother kmm_pi0;
			makeMother(kmm_pi0,kmf_pi0,neutralPi0s[i],0);
		} else {
			refitPi0 = 0;
			p.usable(0);
			break;
		}
	}

	if(!refitPi0)
		return;

	Ext::kmassfitter kf;
	Particle dummyDS;

	// add "-BEAM"
	addTrack2fit((kmassfitter &)kf, ipParticle);
	dummyDS.relation().append(ipParticle);

	// add DTAG
	kf.addTrack(dynamic_cast<UserInfo&>(fitOne[0].userInfo()).getFittedMomentum().p(),
			dynamic_cast<UserInfo&>(fitOne[0].userInfo()).getFittedMomentum().x(),
			dynamic_cast<UserInfo&>(fitOne[0].userInfo()).getFittedMomentum().dpx(),
			fitOne[0].pType().charge(),
			fitOne[0].pType().mass());
	dummyDS.relation().append(fitOne[0]);

	for(unsigned i=1; i<fitOne.size(); i++) {
		addTrack2fit((kmassfitter &)kf,fitOne[i]);
		dummyDS.relation().append(fitOne[i]);
	}
	for(unsigned i=0; i<neutralPi0s.size(); i++) {
		addTrack2fit((kmassfitter &)kf,neutralPi0s[i]);
		dummyDS.relation().append(neutralPi0s[i]);
	}
	if(gamma_from_DTagStar==1) {
		setGammaError(p.child(0).child(1), p.momentum().decayVertex(), p.momentum().dDecayVertex());
		addTrack2fit((kmassfitter &)kf,p.child(0).child(1));
		dummyDS.relation().append(p.child(0).child(1));
	}
	kf.invariantMass(DSTP_MASS);
    kf.vertex(p.momentum().decayVertex());
    kf.errVertex(p.momentum().dDecayVertex());
    unsigned err = kf.fit();
    if(err == 0){
    	p.usable(!err);

    	unsigned n = kf.tracks();
    	kmakemother kmm;
    	for(unsigned i=0;i<n;++i){
    		kmm.addTrack(kf.momentum(i),
    				kf.position(i),
    				kf.error(i),
    				dummyDS.relation().child(i).pType().charge());
    		if(kf.fitWithVertex())kmm.errVertexTrack(kf.errVertexTrack(i));
    		for(unsigned j=i+1;j<n;++j){
    			kmm.correlation(kf.correlation(i,j));
    		}
    	}
    	kmm.vertex(kf.vertex());
    	if(kf.fitWithVertex()){
    		kmm.errVertex(kf.errVertex());
    	}
    	unsigned err = kmm.make();
    	dummyDS.momentum().momentumPosition(kmm.momentum(),
    			kmm.position(),
    			kmm.error());
    	if(kf.fitWithVertex())
    		dummyDS.momentum().decayVertex(kf.vertex(), kf.errVertex());

    	p.momentum().momentumPosition(-kmm.momentum(),kmm.position(),kmm.error());

    	p.usable(!err);
    } else {
    	p.usable(0);
    }
}


void massFitDKXSystem(std::vector<Particle> &plist, HepLorentzVector BEAM) {
	for(unsigned i=0;i<plist.size();++i) massFitDKXSystem(plist[i], BEAM);
	rm_unusable( plist );
}

void massFitDKXSystem(Particle &p, HepLorentzVector BEAM) {
	if(! &p.userInfo() ) setUserInfo(p);

	int ddmR = dynamic_cast<UserInfo&>(p.userInfo()).decayMode();
	int withXFrag  = 0;
	int withProton = 0;
	if(ddmR==511 || ddmR==512 || ddmR==513 || ddmR==514)
		withXFrag  = 1;
	if(ddmR==503 || ddmR==504 || ddmR==513 || ddmR==514)
		withProton = 1;

	dynamic_cast<UserInfo&>(p.userInfo()).setMassBeforeVertexFit(p.mass());
	std::vector<Particle> neutralPi0s;
	std::vector<Particle> fitOne;
	std::vector<Particle> fitAll;

	Particle ipParticle(-BEAM, Ptype(22));

	//std::cout << " ****** massFitDKXSystem ********* " << std::endl;
	//std::cout << " Information of prior vertex fit: " << std::endl;
	//std::cout << " Ds*.p()            = " << p.p() << std::endl;
	//std::cout << " Ds*.decayVertex()  = " << p.momentum().decayVertex() << std::endl;
	//std::cout << " Ds*.dDecayVertex() = " << p.momentum().dDecayVertex() << std::endl;
	//std::cout << "----------------------------------------------------------" << std::endl;
	//std::cout << " BEAM.p()       = " << BEAM << std::endl;
	//std::cout << " ipParticle.p() = " << ipParticle.p() << std::endl;
	//std::cout << "----------------------------------------------------------" << std::endl;
	//std::cout << " Dtag  = " << dynamic_cast<UserInfo&>(p.child(0).userInfo()).decayMode() << std::endl;
	//std::cout << " Kaon  = " << p.child(1).lund() << std::endl;

	// add Dtag (mass vertex fitted momentum)
	int ddm = dynamic_cast<UserInfo&>(p.child(0).userInfo()).decayMode();
	int gamma_from_DTagStar = 0;
	if(ddm<100) { // D0 or D+
		Particle Dtag = p.child(0);
		fitOne.push_back(Dtag);
	} else { // D*
		Particle Dtag = p.child(0).child(0);
		fitOne.push_back(Dtag);

		if(ddm<430) {
			Particle Pis = p.child(0).child(1);
			if(abs(Pis.lund())==211) {
				fitOne.push_back(Pis);
			} else if(abs(Pis.lund())==111) {
				neutralPi0s.push_back(Pis);
			} else if(abs(Pis.lund())==22) {
				gamma_from_DTagStar = 1;
			}
		} else {
			Particle Pis1 = p.child(0).child(1);
			fitOne.push_back(Pis1);

			Particle Pis2 = p.child(0).child(2);
			fitOne.push_back(Pis2);
		}
	}

	// add Kaon
	fitOne.push_back(p.child(1));

	// add proton
	if(withProton) {
		fitOne.push_back(p.child(2));
	}

	// add X
	if(withXFrag) {
		Particle xFrag;
		if(withProton)
			xFrag = p.child(3);
		else
			xFrag = p.child(2);

		//std::cout << " Xfrag = " << dynamic_cast<UserInfo&>(xFrag.userInfo()).decayMode() << std::endl;
		std::vector<Particle> fragChildren;
		getFSParticles2(xFrag, fragChildren);

		const int xfrag_nptcle =fragChildren.size();
		for(int  l = 0; l < xfrag_nptcle; l++) {
			Particle xpi = fragChildren[l];
			if(abs(xpi.lund())==211) {
				fitOne.push_back(xpi);
			} else if (abs(xpi.lund())==111) {
				neutralPi0s.push_back(xpi);
			}
		}
	}

	//std::cout << " fitOne      = " << fitOne.size() << std::endl;
	//std::cout << " neutralPi0s = " << neutralPi0s.size() << std::endl;

	int refitPi0 = 1;
	//std::cout << "-------- Refiting pi0s ---------" << std::endl;

	for(unsigned  i = 0; i < neutralPi0s.size(); i++) {
		Particle g0 = neutralPi0s[i].child(0);
		Particle g1 = neutralPi0s[i].child(1);
		setGammaError(g0, p.momentum().decayVertex(), p.momentum().dDecayVertex());
		setGammaError(g1, p.momentum().decayVertex(), p.momentum().dDecayVertex());

		kmassfitter kmf_pi0;
		addTrack2fit(kmf_pi0, g0);
		addTrack2fit(kmf_pi0, g1);
		kmf_pi0.invariantMass(neutralPi0s[i].pType().mass());
		kmf_pi0.vertex(p.momentum().decayVertex());
		kmf_pi0.errVertex(p.momentum().dDecayVertex());
		unsigned err = kmf_pi0.fit();

		p.usable(!err);

		if(!err) {
			kmakemother kmm_pi0;
			makeMother(kmm_pi0,kmf_pi0,neutralPi0s[i],0);
			//std::cout << " - " << i << " Pi0 refit fit success" << std::endl;
		} else {
			refitPi0 = 0;
			p.usable(0);
			//std::cout << " - " << i << " Pi0 refit fit failed" << std::endl;
			break;
		}
	}

	if(!refitPi0)
		return;

	//std::cout << "*************** MASS FIT START ************* " << std::endl;
	Ext::kmassfitter kf;
	Particle dummyDS;

	//std::cout << "Point 1" << std::endl;
	// add "-BEAM"
	addTrack2fit((kmassfitter &)kf, ipParticle);
	dummyDS.relation().append(ipParticle);

	//std::cout << "Point 2" << std::endl;
	// add DTAG
	kf.addTrack(dynamic_cast<UserInfo&>(fitOne[0].userInfo()).getFittedMomentum().p(),
			dynamic_cast<UserInfo&>(fitOne[0].userInfo()).getFittedMomentum().x(),
			dynamic_cast<UserInfo&>(fitOne[0].userInfo()).getFittedMomentum().dpx(),
			fitOne[0].pType().charge(),
			fitOne[0].pType().mass());
	// add all other charged tracks
	dummyDS.relation().append(fitOne[0]);

	//std::cout << "Point 3" << std::endl;
	for(unsigned i=1; i<fitOne.size(); i++) {
		addTrack2fit((kmassfitter &)kf,fitOne[i]);
		dummyDS.relation().append(fitOne[i]);
	}
	//std::cout << "Point 4" << std::endl;
	// add all neutral pi0s
	for(unsigned i=0; i<neutralPi0s.size(); i++) {
		addTrack2fit((kmassfitter &)kf,neutralPi0s[i]);
		dummyDS.relation().append(neutralPi0s[i]);
	}
	//std::cout << "Point 5" << std::endl;
    // add gamma from Dtag* decays if exists
	if(gamma_from_DTagStar==1) {
		setGammaError(p.child(0).child(1), p.momentum().decayVertex(), p.momentum().dDecayVertex());
		addTrack2fit((kmassfitter &)kf,p.child(0).child(1));
		dummyDS.relation().append(p.child(0).child(1));
	}
	//std::cout << "Point 6" << std::endl;
	kf.invariantMass(DS_STAR_MASS);
    kf.vertex(p.momentum().decayVertex());
    kf.errVertex(p.momentum().dDecayVertex());
    //std::cout << "Point 7" << std::endl;
    unsigned err = kf.fit();
    //std::cout << "Point 8" << std::endl;
    if(err == 0){
    	p.usable(!err);

    	//std::cout << "-----Make Mother------" << std::endl;
    	//std::cout << "Before m = " << p.mass() << std::endl;
    	//std::cout << "Before p = " << p.p() << std::endl;
    	//int ok = makeMother((kmassfitter &)kf, p);// change Momentum Info.

    	unsigned n = kf.tracks();
    	//std::cout << " kf.tracks() = " << n << std::endl;
    	kmakemother kmm;
    	for(unsigned i=0;i<n;++i){
    		//std::cout << " i = " << i << std::endl;
    		kmm.addTrack(kf.momentum(i),
    				kf.position(i),
    				kf.error(i),
    				dummyDS.relation().child(i).pType().charge());
    		//std::cout << " with vertex " << std::endl;
    		if(kf.fitWithVertex())kmm.errVertexTrack(kf.errVertexTrack(i));
    		//std::cout << " correlations " << std::endl;
    		for(unsigned j=i+1;j<n;++j){
    			//std::cout << " j = " << j << std::endl;
    			kmm.correlation(kf.correlation(i,j));
    		}
    	}
    	kmm.vertex(kf.vertex());
    	if(kf.fitWithVertex()){
    		kmm.errVertex(kf.errVertex());
    	}
    	//std::cout << "Mother make" << std::endl;
    	unsigned err = kmm.make();
    	//std::cout << "Mother make done" << std::endl;
    	//if(err != 0)return 0;
    	dummyDS.momentum().momentumPosition(kmm.momentum(),
    			kmm.position(),
    			kmm.error());
    	if(kf.fitWithVertex())
    		dummyDS.momentum().decayVertex(kf.vertex(), kf.errVertex());

    	//std::cout << "After dummyDS  m = " << dummyDS.mass() << std::endl;
    	//std::cout << "After dummyDS  p = " << dummyDS.p() << std::endl;

    	p.momentum().momentumPosition(-kmm.momentum(),kmm.position(),kmm.error());

    	//std::cout << "After  m = " << p.mass() << std::endl;
    	//std::cout << "After  p = " << p.p() << std::endl;
    	//std::cout << "OK  = " << err << std::endl;
    	p.usable(!err);
    } else {
    	p.usable(0);
    }
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// Vertex FIT
// flag == 1 -> fit to IP, fit == 2 set initial point
// chi2_cut -> cut on chi2/ndf

void doKvFit(Particle &p, float chi2_cut, unsigned flag) {
  if(! &p.userInfo() ) setUserInfo(p);

  dynamic_cast<UserInfo&>(p.userInfo()).setMassBeforeVertexFit(p.mass());

  kvertexfitter kv;

  int nfit=0;
  for(unsigned i=0; i<p.nChildren(); i++) {
//    std::cout<< "B nChild" <<p.child(i).nChildren() << std::endl;
//    if(p.child(i).nChildren()==0){

    addTrack2fit(kv,p.child(i));
    nfit++;
//    }
//    else{
//        for(unsigned j=0; j<p.child(i).nChildren(); j++)
//           { addTrack2fit(kv,p.child(i).child(j));
//            nfit++;
//           }
//        }
  }

  if(flag == 1) {
    // IpProfile::begin_run(); <--- should be called before doKvFit !
    HepPoint3D bp = IpProfile::position(EVENT_DEPEND_IP);
    HepSymMatrix berr = IpProfile::position_err(EVENT_DEPEND_IP);
    kv.initialVertex(bp);
    kv.vertexProfile(berr);
  }
  
  if(flag == 2) {
    kv.initialVertex(p.momentum().decayVertex());
  }

  int err=kv.fit();

  // set usable flag; if fit is OK
  p.usable(!err);
  if(!err) {
//    VectorL PP;
//    for(int i=0; i<nfit; i++)
//      PP += kv.momentum(i);

    dynamic_cast<UserInfo&>(p.userInfo()).setType(p.lund());
    dynamic_cast<UserInfo&>(p.userInfo()).setVChiSq(kv.chisq()/kv.dgf());
    dynamic_cast<UserInfo&>(p.userInfo()).setNDF(kv.dgf());
    p.momentum().decayVertex(kv.vertex(), kv.errVertex());

    float chi = kv.chisq();
    int ndf = kv.dgf();
//    float prob = prob_chi2( chi, ndf );
    p.usable(chi/ndf < chi2_cut);
    if(p.usable()) {
//      dynamic_cast<UserInfo&>(p.userInfo()).setVProb(prob);
     //if(flag == 1) dynamic_cast<UserInfo&>(p.userInfo()).setCorrErrMatrix(kv.errVertexTrack(0));
      int ok=makeMother(kv, p);
      dynamic_cast<UserInfo&>(p.userInfo()).setvtex(kv.vertex());
      p.usable(ok);
    }
  } else {
//    std::cout << "Error doKvFit: " << err << std::endl;
  }
}

void doKvFit(std::vector<Particle> &plist, float chi2_cut, unsigned flag) {
  for(unsigned i=0;i<plist.size();++i) 
    doKvFit(plist[i], chi2_cut, flag);
  rm_unusable( plist );
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// Production Vertex FIT (charm to IP)
// chi2_cut -> cut on probability

void doPvFit(Particle &p, float chi2_cut) {
  if(! &p.userInfo() ) setUserInfo(p);

  kvertexfitter kv;

  addTrack2fit(kv,p);

  // this
  addTube2fit(kv);

  // or this
  //HepPoint3D bp = IpProfile::position(EVENT_DEPEND_IP);
  //HepSymMatrix berr = IpProfile::position_err(EVENT_DEPEND_IP);
  //kv.initialVertex(bp);
  //kv.vertexProfile(berr);

  int err=kv.fit();

  // set usable flag; if fit is OK
  p.usable(!err);
  if(!err) {
    p.momentum().vertex(kv.vertex(), kv.errVertex());

    float chi = kv.chisq();
    int ndf = kv.dgf();
    float prob = prob_chi2( chi, ndf );
    p.usable(chi < chi2_cut);
    if(p.usable()) {
      dynamic_cast<UserInfo&>(p.userInfo()).setPProb(prob);
      dynamic_cast<UserInfo&>(p.userInfo()).setCorrErrMatrix(kv.errVertexTrack(0));
    }
  } else {
    //std::cout << "Error doKvFit: " << err << std::endl;
  }
}

void doPvFit(std::vector<Particle> &plist, float chi2_cut) {
  for(unsigned i=0;i<plist.size();++i) 
    doPvFit(plist[i], chi2_cut);
  
  rm_unusable( plist );
}

int doKmvFit(Particle &p){
	if(! &p.userInfo() ) setUserInfo(p);

	kmassvertexfitter kmv;
	kmv.invariantMass(p.pType().mass());

	for(unsigned i=0; i<p.nChildren(); ++i)
		addTrack2fit(kmv,p.child(i));

	int err=kmv.fit();
	// set usable flag; if fit is OK
	p.usable(!err);
	if(!err){

		dynamic_cast<UserInfo&>(p.userInfo()).setType(p.lund());
		dynamic_cast<UserInfo&>(p.userInfo()).setVChiSq(kmv.chisq());
		dynamic_cast<UserInfo&>(p.userInfo()).setNDF(kmv.dgf());
		p.momentum().decayVertex(kmv.vertex(), kmv.errVertex());

		float chi = kmv.chisq();
		int ndf = kmv.dgf();
		float prob = prob_chi2( chi, ndf );
		dynamic_cast<UserInfo&>(p.userInfo()).setVProb(prob);
		int ok=makeMother(kmv, p);
		p.usable(ok);


	} else {
		//std::cout << "Error doKmvFit: " << err << std::endl;
	}
	return err;
}

void doKmvFit(std::vector<Particle> &plist) {
  for(unsigned i=0;i<plist.size();++i) doKmvFit(plist[i]);
    rm_unusable( plist );
}

/*************** Make Vertex ReFit for D* ****************************/
void refitPislow(Particle &p, int pislow) {
  if(! &p.userInfo() ) setUserInfo(p);

  if(p.child(pislow).lund()!=111) {
    kvertexfitter kbf;
    addTrack2fit(kbf, p.child(pislow));
    kbf.initialVertex(p.child(0).momentum().vertex());
    kbf.knownVertex();

    int err=kbf.fit();
    p.usable(!err);

    if(!err){
      p.momentum(p.child(0).p()+kbf.momentum(0));
      p.momentum().decayVertex(p.child(0).momentum().vertex(), p.child(0).momentum().dVertex());
    }
  } else {
    Particle g0 = p.child(pislow).child(0);
    Particle g1 = p.child(pislow).child(1);
    setGammaError(g0, p.child(0).momentum().vertex(), p.child(0).momentum().dVertex());
    setGammaError(g1, p.child(0).momentum().vertex(), p.child(0).momentum().dVertex());
    
    kmassfitter kmf_pi0;
    addTrack2fit(kmf_pi0, g0);
    addTrack2fit(kmf_pi0, g1);
    kmf_pi0.invariantMass(p.child(pislow).pType().mass());
    kmf_pi0.vertex(p.child(0).momentum().vertex());
    kmf_pi0.errVertex(p.child(0).momentum().dVertex());
    
    int err = kmf_pi0.fit();
    
    p.usable(!err);
    if(!err) {
      Particle pi0_part = p.child(pislow);
      kmakemother kmm_pi0;
      makeMother(kmm_pi0,kmf_pi0,pi0_part,0);
      
      p.momentum(p.child(0).p()+pi0_part.p());
      p.momentum().decayVertex(p.child(0).momentum().vertex(), p.child(0).momentum().dVertex());
    }
  }
}

void refitPislow(Particle &p, HepPoint3D vertex, HepSymMatrix errVertex) {
	if(! &p.userInfo() ) setUserInfo(p);

	if(p.lund()!=111) {
		kvertexfitter kbf;
		addTrack2fit(kbf, p);
		kbf.initialVertex(vertex);
		kbf.knownVertex();

		int err=kbf.fit();
		p.usable(!err);

		if(!err){
			p.momentum(kbf.momentum(0));
			p.momentum().decayVertex(vertex, errVertex);
		}
	} else {
		Particle g0 = p.child(0);
		Particle g1 = p.child(1);
		setGammaError(g0, vertex, errVertex);
		setGammaError(g1, vertex, errVertex);

		kmassfitter kmf_pi0;
		addTrack2fit(kmf_pi0, g0);
		addTrack2fit(kmf_pi0, g1);
		kmf_pi0.invariantMass(p.pType().mass());
		kmf_pi0.vertex(vertex);
		kmf_pi0.errVertex(errVertex);

		int err = kmf_pi0.fit();

		p.usable(!err);
		if(!err) {
			kmakemother kmm_pi0;
			makeMother(kmm_pi0,kmf_pi0,p,0);
		}
	}
}


void refitPislow(std::vector<Particle> &plist, int pislow) {
	for(unsigned i=0;i<plist.size();++i) refitPislow(plist[i], pislow);
	rm_unusable( plist );
}

/************ Set error matrix (dpx) for pi0 ***************************/
int setPi0Error(Particle &p){
	if( !p.child(0) || !p.child(1) ) return 1;
	if( !p.child(0).mdstGamma() ||
			!p.child(1).mdstGamma()    ) return 2;
	for(unsigned i=0; i<2; i++){
		double E     = p.child(i).mdstGamma().ecl().energy();
		double phi   = p.child(i).mdstGamma().ecl().phi();
		double theta = p.child(i).mdstGamma().ecl().theta();

		double E2  = E*E;
		double E4  = E2*E2;
		double ct2 = cos(theta)*cos(theta);
		double st2 = sin(theta)*sin(theta);
		double st4 = st2*st2;

		const HepPoint3D pivot(0.,0.,0.);
		HepMatrix  tmp_a(5,1);
		tmp_a[0][0] = 0.;
		tmp_a[1][0] = phi-M_PI/2;
		tmp_a[2][0] = 1/E/sin(theta);
		tmp_a[3][0] = 0.;
		tmp_a[4][0] = tan(M_PI/2-theta);
		HepVector  a(tmp_a);

		double errE     = p.child(i).mdstGamma().ecl().error(0);
		double errPHI   = p.child(i).mdstGamma().ecl().error(2);
		double errTHETA = p.child(i).mdstGamma().ecl().error(5);

		HepSymMatrix Ea(5,0);
		Ea[0][0] = 1.;
		Ea[1][1] = errPHI;
		Ea[2][2] = errE/E4/st2 + errTHETA*ct2/E2/st4;
		Ea[3][3] = 1.;
		Ea[4][2] = errTHETA*cos(theta)/E/st4;
		Ea[4][4] = errTHETA/st4;

		Helix helix(pivot, a, Ea);

		HepLorentzVector momentum;
		HepPoint3D position;
		HepSymMatrix error(7,0);

		momentum = helix.momentum(0.,0.,position,error);
		p.child(i).momentum().momentumPosition(momentum,position,error);
	}
	return doKmvFit(p);
}

void setPi0Error(std::vector<Particle> &p){
	for(std::vector<Particle>::iterator i=p.begin(); i!=p.end(); ++i)
		setPi0Error(*i);
}

/************ Set error matrix (dpx) for gamma ***************************/
void setMyGammaError(Particle &p){
	if( !p.mdstGamma() )
		return;

	double E     = p.mdstGamma().ecl().energy();
	double phi   = p.mdstGamma().ecl().phi();
	double theta = p.mdstGamma().ecl().theta();

	double E2  = E*E;
	double E4  = E2*E2;
	double ct2 = cos(theta)*cos(theta);
	double st2 = sin(theta)*sin(theta);
	double st4 = st2*st2;

	const HepPoint3D pivot(0.,0.,0.);
	HepMatrix  tmp_a(5,1);
	tmp_a[0][0] = 0.;
	tmp_a[1][0] = phi-M_PI/2;
	tmp_a[2][0] = 1/E/sin(theta);
	tmp_a[3][0] = 0.;
	tmp_a[4][0] = tan(M_PI/2-theta);
	HepVector  a(tmp_a);

	double errE     = p.mdstGamma().ecl().error(0);
	double errPHI   = p.mdstGamma().ecl().error(2);
	double errTHETA = p.mdstGamma().ecl().error(5);

	HepSymMatrix Ea(5,0);
	Ea[0][0] = 1.;
	Ea[1][1] = errPHI;
	Ea[2][2] = errE/E4/st2 + errTHETA*ct2/E2/st4;
	Ea[3][3] = 1.;
	Ea[4][2] = errTHETA*cos(theta)/E/st4;
	Ea[4][4] = errTHETA/st4;

	Helix helix(pivot, a, Ea);

	HepLorentzVector momentum;
	HepPoint3D position;
	HepSymMatrix error(7,0);

	momentum = helix.momentum(0.,0.,position,error);
	p.momentum().momentumPosition(momentum,position,error);
}

void setMyGammaError(std::vector<Particle> &p){
	for(std::vector<Particle>::iterator i=p.begin(); i!=p.end(); ++i)
		setMyGammaError(*i);
}

double KsVtxSigmas(Particle& p){
	double Sigmax = p.momentum().dVertex()[0][0];
	double Sigmay = p.momentum().dVertex()[1][1];
	double Sigmaz = p.momentum().dVertex()[2][2];
	double Sigma  = sqrt(Sigmax*Sigmax+Sigmay*Sigmay+Sigmaz*Sigmaz);
	Vector3 IPvtx(IpProfile::position());
	Vector3 KSvtx(p.momentum().vertex());
	double dist   = (KSvtx - IPvtx).mag();
	return dist/Sigma;
}


double dir2par(Particle &par, int flag)
{
	Vector3 vIp(IpProfile::position());
	Vector3 vPar(par.momentum().vertex());
	Vector3 pPar(par.momentum().p().vect());
	Vector3 dIpPar = vPar - vIp;

	pPar   = Vector3(pPar.x(),   pPar.y(),   0.);
	dIpPar = Vector3(dIpPar.x(), dIpPar.y(), 0.);

	if(flag==1)
		return dIpPar.mag();
	else
		return cos(dIpPar.angle(pPar));
}

double dir2Dsigmas(class Particle &D)
{
	Vector3 vIp(IpProfile::position().x(), IpProfile::position().y(), 0.);
	Vector3 vD (D.momentum().vertex().x(), D.momentum().vertex().y(), 0.);
	Vector3 pD (D.momentum().p().px(),     D.momentum().p().py(),     0.);
	Vector3 vflightD = vD - vIp;

	double flightD = vflightD.mag();
	double sinIpD  = abs(sin(vflightD.angle(pD)));
	double cosIpD  =     cos(vflightD.angle(pD));

	double sigmaIp = sqrt(IpProfile::position_err_b_life_smeared()(1,1) +
			IpProfile::position_err_b_life_smeared()(2,2) );
	double sigmaD  = sqrt(D.momentum().dVertex()(1,1) +
			D.momentum().dVertex()(2,2) );
	double sigma   = sigmaIp + sigmaD;

	double kSigmas(0.);
	if(cosIpD > 0.) kSigmas = flightD * sinIpD / sigma;
	else            kSigmas = flightD          / sigma;

	return kSigmas;
}



#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
