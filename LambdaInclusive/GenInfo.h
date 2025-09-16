#ifndef GENINFO_H
#define GENINFO_H

#include "belle.h"
#include MDST_H
#include EVTCLS_H
#include MDST_OBS_H
#include HEPEVT_H
#include TRK_H
#include <mdst/Muid_mdst.h>
#include "tuple/BelleTupleManager.h"
#include "handAna/AnaConsts.h"
#include "belle.h"
#include "TTree.h"
#include <string>
#include <iomanip>
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>
#include "toolbox/Thrust.h"
#include "toolbox/FuncPtr.h"
#include "handAna/HadronPair.h"
#include "handAna/HadronQuadruple.h"
#include "handAna/Savable.h"
#include "handAna/ParticleInfoMass.h"
#include "handAna/DebugHistos.h"


#include "fastjet/ClusterSequence.hh"
#include <iostream>

using namespace fastjet;
using namespace std;

//taken from Jimmies code...
#define ID_GAMMA 22
#define ID_PI 211
#define ID_PI0 111
#define ID_K 321
#if defined(BELLE_NAMESPACE)
namespace Belle {

#endif
  /*
    class that encapsulates generator info and methods to fill it
  */
  Hep3Vector& gretSelf(Hep3Vector& vec);

  class GenInfo
  {
  public:
    GenInfo():cPiPlus("PI+"), cPiNeg("PI-"),cPiZero("PI0"),cKPlus("K+"),cKNeg("K-")
      {
	//mGZBins.push_back(0.2);
	mGZBins.push_back(0.3);
	mGZBins.push_back(0.4);
	mGZBins.push_back(0.55);
	mGZBins.push_back(0.75);
	mGZBins.push_back(1.1);
      }
      int getBin(vector<float>& b1, float value)
      {
	int coo1=-1;
	for(int i=0;i<b1.size();i++)
	  {
	    if(value<=b1[i])
	      {
		coo1=i;
		break;
	      }
	  }
	return coo1;
      }

      vector<float> mGZBins;
      Hep3Vector cmThrust;
      Hep3Vector jet1;
      Hep3Vector jet2;
      Hep3Vector dijet[2];
      float fluffiness1;
      float fluffiness2;
      float jetE1;
      float jetE2;
      float jetPhi1;
      float jetPhi2;
      float jetTheta1;
      float jetTheta2;


      float cmThrustMag;
      Hep3Vector labThrust;
      float thrustPR;
      float vpEnergy;
      float vpPx;
      float vpPy;
      float vpPz;
      int numQuarks;
      float quarkAngle;
      float thetaEThrust;


      //the adresses from which the tree should read its data
      static Savable tData;
      int num_f_data;
      int num_fA_data;
      int num_iA_data;
      int num_i_data;

      void setDebugHistos(DebugHistos* m_histos)
      {
	this->m_histos=m_histos;
      }

      void fillInf()
      {
	computeGenThrust();
	thrustPR=-log(tan(cmThrust.theta()/2));
      }

      void doAll()
      {
	//hopefully after thust computation
	//  if(kinematics::thrustMag<cuts::minThrust || abs(kinematics::thrustDirCM.z())/kinematics::thrustDirCM.mag()>cuts::maxThrustZ|| visEnergyOnFile<cuts::minVisEnergy || iChTrks < cuts::minNTracks)
	//if(cmThrustMag<cuts::minThrust|| abs(cmThrust.z()/cmThrust.mag()> cuts::maxThrustZ))
	//return;
	vector<Particle*> v_allParticles;
	vector<Particle*> v_firstHemiPos;
	vector<Particle*> v_firstHemiNeg;
	vector<Particle*> v_secondHemiPos;
	vector<Particle*> v_secondHemiNeg;
	vector<HadronPair*> v_hadronPairsFirstH_PN;
	vector<HadronPair*> v_hadronPairsSecondH_PN;
	vector<HadronQuadruple*> v_hadronQuadruplesPN;
	getHadronLists(v_allParticles);
	//cout <<" w/o acceptance " << endl;
	//cout <<"we have " << v_allParticles.size() <<" particles " <<endl;
	setParticleProperties(v_allParticles, cmThrust,v_firstHemiPos, v_firstHemiNeg, v_secondHemiPos, v_secondHemiNeg);
	findHadronPairs(v_firstHemiPos, v_firstHemiNeg,v_secondHemiPos,v_secondHemiNeg, v_hadronPairsFirstH_PN, v_hadronPairsSecondH_PN);
	//      cout <<"we have " << v_firstHemiPos.size() << " pos "<<v_firstHemiNeg.size() <<" neg in first hemi " << v_secondHemiPos.size() << " pos in sec " << v_secondHemiNeg.size()<<" neg in second"<<endl;
	//      cout <<" we have " << v_hadronPairsFirstH_PN.size() <<" pairs in first hemi " << v_hadronPairsSecondH_PN.size() <<" in second hemi " <<endl;
	combineHadronPairs(v_hadronPairsFirstH_PN, v_hadronPairsSecondH_PN, v_hadronQuadruplesPN);
	//      cout << " we have " << v_hadronQuadruplesPN.size() <<" quads " <<endl<<endl;
	//                  cout <<" pairs: " << v_hadronPairsFirstH_PN.size() <<" sec: " << v_hadronPairsSecondH_PN.size()<<endl;
	//                  cout <<" num quads: " << v_hadronQuadruplesPN.size() <<endl;

	fillWQuadrupleData(v_hadronQuadruplesPN);
	//      cout <<"push back event data.. " <<endl; 
	int numF=num_fA_data;
	int numI=num_iA_data;

	//insert here


	tData.dataF.push_back(cmThrust.theta());
	tData.dataF.push_back(labThrust.theta());
	tData.dataF.push_back(cmThrustMag);
	tData.dataF.push_back(thetaEThrust);
	//      cout <<"save data.. " <<endl;
	if(v_hadronQuadruplesPN.size()>0)
	  saveData(tData.dataF,tData.dataI,2*(numI+numF));
	//      cout <<"now clean up. " <<endl;
	//get ready for next call where dataF has to be clean (this is static)
	tData.dataF.clear();
	tData.dataI.clear();
	//and the event data has to be filled afterwards:
	cleanUp(v_allParticles,v_firstHemiPos, v_firstHemiNeg,v_secondHemiPos, v_secondHemiNeg, v_hadronPairsFirstH_PN, v_hadronPairsSecondH_PN, v_hadronQuadruplesPN);
      }
      //reference particle types
      Ptype cPiPlus;
      Ptype cPiNeg;
      Ptype cPiZero;
      Ptype cKPlus;
      Ptype cKNeg;

      void initializeTree()
      {
	if(tData.initialized)
	  return;
	num_f_data=0;
	num_fA_data=0;
	num_iA_data=0;
	num_i_data=0;
	//the first time the class is initialized, construct the tree
	tData.pDataTree=new TTree("GenTree","Generated Tree");
	addArrayF("z1_mcWoA");
	addArrayF("z2_mcWoA");
	addArrayF("z1Ratio_mcWoA");
	addArrayF("z2Ratio_mcWoA");
	addArrayF("mass1_mcWoA");
	addArrayF("mass2_mcWoA");
	addArrayF("phiRSum_mcWoA");
	addArrayF("phiR1_mcWoA");
	addArrayF("twoPhiZero_mcWoA");
	addArrayF("phiZero1_mcWoA");

	addArrayF("theta11_mcWoA");
	addArrayF("theta12_mcWoA");
	addArrayF("theta21_mcWoA");
	addArrayF("theta22_mcWoA");
	addArrayF("decayTheta1_mcWoA");
	addArrayF("decayTheta2_mcWoA");
	addArrayF("thrustProj11_mcWoA");
	addArrayF("thrustProj12_mcWoA");
	addArrayF("thrustProj21_mcWoA");
	addArrayF("thrustProj22_mcWoA");
	//hadron quad level
	//qT
	addArrayF("qT_mcWoA");

	//hadron pair level seems to be fine

	addArrayI("motherGenId_11_mcWoA");
	addArrayI("motherGenId_12_mcWoA");
	addArrayI("motherGenId_21_mcWoA");
	addArrayI("motherGenId_22_mcWoA");

	addArrayI("chargeType_mcWoA");
	addArrayI("particleType_mcWoA");
	addArrayI("chargeType1_mcWoA");
	addArrayI("particleType1_mcWoA");
	addArrayI("chargeType2_mcWoA");
	addArrayI("particleType2_mcWoA");
    
	//missing on event level:
	//thrustTheta
	//thrustThetaLab
	//thetaEThrust

	//fields after arrays to make the filling work....
	addFieldF("thrustTheta_mcWoA");
	addFieldF("thrustThetaLab_mcWoA");
	addFieldF("thrustMag_mcWoA");
	addFieldF("thetaEThrust_mcWoA");


  
	tData.initialized=true;
      };
  private:
      void addFieldF(char* fieldname)
      {
	num_f_data++;
	//construct the memory location from which the tree should read the new data field
	float* memLoc=new float;
	tData.treeData.push_back(memLoc);
	tData.pDataTree->Branch(fieldname, memLoc, (fieldname+string("/F")).c_str());
	tData.fieldNamesF.push_back(fieldname);
      };


      void addFieldI(char* fieldname)
      {
	num_i_data++;
	int* memLoc=new int;
	tData.treeData.push_back(memLoc);
	tData.pDataTree->Branch(fieldname,memLoc,(fieldname+string("/I")).c_str());
	tData.fieldNamesI.push_back(fieldname);
      };

      void addArrayI(char* fieldname)
      {
	num_iA_data++;
	//standard lenth, shouldn't be more than that
	int* memLoc=new int[1200];
	int* memLocCounter=new int;
	tData.treeData.push_back(memLocCounter);
	tData.treeData.push_back(memLoc);
	string counterName=string(fieldname)+string("Counter");
	tData.pDataTree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
	tData.pDataTree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/I")).c_str());
      };

      void addArrayF(char* fieldname)
      {
	num_fA_data++;
	float* memLoc=new float[1200];
	int* memLocCounter=new int;
	tData.treeData.push_back(memLocCounter);
	tData.treeData.push_back(memLoc);
	string counterName=string(fieldname)+string("Counter");
	tData.pDataTree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
	tData.pDataTree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/F")).c_str());
      };

    

      void fillWQuadrupleData(vector<HadronQuadruple*>& v_hq, bool mcPart=false)
      {
	int numI=-1;
	int numF=-1;
	int counter=0;//offset from one quadruple to the next
	for(vector<HadronQuadruple*>::iterator it=v_hq.begin();it!=v_hq.end();it++)
	  {
	    HadronQuadruple* p_hq=*it;
	    if(p_hq==0)
	      {
		for(int i=0;i<num_fA_data;i++)
		  {
		    tData.dataF.push_back(-1);
		  }
		for(int i=0;i<num_iA_data;i++)
		  {
		    tData.dataI.push_back(-1);
		  }
	      }
	    else
	      {
		tData.dataF.push_back(p_hq->firstHPair->z);//z1
		tData.dataF.push_back(p_hq->secondHPair->z); //z2
		tData.dataF.push_back(p_hq->firstHPair->z/dynamic_cast<ParticleInfo&>(p_hq->firstHPair->firstHadron->userInfo()).z);//z1Ratio
		tData.dataF.push_back(p_hq->secondHPair->z/dynamic_cast<ParticleInfo&>(p_hq->secondHPair->firstHadron->userInfo()).z); //z2Ratio

		tData.dataF.push_back(p_hq->firstHPair->mass);      //mass1
		tData.dataF.push_back(p_hq->secondHPair->mass);      //mass2
		tData.dataF.push_back(p_hq->phiSum);            //phiR_Sum
		tData.dataF.push_back(p_hq->firstHPair->phiR);
		tData.dataF.push_back(p_hq->phiZero1);
		tData.dataF.push_back(p_hq->phiZeroR);

		float theta1=dynamic_cast<ParticleInfo&>(p_hq->firstHPair->firstHadron->userInfo()).theta;
		float theta2=dynamic_cast<ParticleInfo&>(p_hq->firstHPair->secondHadron->userInfo()).theta;
		float theta3=dynamic_cast<ParticleInfo&>(p_hq->secondHPair->firstHadron->userInfo()).theta;
		float theta4=dynamic_cast<ParticleInfo&>(p_hq->secondHPair->secondHadron->userInfo()).theta;
		if(theta1==0 || theta2==0 || theta3 == 0 || theta4==0)
		  {
		    tData.dataF.push_back(0);
		    tData.dataF.push_back(0);
		    tData.dataF.push_back(0);
		    tData.dataF.push_back(0);
		  }
		else
		  {
		    tData.dataF.push_back(theta1);
		    tData.dataF.push_back(theta2);
		    tData.dataF.push_back(theta3);
		    tData.dataF.push_back(theta4);
		  }
		//mass makes only sense for reconstructed particles, in mc they should always have the Pi0 amss
		//	    tData.dataF.push_back(dynamic_cast<ParticleInfoMass&>(p_hq->firstHPair->firstHadron->userInfo()).mass);
		//	    tData.dataF.push_back(dynamic_cast<ParticleInfoMass&>(p_hq->firstHPair->secondHadron->userInfo()).mass);
		//	    tData.dataF.push_back(dynamic_cast<ParticleInfoMass&>(p_hq->secondHPair->firstHadron->userInfo()).mass);
		//	    tData.dataF.push_back(dynamic_cast<ParticleInfoMass&>(p_hq->secondHPair->secondHadron->userInfo()).mass);

	  
		//	cout <<"1" <<endl;
		tData.dataF.push_back(p_hq->firstHPair->decayTheta);
		tData.dataF.push_back(p_hq->secondHPair->decayTheta);
		tData.dataF.push_back(dynamic_cast<ParticleInfo&>(p_hq->firstHPair->firstHadron->userInfo()).thrustProj);
		tData.dataF.push_back(dynamic_cast<ParticleInfo&>(p_hq->firstHPair->secondHadron->userInfo()).thrustProj);
		tData.dataF.push_back(dynamic_cast<ParticleInfo&>(p_hq->secondHPair->firstHadron->userInfo()).thrustProj);
		tData.dataF.push_back(dynamic_cast<ParticleInfo&>(p_hq->secondHPair->secondHadron->userInfo()).thrustProj);
		tData.dataF.push_back(p_hq->qT);
		tData.dataI.push_back(dynamic_cast<ParticleInfo&>(p_hq->firstHPair->firstHadron->userInfo()).motherGenId);
		tData.dataI.push_back(dynamic_cast<ParticleInfo&>(p_hq->firstHPair->secondHadron->userInfo()).motherGenId);
		tData.dataI.push_back(dynamic_cast<ParticleInfo&>(p_hq->secondHPair->firstHadron->userInfo()).motherGenId);
		tData.dataI.push_back(dynamic_cast<ParticleInfo&>(p_hq->secondHPair->secondHadron->userInfo()).motherGenId);
		tData.dataI.push_back(p_hq->hadCharge);      
		tData.dataI.push_back(p_hq->hadPType);//particle type

		/*	  if(!mcPart)
		  {
		  }*/
		//	cout <<"2" <<endl;
		int firstPType=p_hq->firstHPair->hadPType;
		int secondPType=p_hq->secondHPair->hadPType;
		int firstCharge=p_hq->firstHPair->hadCharge;
		int secondCharge=p_hq->secondHPair->hadCharge;
		/* see in Tree Saver, don't know why this should make sense
		   if(firstPType<secondPType)
		   {
		   int tmp=firstPType;
		   firstPType=secondPType;
		   secondPType=tmp;
		   }
		   if(firstCharge<secondCharge)
		   {
		   int tmp=firstCharge;
		   firstCharge=secondCharge;
		   secondCharge=tmp;
		   }
		*/
		//	cout <<"3" <<endl;



		tData.dataI.push_back(firstCharge);      
		tData.dataI.push_back(firstPType);
		tData.dataI.push_back(secondCharge);      
		tData.dataI.push_back(secondPType);
	      }
	    //	cout <<"4"<<endl;      
	    numF=tData.dataF.size();  //always the same, 
	    numI=tData.dataI.size(); 
	    //how many array values do we have? At this point we haven't added event data yet...
	    for(int j=0;j<tData.dataF.size();j++)
	      {
		//	    cout <<"looking at field nr " << j <<endl;
		if(counter >=1200)
		  {
		    cout <<"index q to large" <<endl;
		    continue;
		  }
		//tree data has all fields
		if(2*j+1 > (tData.treeData.size()-num_f_data))
		  {
		    cout <<"index td to large" << tData.treeData.size() << ", " << 2*j+1<<" num scalar fields: "<< num_f_data <<endl;
		    continue;
		  }
		((float*)tData.treeData[2*j+1])[counter]=tData.dataF[j];
	      }
	    //tData.dataI has size 4
	    //	cout <<" off to integeres.. " <<endl;
	    for(int j=0;j<tData.dataI.size();j++)
	      {
		if(2*(j+numF)+1 > tData.treeData.size())
		  {
		    cout <<"index2 td to large" << tData.treeData.size() << ", " << 2*(j+numF)+1 <<endl;
		    continue;
		  }

		((int*)tData.treeData[2*(j+numF)+1])[counter]=tData.dataI[j];
	      }
	    //has to be cleared to save next quad!
	    tData.dataF.clear();
	    tData.dataI.clear();
	    counter++;
	  }
	//	cout <<"5" <<endl;
	counter=v_hq.size();
	//save counter info
	for(int j=0;j<numF;j++)
	  {
	    *(int*)tData.treeData[2*j]=counter;
	  }
	for(int j=0;j<numI;j++)
	  {
	    *(int*)tData.treeData[2*(j+numF)]=counter;
	  }
    
	if(numF < 0) //no events, all quad vector empty
	  return;


      }

      void saveData(vector<float>& dataF, vector<int>& dataI, int offset)
      {
	int dataSize=dataF.size()+dataI.size();
	//      cout <<"data size: "<< dataSize <<endl;
	if(dataSize !=tData.treeData.size()-offset)
	  {
	    cout << "data size does not match number of branches " <<dataSize <<"+ " << offset << " = " <<tData.treeData.size() <<endl <<flush;
	    exit(0);
	  }
	//      cout <<"puting " << dataF.size() <<" floats " << endl;
	for(int i=0;i<dataF.size();i++)
	  {
	    *(float*)tData.treeData[i+offset]=dataF[i];
	  }
	//      cout <<" putting " << dataI.size() << " integers " << endl;
	for(int i=0;i<dataI.size();i++)
	  {
	    *(int*)tData.treeData[i+dataF.size()+offset]=dataI[i];
	  }
	//      cout <<"fill..." <<endl;
	tData.pDataTree->Fill();
	//      cout <<"done filling " <<endl;
      };
      //unfortunately just a copy of the handAna functions
      void getHadronLists(vector<Particle*>& v_allParticles)
      {
	int omegaCounter=0;
	//      cout <<"getHadronlist" <<endl;
	int motherGenId=0;
	Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();
	//      cout <<"iterating" <<endl;
	for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
	  {
	    //	  cout << "in loop " <<endl;
	    //here I can probably hack bacchettas stuff in, just require that the pion parents are the ones that he wants...
	    int geantID=abs(gen_it->idhep());//plus is ok, since it is the abs value
	    //normal w/o k0:

	    if(gen_it->mother())
	      {
		motherGenId=gen_it->mother().idhep();
		/*	      	      cout <<"found mother with id: " << gen_it->mother().idhep() <<endl;
		  if(gen_it->mother().mother())
		  cout <<"mothers mother: " << gen_it->mother().mother().idhep()<<endl;
		  else
		  cout <<"no grandmother " <<endl;*/
		//311: K0
		//310: K_S0
		//313 K*
		//130: K_L0
		//323 K*+
		//321: K+
		// 113 rho
		if(abs(gen_it->mother().idhep()) == 311 || abs(gen_it->mother().idhep())==310|| abs(gen_it->mother().idhep())==313 || abs(gen_it->mother().idhep())==130)//k0
		  {
		    //		continue;
		  }
		if(abs(gen_it->mother().idhep()) == 323 || abs(gen_it->mother().idhep())==321)//k0
		  {
		    //		continue;
		  }
		//rho: 
		//	      if(abs(gen_it->mother().idhep())!=223) //omega
		if(abs(gen_it->mother().idhep())!=113) //rho
		  {
		    //	  omegaCounter++;
		    //		  cout <<"found omega " << gen_it->mother()<<" E: " << gen_it->mother().E() << "  <<omegaCounter: " << omegaCounter <<endl;
		    //  continue;
		  }
	      }
	    else
	      {
		//	      continue;
	      }


	    if(!(geantID==lc_pi0 || geantID==lc_piPlus || geantID==lc_kPlus))
	      continue;
	    //	   cout <<"valid " <<endl;
	    Particle* np=new Particle(*gen_it);
	    HepLorentzVector boostedVec(gen_it->PX(),gen_it->PY(),gen_it->PZ(),gen_it->E());
	    float m_theta=boostedVec.theta();
	    //	  cout <<"got tehta: " << m_theta << " is there a diff with vec3? " << boostedVec.vect().theta()<<endl;
	    boostedVec.boost(kinematics::CMBoost);
	    np->userInfo(*(new ParticleInfoMass())); //gets deleted in destructor of Particle
	    ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(np->userInfo());
	    float m_z=2*boostedVec.t()/kinematics::Q;
	    pinf.motherGenId=motherGenId;
	    pinf.z=m_z;
	    //theta should be the one in the lab system as before...
	    pinf.theta=m_theta;
	    if(jet1.dot(boostedVec.vect())>0)
	      pinf.thrustProj=jet1.dot(boostedVec.vect())/(jet1.mag()*boostedVec.vect().mag());
	    else
	      pinf.thrustProj=jet2.dot(boostedVec.vect())/(jet2.mag()*boostedVec.vect().mag());

	    np->momentum().momentum(boostedVec);
	    v_allParticles.push_back(np);
	  }
      }

      void setParticleProperties(vector<Particle*>& v_allParticles, Hep3Vector& mThrustDir,vector<Particle*>& v_firstHemiPos, vector<Particle*>& v_firstHemiNeg, vector<Particle*>& v_secondHemiPos, vector<Particle*>& v_secondHemiNeg)
      {
	//      cout <<" there are " << v_allParticles.size() <<" MC particles " <<endl;


	for(vector<Particle*>::const_iterator it=v_allParticles.begin();it!=v_allParticles.end();it++)
	  {
	    Hep3Vector axis=jet1;
	    if(jet1.dot((*it)->p().vect())<0)
	      axis=jet2;
	    ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	    float phi=AuxFunc::getPhi(axis,**it);

	    //      float phi=AuxFunc::getPhi(mThrustDir,**it);
	    //      cout <<"phi is: " << phi<< endl;
	    //apparently in getHadronlist, the cm theta is saved here:
	    float cmTheta=pinf.theta;
	    float theta=AuxFunc::getTheta(axis,**it);
	    pinf.phi=phi;
	    //	  pinf.theta=theta;
	    //roughly cuts as for the analys, for wulfrim <<-- Disabled now...
	    if(pinf.z > 0.1)//was 0.1
	      {
		if(!(cos(cmTheta)<cuts::minCosTheta||cos(cmTheta)>cuts::maxCosTheta))
		  {
		    if(((*it)->pType()==cPiPlus)||((*it)->pType()==cPiNeg))
		      {
			//		      m_histos->hEFlowMC->Fill(theta,pinf.z);
		      }
		  }
	      }
	    else
	      {
		//		continue;
	      }
	    pinf.thrustProj=axis.dot((*it)->p().vect())/(axis.mag()*(*it)->p().vect().mag());	  

	    if(fabs(pinf.thrustProj)<cuts::minThrustProj)
	      {
		//no cuts...? doch...
		//		continue;
	      }
	    if(pinf.thrustProj>0)//one hemisphere
	      {
		if((*it)->pType().charge()>0)
		  v_firstHemiPos.push_back(*it);
		else
		  v_firstHemiNeg.push_back(*it);
	      }
	    else//particle is in the other hemisphere
	      {
		if((*it)->pType().charge()>0)
		  v_secondHemiPos.push_back(*it);
		else
		  v_secondHemiNeg.push_back(*it);
	      }
	  }
      }

      void findHadronPairs(vector<Particle*>& v_firstHemiPos, vector<Particle*>& v_firstHemiNeg,vector<Particle*>& v_secondHemiPos,vector<Particle*>& v_secondHemiNeg, vector<HadronPair*>& v_hadronPairsFirstH_PN, vector<HadronPair*>& v_hadronPairsSecondH_PN)
      {
	//collect all +/- pairs in first Hemisphere
	//      cout <<"find hadron pairs with " << v_firstHemiPos.size() << " in pos and " << v_firstHemiNeg.size() <<endl;
	for(vector<Particle*>::const_iterator it=v_firstHemiPos.begin();it!=v_firstHemiPos.end();it++)
	  {
	    //don't continue if it is not pion or kaon
	    //  cout << (*it)->pType() <<" pion: " << cPiPlus <<" kaon: " << cKPlus <<endl;
	    if(!((*it)->pType()==cPiPlus || (*it)->pType()==cKPlus))
	      {
		continue;
	      }
	    ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	    for(vector<Particle*>::const_iterator it2=v_firstHemiNeg.begin();it2!=v_firstHemiNeg.end();it2++)
	      {
		ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());
		if(pinf.z+pinf2.z < cuts::min2H_Z)
		  {
		    //continue;
		  }
		if(pinf.z+pinf2.z < 0.001)
		  {
		    continue;
		  }
		if(!((*it2)->pType()==cPiNeg || (*it2)->pType()==cKNeg))
		  {
		    continue;
		  }
		//now unknowns...
		HadronPair* hp=new HadronPair();
		hp->firstHadron=*it;
		hp->secondHadron=*it2;


		hp->hadCharge=AnaDef::PN;
		hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType());
		//first hemi
		hp->computeR(jet1);
		hp->computeThrustTheta(jet1);
		if(isnan(hp->phiR))
		  {
		    cout <<" hp r is not a number!!!, z1: " << pinf.z <<" z2: " << pinf2.z << " cm Thrust x: " << cmThrust.x() <<" y: " << cmThrust.y() << " z: " << cmThrust.z();
		    cout <<" first had x: " << hp->firstHadron->p().vect().x() <<" y: " <<hp->firstHadron->p().vect().y()<< " z: "<< hp->firstHadron->p().vect().z()  <<endl;
		    cout <<"second had: " <<  hp->secondHadron->p().vect().x() <<" y: " <<hp->secondHadron->p().vect().y()<< " z: "<< hp->secondHadron->p().vect().z()  <<endl;
		  }
		//	  hp->computeR(mThrustDir);
		if((*it)->pType()==cPiPlus)
		  {
		    if((*it2)->pType()==cPiNeg)
		      {
			//		      if(hp->mass<0.51 && hp->mass > 0.49)
			//			cout <<"have " << hp->mass<<" mass with ph: " << hp->P_h<<" decTheta: " << hp->decayTheta <<" z: " << hp->z <<" rvect: " << hp->Rvect<<endl;
			//the 4pi thing for bacchetta//gen_it->mother().idhep()
			//		      cout <<"found pair " <<endl;
			//		      cout <<"mass i: " << hp->mass <<endl;
			//		      if(hp->mass< 0.495  || hp->mass > 0.499) //cut out the k0
			{
			  //  cout <<"putting in histo with m: " << hp->mass <<endl;
			  //to record the number of pairs at a specific z
			  int iZBin=getBin(mGZBins,hp->z);
			  m_histos->hZBinValues->Fill(iZBin);
			  if(iZBin<4 && iZBin > -1)
			    {
			      m_histos->hHPairMassesMC[iZBin]->Fill(hp->mass);
			    }
			  m_histos->hHPairMassMC->Fill(hp->mass);

			}
		      }
		  }
		v_hadronPairsFirstH_PN.push_back(hp);
	      }
	  }
	for(vector<Particle*>::const_iterator it=v_secondHemiPos.begin();it!=v_secondHemiPos.end();it++)
	  {
	    if(!((*it)->pType()==cPiPlus || (*it)->pType()==cKPlus))
	      {
		continue;
	      }
	    ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	    for(vector<Particle*>::const_iterator it2=v_secondHemiNeg.begin();it2!=v_secondHemiNeg.end();it2++)
	      {
		ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());
		if(pinf.z+pinf2.z < cuts::min2H_Z)
		  {
		    // continue;
		  }
		if(pinf.z+pinf2.z < 0.001)
		  {
		    continue;
		  }
		if(!((*it2)->pType()==cPiNeg || (*it2)->pType()==cKNeg))
		  {
		    continue;
		  }
		HadronPair* hp=new HadronPair();
		hp->firstHadron=*it;
		hp->secondHadron=*it2;

		hp->hadCharge=AnaDef::PN;
		hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType());
		hp->computeR(jet2);
		//	  hp->computeR(mThrustDir);
		hp->computeThrustTheta(jet2);
		if(isnan(hp->phiR))
		  {
		    cout <<" hp r is not a number!!!, z1: " << pinf.z <<" z2: " << pinf2.z << " cm Thrust x: " << cmThrust.x() <<" y: " << cmThrust.y() << " z: " << cmThrust.z();
		    cout <<" first had x: " << hp->firstHadron->p().vect().x() <<" y: " <<hp->firstHadron->p().vect().y()<< " z: "<< hp->firstHadron->p().vect().z()  <<endl;
		    cout <<"second had: " <<  hp->secondHadron->p().vect().x() <<" y: " <<hp->secondHadron->p().vect().y()<< " z: "<< hp->secondHadron->p().vect().z()  <<endl;
		  }

		if((*it)->pType()==cPiPlus)
		  {
		    if((*it2)->pType()==cPiNeg)
		      {
			//		      if(hp->mass< 0.495  || hp->mass > 0.499) //cut out the k0
			m_histos->hHPairMassMC->Fill(hp->mass);
		      }
		  }

		v_hadronPairsSecondH_PN.push_back(hp);
	      }
	  }
      }

      void combineHadronPairs(vector<HadronPair*>& v_hadronPairsFirstH_PN, vector<HadronPair*>& v_hadronPairsSecondH_PN, vector<HadronQuadruple*>& v_hadronQuadruplesPN)
      {
	//      cout <<v_hadronPairsFirstH_PN.size() <<" had pairs in firt " << v_hadronPairsSecondH_PN.size() <<" in second hemi " << endl;
	//  int  evtNr=Belle_event_Manager::get_manager().begin()->EvtNo();
	for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_PN.begin();it!=v_hadronPairsFirstH_PN.end();it++)
	  {
	    for(vector<HadronPair*>::iterator it2=v_hadronPairsSecondH_PN.begin();it2!=v_hadronPairsSecondH_PN.end();it2++)
	      {
		AnaDef::TwoHadCharge hadCharge;
		AnaDef::TwoHadPType hadPType;
		if(((*it)->hadPType==AnaDef::KPi && (*it2)->hadPType==AnaDef::PiK)  ||((*it2)->hadPType==AnaDef::KPi && (*it)->hadPType==AnaDef::PiK))
		  {
		    hadCharge=AnaDef::PNNP;
		    hadPType=AnaDef::PiK;
		  }
		else
		  {
		    if((*it)->hadPType!=(*it2)->hadPType)
		      {
			continue; //not really efficient, but easy fix.... should have separate vectors for all particle types
		      }
		    hadCharge=AnaDef::PN;
		    hadPType=(*it)->hadPType;
		  } // now we take all types

		HadronQuadruple* hq=new HadronQuadruple();
		hq->hadCharge=hadCharge;
		hq->hadPType=hadPType;
		hq->firstHPair=*it;
		hq->secondHPair=*it2; 
		hq->phiSum=(*it)->phiR+(*it2)->phiR;
		if(isnan(hq->phiSum))
		  cout <<"sum not a number!!: " << (*it)->phiR <<" r2: " << (*it2)->phiR <<endl;
		hq->compQT();
		if(hq->phiSum>pi)
		  hq->phiSum-=2*pi;
		if(hq->phiSum<-pi)
		  hq->phiSum+=2*pi;
		if(hq->qT>cuts::maxQt)
		  {
		    delete hq;
		    continue;
		  }
		v_hadronQuadruplesPN.push_back(hq);
	      }
	  }
  


      }

      void cleanUp(vector<Particle*>& v_allParticles,vector<Particle*>& v_firstHemiPos, vector<Particle*>& v_firstHemiNeg,vector<Particle*>& v_secondHemiPos, vector<Particle*>& v_secondHemiNeg, vector<HadronPair*>& v_hadronPairsFirstH_PN, vector<HadronPair*>& v_hadronPairsSecondH_PN, vector<HadronQuadruple*>& v_hadronQuadruplesPN)
      {

	for(vector<HadronPair*>::iterator it=v_hadronPairsFirstH_PN.begin();it!=v_hadronPairsFirstH_PN.end();it++)
	  {
	    delete *it;
	  }
	for(vector<HadronPair*>::iterator it=v_hadronPairsSecondH_PN.begin();it!=v_hadronPairsSecondH_PN.end();it++)
	  {
	    delete *it;
	  }

	for(vector<HadronQuadruple*>::iterator it=v_hadronQuadruplesPN.begin();it!=v_hadronQuadruplesPN.end();it++)
	  {delete *it;}
	v_hadronQuadruplesPN.clear();
	//      cout <<"c1"<<endl;
	v_hadronPairsFirstH_PN.clear();

	//      cout <<"c2"<<endl;
	v_hadronPairsSecondH_PN.clear();

	v_firstHemiPos.clear();
	v_firstHemiNeg.clear();

	v_secondHemiPos.clear();
	v_secondHemiNeg.clear();

	for(int i=0;i<v_allParticles.size();i++)
	  {
	    delete v_allParticles[i];
	  }
	v_allParticles.clear();
      }


      //mirrors thrust computation for real particles
      void computeGenThrust()
      {
	numQuarks=0;
	vector<Hep3Vector> allParticlesBoosted;  
	vector<PseudoJet> fjParticles;
	//in principle easier to just boot the thrust vector I guess...
	vector<Hep3Vector> allParticlesNonBoosted; 
	vector<float> allPB_E;
	Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();
	//      	cout <<"-------------------------new evt--------------------------"<<endl;
	//	cout <<setw(9)<<" ID "<<setw(9)<<"ISTHEP"<<setw(9)<<" LUND_ID  "<<setw(9)<<"Mo first"<<setw(9)<<"    MoLast"<<setw(9)<<"   daFirst"<<setw(9)<<" daLast "<<setw(9)<<" PX  "<<setw(9)<<"PY"<<setw(9)<<"  PZ "<<setw(9)<<"  E "<<setw(9)<<"  M  "<<setw(9)<<" VX "<<setw(9)<<"  VY "<<setw(9)<<"  VZ"<<setw(9)<<"   T " <<endl;
	HepLorentzVector* lv_q1;
	HepLorentzVector* lv_q2;
	float m_px=0;
	float m_py=0;
	float m_pz=0;
	float m_e=0;

	float m_px2=0;
	float m_py2=0;
	float m_pz2=0;
	float m_e2=0;

	float m_pxVp=0;
	float m_pyVp=0;
	float m_pzVp=0;
	float m_eVp=0;


	HepLorentzVector boostVP(gen_hep_Mgr.begin()->PX(),gen_hep_Mgr.begin()->PY(),gen_hep_Mgr.begin()->PZ(),gen_hep_Mgr.begin()->E());
	boostVP.boost(kinematics::CMBoost);
	vpPx=boostVP.px();
	vpPy=boostVP.py();
	vpPz=boostVP.pz();
	vpEnergy=boostVP.e();

	float pxSum=0.0;
	float pySum=0.0;
	float pzSum=0.0;
	float eSum=0.0;

	for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
	  {

	    //	  cout << setw(6)<< gen_it->get_ID() << " |  "  << setw(6)<<gen_it->isthep() << " | " << setw(6)<<gen_it->idhep()  <<" | " <<  setw(6)<< gen_it->moFirst() <<" | " <<setw(6)<< gen_it->moLast() << " | " <<  setw(6)<<gen_it->daFirst() <<" | " <<setw(6)<< gen_it->daLast() <<" | " <<setw(6)<< gen_it->PX() << " | " << setw(6)<< gen_it->PY() <<" | " << setw(6)<< gen_it->PZ() <<" | " <<setw(6)<< gen_it->E() << " | " << setw(6)<< gen_it->M() <<" | " << setw(6)<< gen_it->VX() << " | " << setw(6)<< gen_it->VY() << " | " << setw(6)<< gen_it->VZ() << " | " << setw(6)<< gen_it->T() <<endl;


	    /*	    if(gen_it->idhep()>94 && gen_it->idhep() <100)
	      cout <<"found special code : " << gen_it->idhep() <<endl;
	      if(gen_it->get_ID()==1)
	      {
	      m_pxVp=gen_it->PX();
	      m_pyVp=gen_it->PY();
	      m_pzVp=gen_it->PZ();
	      m_eVp=gen_it->E();
	      }
	      ///if(gen_it->idhep()!=0 && gen_it->moFirst()==1 && (gen_it->idhep() < 20 || (gen_it->idhep()>20 && gen_it->idhep() <40)  ))//either lepton or gauge boson that came from gamma*
	      //	    	/*       if(gen_it->idhep()!=0 && gen_it->moFirst()==1 && (gen_it->idhep() > 110 || (gen_it->idhep()>1 && gen_it->idhep() <40)  ))//either lepton or gauge boson that came from gamma* */
	    /* 	    	    if(gen_it->idhep()!=0 && gen_it->moFirst()==1) */
	    /* 	      { */
	    /* 		if(gen_it->daLast()==-1 )//either lepton or gauge boson that came from gamma* */
	    /* 		  { */
	    /* 		    m_px+=gen_it->PX(); */
	    /* 		    m_py+=gen_it->PY(); */
	    /* 		    m_pz+=gen_it->PZ(); */
	    /* 		    m_e+=gen_it->E(); */
	    /* 		  } */
	    /* 		else */
	    /* 		  { */
	    /* 		    m_px2+=gen_it->PX(); */
	    /* 		    m_py2+=gen_it->PY(); */
	    /* 		    m_pz2+=gen_it->PZ(); */
	    /* 		    m_e2+=gen_it->E(); */
	    /* 		  } */
	    /* 		  }*\/ */


	    if(abs(gen_it->idhep()) <10 && gen_it->idhep()!=0) 
	      { 
		if(numQuarks==0) 
		  lv_q1=new HepLorentzVector(gen_it->PX(),gen_it->PY(),gen_it->PZ(), gen_it->E()); 
		if(numQuarks==1) 
		  { 
		    lv_q2=new HepLorentzVector(gen_it->PX(),gen_it->PY(),gen_it->PZ(), gen_it->E()); 
		    lv_q1->boost(kinematics::CMBoost); 
		    lv_q2->boost(kinematics::CMBoost); 
		  } 
		numQuarks++; 
		  
	      }
	    /* 	    /\*	    	    	    if(abs(gen_it->idhep()) <10 && gen_it->idhep()!=0) */
	    /* 	      { */
	    /* 		//		cout <<gen_it->idhep()<<endl; */

	    /* 		//		cout <<"daFirst: " << gen_it->daFirst()<< ", daSecond: " << gen_it->daLast() <<endl; */
	    /* 		//		cout <<"moFirst: " << gen_it->moFirst()<< ", moSecond: " << gen_it->moLast() <<endl; */
	    /* 		//		cout <<"mother: " << gen_it->mother() <<endl; */
	    /* 		if(gen_it->mother()) */
	    /* 		  { */
	    /* 		    //10022 is virtual photon */
	    /* 		    cout <<"motherID " << gen_it->mother().idhep() << " mother addr: " << gen_it->get_ID()<<endl; */
	    /* 		    cout <<"MotherE: " << gen_it->mother().E() <<", px : "<< gen_it->mother().PX() << " , py: " << gen_it->mother().PY() << " pz: " << gen_it->mother().PZ() <<endl; */
	    /*	    cout <<"mother addr: " << gen_it->mother() <<" dafirst : " << gen_it->mother().daFirst() << " last: " << gen_it->mother().daLast() <<endl;
	      if(gen_it->mother().mother()!=0)
	      cout <<"mother has a mother ;-) " <<endl;

	      cout <<"dafirst: " << gen_it->daFirst() <<" last: " << gen_it->daLast()<<endl;
	      }
	      cout <<"E: " << gen_it->E() <<", vx : "<< gen_it->VX() << " , vy: " << gen_it->VY() << " vz: " << gen_it->VZ() << " T: " <<gen_it->T() <<endl;
	      cout <<"px : "<< gen_it->PX() << " , py: " << gen_it->PY() << " pz: " << gen_it->PZ() <<endl;
	      cout <<"quarkCount: " <<quarkCount<< endl;
	      }*/
	    int gId=abs(gen_it->idhep());
	    //	    if(gId==ID_GAMMA || gId==ID_PI || gId==ID_K)
	    //	    if(gen_it->mother()==1 && gen_it->daLast()==-1)
	    //	    if(gen_it->isthep()==2 || gen_it->mother()==1 && gen_it->idhep()<30)
	    //	    if(gen_it->isthep()==1)
	    //	    if(gen_it->mother()==1 && (gen_it->idhep()<10 || (gen_it->idhep()>21)))// ||gen_it->daLast()==0))//gives around 124 mrad, with the ==0 it gives 167
	    //	    if(gen_it->mother()==1 && (gen_it->idhep()<10 || (gen_it->idhep()>21)) && gen_it->daLast()==-1)// ||gen_it->daLast()==0))//gives around 124 mrad, with the ==0 it gives 167



	    //	    if( gen_it->daLast()==-1)
	    
	    HepLorentzVector boostedVec(gen_it->PX(),gen_it->PY(),gen_it->PZ(),gen_it->E());
	    allParticlesNonBoosted.push_back(boostedVec.vect());
	    boostedVec.boost(kinematics::CMBoost);
	    //replace with the question if final state particle...

	    // cout <<"isthep is: " << gen_it->isthep()<<endl;
	    if(gen_it->isthep()==1) //isthep==1: stable, ==2, unstable
	      {
		//	      if(gen_it->idhep()>6 &&(gId==ID_GAMMA || gId==ID_PI || gId==ID_K))
		fjParticles.push_back(PseudoJet(boostedVec.px(),boostedVec.py(),boostedVec.pz(),boostedVec.e()));

		pxSum+=boostedVec.px();
		pySum+=boostedVec.py();
		pzSum+=boostedVec.pz();
		eSum+=boostedVec.e();

	      }
	    if(gen_it->mother()==1 && (gen_it->idhep()<10 || (gen_it->idhep()>20&& gen_it->idhep() < 26)) && gen_it->daLast()==-1)// ||gen_it->daLast()==0))//gives around 124 mrad, with the ==0 it gives 167
	      {
		allParticlesBoosted.push_back(boostedVec.vect());
		allPB_E.push_back(gen_it->E());
	      }

	  }

	//      cout <<" sum Px: "<< pxSum << " sum Py: " << pySum <<" sumPz: "<< pzSum << " sumPe: "<< eSum<<endl;
	//	cout << lv_q1->px() <<" " <<lv_q1->py() << " " << lv_q1->pz() <<" " << lv_q1->e() <<endl;
	//	cout << lv_q2->px() <<" " <<lv_q2->py() << " " << lv_q2->pz() <<" " << lv_q2->e() <<endl;
	//	cout <<"p Sum: " << m_px << " " << m_py << " " << m_pz << " " << m_e <<endl;
	//	cout <<"pVp Sum: " << m_pxVp << " " << m_pyVp << " " << m_pzVp << " " << m_eVp <<endl;
	//	    cout <<"num Quarks: " <<numQuarks<<endl;

	quarkAngle= lv_q1->angle(*lv_q2);

	HepLorentzVector boostP(m_px,m_py,m_pz,m_e);
	HepLorentzVector boostP2(m_px2,m_py2,m_pz2,m_e2);
	HepLorentzVector boostPVp(m_pxVp,m_pyVp,m_pzVp,m_eVp);

	boostP.boost(kinematics::CMBoost);
	boostP2.boost(kinematics::CMBoost);
	boostPVp.boost(kinematics::CMBoost);

	//	cout <<"p Sum Boosted: " << boostP.px() << " " << boostP.py() << " " << boostP.pz() << " " << boostP.e()<< " " << abs(boostP.pz())+boostP.e()<<endl;
	//	cout <<"p Sum Boosted2: " << boostP2.px() << " " << boostP2.py() << " " << boostP2.pz() << " " << boostP2.e()<< " " << abs(boostP2.pz())+boostP2.e()<<endl;
	//	cout <<"p SumVp Boosted: " << boostPVp.px() << " " << boostPVp.py() << " " << boostPVp.pz() << " " << boostPVp.e()<< " " << abs(boostPVp.pz())+boostPVp.e()<<endl;
	delete lv_q1;
	delete lv_q2;


	////////////////jet algos....

	JetDefinition jet_def(ee_genkt_algorithm,kinematics::R,-1);
	ClusterSequence cs(fjParticles,jet_def);
	vector<PseudoJet> jets=sorted_by_E(cs.inclusive_jets());
	int numHighEJets=0;
	for(unsigned int i=0;i<jets.size();i++)
	  {
	    if(jets[i].modp()>2.75)
	      numHighEJets++;
	  }

	//need dijet event...
	if(numHighEJets!=2)
	  {
	    //      return;
	  }
	if(kinematics::thrustZReverted)
	  {
	    jet2=Hep3Vector(jets[0].px(),jets[0].py(),jets[0].pz());
	    //to have the same convention as with thrust we have to flip the direction of the second jet...
	    jet1=Hep3Vector((-1)*jets[1].px(),(-1)*jets[1].py(),(-1)*jets[1].pz());
	  }
	else
	  {
	    jet1=Hep3Vector(jets[0].px(),jets[0].py(),jets[0].pz());
	    //to have the same convention as with thrust we have to flip the direction of the second jet...
	    jet2=Hep3Vector((-1)*jets[1].px(),(-1)*jets[1].py(),(-1)*jets[1].pz());
	  }

	fluffiness1=AuxFunc::computeFluffiness(jets[0]);
	fluffiness2=AuxFunc::computeFluffiness(jets[1]);

	jetE1=jets[0].modp();
	jetE2=jets[1].modp();

	jetPhi1=jets[0].phi();
	jetPhi2=jets[1].phi();

	jetTheta1=jets[0].phi();
	jetTheta2=jets[1].phi();


	////?
	dijet[0]=kinematics::jet1;
	dijet[1]=kinematics::jet2;


	////////////
	cmThrust=thrust(allParticlesBoosted.begin(),allParticlesBoosted.end(),gretSelf);
	Thrust t=thrustall(allParticlesBoosted.begin(),allParticlesBoosted.end(),gretSelf);
	cmThrustMag=t.thru;

	labThrust=thrust(allParticlesNonBoosted.begin(),allParticlesNonBoosted.end(),gretSelf);

	///well, randomizing the direciton is not a good idea, lets set it in the same dir as the original thrust...
	//      if(rand() % 100 <50)
	if(kinematics::thrustDirCM.angle(cmThrust)>pi/2)
	  {
	    cmThrust.setZ((-1)*cmThrust.z());
	    cmThrust.setY((-1)*cmThrust.y());
	    cmThrust.setX((-1)*cmThrust.x());
	  }
	thetaEThrust=cmThrust.angle(kinematics::firstElectronCM.vect());
	for(int i=0;i<allParticlesBoosted.size();i++)
	  {
	    float m_z=2*allPB_E[i]/kinematics::Q;
	    float theta=AuxFunc::getTheta(cmThrust,allParticlesBoosted[i]);
	    m_histos->hEFlowMC->Fill(theta,m_z);
	  }
      }
      DebugHistos* m_histos;
      //    jet1.setX(cmThrust.GetX());
      //    jet2=cmThrust;
  };



#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
