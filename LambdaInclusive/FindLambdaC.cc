#include "particle/Particle.h"
#include "FindLambdaC.h"
#include "./AnaConsts.h"
#include <vector>
#include "mdst/mdst.h"


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//void FindLambdaC.c::find(Mdst_vee2_Manager::iterator, Vind ipion){
int FindLambdaC::find(HepLorentzVector pLambda, Vint ipion){

int lambdac_cout=0;
m_lambdac.clear();
m_lcmode.clear();
m_charge.clear();
if(ipion.size()<1) return lambdac_cout;

Mdst_charged_Manager& Mdst_charged_Mgr=Mdst_charged_Manager::get_manager();
Mdst_charged_Manager::iterator chr_it=Mdst_charged_Mgr.begin();
Mdst_charged_Manager::iterator chr_it_begin=Mdst_charged_Mgr.begin();
 Mdst_pi0_Manager& Mdst_pi0_Mgr = Mdst_pi0_Manager::get_manager();

//find Xi first
 chr_it=Mdst_charged_Mgr.begin();

for(Mdst_pi0_Manager::iterator pi0_it=Mdst_pi0_Mgr.begin();pi0_it!=Mdst_pi0_Mgr.end();pi0_it++)
     {
    if((*pi0_it).chisq()>20) continue;
    HepLorentzVector p4pi0( (*pi0_it).px(), (*pi0_it).py(),(*pi0_it).pz(),(*pi0_it).energy());
    Particle Pi0(*pi0_it);
    Particle &child0 = Pi0.relation().child(0);
    Particle &child1 = Pi0.relation().child(1);
    HepLorentzVector ch0_fit = child0.p();
    HepLorentzVector ch1_fit = child1.p();
    double egam1=ch0_fit.e();
    double egam2=ch1_fit.e();
    if(p4pi0.m()<pi0mass_lo || p4pi0.m()>pi0mass_up) continue;
    HepLorentzVector p4lc =  p4pi0 + pLambda ;
    float mass = p4lc.m();
    if(mass<Ximass_lo||mass>Ximass_up)  continue;
    lambdac_cout++;
    //Particle Plambdac(p4lc, Ptype(charge>0 ? "D+":"D-"));
    Particle Plambdac(p4lc, Ptype("D0"));
    //Plambdac.relation().append(p);
    Plambdac.relation().append(Pi0);
    m_lambdac.push_back(Plambdac);
    m_lcmode.push_back(4);
    m_charge.push_back(0);
 }//end of loop pi0

//for(Mdst_charged_Manager::iterator chr_it=Mdst_charged_Mgr.begin();chr_it!=Mdst_charged_Mgr.end();chr_it++,ii++)
for(int ii=0;ii<ipion.size();ii++)
    { 
    chr_it= chr_it_begin + ipion[ii];
    const Mdst_charged &chg = *chr_it;
    Hep3Vector h3Vect((*chr_it).p(0),(*chr_it).p(1),(*chr_it).p(2));
    double E=sqrt(m_pi*m_pi+h3Vect.mag2());
    HepLorentzVector LabVec(h3Vect,E);
    int charge=(*chr_it).charge();
    //float mom= h3Vect.mag();
    //if(h3Vect.mag()<0.2) continue;
    Particle p(chg,charge>0 ? "PI+":"PI-"); 

    HepLorentzVector p4lc = LabVec + pLambda; 
    
    double openAngle = LabVec.vect().angle(pLambda.vect());
    //if(openAngle>60/180.*3.14159) continue;
    float mass = p4lc.m();
    //if(mass<lcmass_lo||mass>lcmass_up)  continue;
    if(mass<Ximass_lo||mass>lcmass_up)  continue;
    lambdac_cout++;
    Particle Plambdac(p4lc, Ptype(charge>0 ? "D+":"D-"));
    Plambdac.relation().append(p);
    m_lambdac.push_back(Plambdac);
    m_lcmode.push_back(0);
    m_charge.push_back(charge);
}

 chr_it=Mdst_charged_Mgr.begin();
 //Mdst_pi0_Manager& Mdst_pi0_Mgr = Mdst_pi0_Manager::get_manager();
for(int ii=0;ii<ipion.size();ii++)
   {
    chr_it= chr_it_begin + ipion[ii];
    const Mdst_charged &chg = *chr_it;
    Hep3Vector h3Vect((*chr_it).p(0),(*chr_it).p(1),(*chr_it).p(2));
    double E=sqrt(m_pi*m_pi+h3Vect.mag2());
    HepLorentzVector LabVec(h3Vect,E);
    //if(h3Vect.mag()<0.2) continue;
    int charge=(*chr_it).charge();
    Particle p(chg,charge>0 ? "PI+":"PI-");
 
 for(Mdst_pi0_Manager::iterator pi0_it=Mdst_pi0_Mgr.begin();pi0_it!=Mdst_pi0_Mgr.end();pi0_it++)
     {
    if((*pi0_it).chisq()>20) continue;
    HepLorentzVector p4pi0( (*pi0_it).px(), (*pi0_it).py(),(*pi0_it).pz(),(*pi0_it).energy());
    Particle Pi0(*pi0_it);
    Particle &child0 = Pi0.relation().child(0);
    Particle &child1 = Pi0.relation().child(1);
    HepLorentzVector ch0_fit = child0.p();
    HepLorentzVector ch1_fit = child1.p();
    double egam1=ch0_fit.e();
    double egam2=ch1_fit.e();
    //constrained fitted?
    //cout<<"pi0 mass "<<p4pi0.m()<<endl;
    if(p4pi0.m()<pi0mass_lo || p4pi0.m()>pi0mass_up) continue;

    //cout<<"ch0_fit + ch1_fit "<<(ch0_fit + ch1_fit).m()<<endl;
    //HepLorentzVector p4lc = LabVec + ch0_fit + ch1_fit + pLambda ;
    HepLorentzVector p4lc = LabVec +  p4pi0 + pLambda ;
    float mass = p4lc.m();
    if(mass<lcmass_lo||mass>lcmass_up)  continue;
    lambdac_cout++;

    //Particle Plambdac(p4lc, Ptype(charge>0 ? "lambdac_p":"lambdac_m"));
    Particle Plambdac(p4lc, Ptype(charge>0 ? "D+":"D-"));
    Plambdac.relation().append(p);
    Plambdac.relation().append(Pi0);
    m_lambdac.push_back(Plambdac);
    m_lcmode.push_back(1);
    m_charge.push_back(charge);    
     }//end of loop pi0

   }

if(ipion.size()<3) return lambdac_cout;
chr_it=Mdst_charged_Mgr.begin();
Mdst_charged_Manager::iterator chr_jt=Mdst_charged_Mgr.begin()+1;
Mdst_charged_Manager::iterator chr_tt=Mdst_charged_Mgr.begin()+2;
for(int ii=0;ii<ipion.size()-2;ii++)
    {
    chr_it= chr_it_begin + ipion[ii];
    const Mdst_charged &chg = *chr_it;
    Hep3Vector h3Vect((*chr_it).p(0),(*chr_it).p(1),(*chr_it).p(2));
    //if(h3Vect.mag()<0.2) continue;
    double E=sqrt(m_pi*m_pi+h3Vect.mag2());
    HepLorentzVector LabVec1(h3Vect,E);
    int charge1=(*chr_it).charge();
    Particle p1(chg,charge1>0 ? "PI+":"PI-");


   for(int jj=ii+1;jj<ipion.size()-1;jj++)
    {
    chr_jt= chr_it_begin + ipion[jj];
    const Mdst_charged &chg = *chr_jt;
    Hep3Vector h3Vect((*chr_jt).p(0),(*chr_jt).p(1),(*chr_jt).p(2));
    //if(h3Vect.mag()<0.2) continue;
    double E=sqrt(m_pi*m_pi+h3Vect.mag2());
    HepLorentzVector LabVec2(h3Vect,E);
    int charge2=(*chr_jt).charge();
    Particle p2(chg,charge2>0 ? "PI+":"PI-");

  for(int tt=jj+1;tt<ipion.size();tt++)
    {
    chr_tt= chr_it_begin + ipion[tt];
    int charge3=(*chr_tt).charge();
    if(charge3==charge1&&charge3==charge2) continue;
    const Mdst_charged &chg = *chr_tt;
    Hep3Vector h3Vect((*chr_tt).p(0),(*chr_tt).p(1),(*chr_tt).p(2));
    //if(h3Vect.mag()<0.2) continue;
    double E=sqrt(m_pi*m_pi+h3Vect.mag2());
    HepLorentzVector LabVec3(h3Vect,E);
    Particle p3(chg,charge3>0 ? "PI+":"PI-");

    HepLorentzVector p4lc = LabVec1+LabVec2 +LabVec3 + pLambda;
    float mass = p4lc.m();
    if(mass<lcmass_lo||mass>lcmass_up)  continue;
    lambdac_cout++;

    Particle Plambdac(p4lc, Ptype((charge1+charge2+charge3)>0 ? "D+":"D-"));
    Plambdac.relation().append(p1);
    Plambdac.relation().append(p2);
    Plambdac.relation().append(p3);
    m_lambdac.push_back(Plambdac);
    m_lcmode.push_back(2);
    m_charge.push_back(charge1+charge2+charge3);
   }
   }
   }

return lambdac_cout;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif


