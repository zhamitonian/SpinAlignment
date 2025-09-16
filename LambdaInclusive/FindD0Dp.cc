#include "particle/Particle.h"
#include "FindD0Dp.h"
#include "./AnaConsts.h"
#include <vector>
#include "mdst/mdst.h"
//#include <mdst/findLambda.h>
#include <mdst/findKs.h>


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//void FindD0Dp.c::find(Mdst_vee2_Manager::iterator, Vind ipion){
int FindD0Dp::find(Vint ipion, Vint ikaon, HepPoint3D ip_position){

//cout<<"FindD0Dp:: "<<ipion.size()<<", "<<ikaon.size()<<endl;

int D_cout=0;
m_Dmeson.clear();
m_Dmode.clear();
m_charge.clear();
//if(ipion.size()<1) return lambdac_cout;
Vp4 Vpion; Vpion.clear(); 
Vp4 Vkaon; Vkaon.clear();
Vp4 Vpi0; Vpi0.clear();
Vp4 VKs; VKs.clear();
Vint chpion; chpion.clear();
Vint chkaon; chkaon.clear();
VParticle vPion; vPion.clear();
VParticle vKaon; vKaon.clear();
VParticle vPi0; vPi0.clear();
VParticle vKs; vKs.clear();

Mdst_charged_Manager& Mdst_charged_Mgr=Mdst_charged_Manager::get_manager();
Mdst_charged_Manager::iterator chr_it=Mdst_charged_Mgr.begin();
Mdst_charged_Manager::iterator chr_it_begin=Mdst_charged_Mgr.begin();

//for(Mdst_charged_Manager::iterator chr_it=Mdst_charged_Mgr.begin();chr_it!=Mdst_charged_Mgr.end();chr_it++,ii++)
for(unsigned int ii=0;ii<ipion.size();ii++)
    { 
    chr_it= chr_it_begin + ipion[ii];
    const Mdst_charged &chg = *chr_it;
    Hep3Vector h3Vect((*chr_it).p(0),(*chr_it).p(1),(*chr_it).p(2));
    double E=sqrt(m_pi*m_pi+h3Vect.mag2());
    HepLorentzVector ppion(h3Vect,E);
    int picharge=(*chr_it).charge();
    Vpion.push_back(ppion);
    chpion.push_back(picharge);
    
    //float mom= h3Vect.mag();
    //if(h3Vect.mag()<0.2) continue;
    Particle pi(chg,picharge>0 ? "PI+":"PI-"); 
    vPion.push_back(pi);

}

  for(unsigned int jj=0;jj<ikaon.size();jj++)
  {
    chr_it= chr_it_begin + ikaon[jj];
    const Mdst_charged &chg = *chr_it;
    Hep3Vector h3Vect((*chr_it).p(0),(*chr_it).p(1),(*chr_it).p(2));
    double E=sqrt(m_k*m_k+h3Vect.mag2());
    HepLorentzVector pkaon(h3Vect,E);
    int kacharge=(*chr_it).charge();
    Vkaon.push_back(pkaon);
    chkaon.push_back(kacharge);

    Particle ka(chg,kacharge>0 ? "K+":"K-"); 
    vKaon.push_back(ka);

}
    
//pi0 list
int npi0=0;
 //chr_it=Mdst_charged_Mgr.begin();
 Mdst_pi0_Manager& Mdst_pi0_Mgr = Mdst_pi0_Manager::get_manager();
 
 for(Mdst_pi0_Manager::iterator pi0_it=Mdst_pi0_Mgr.begin();pi0_it!=Mdst_pi0_Mgr.end();pi0_it++)
     {
    if((*pi0_it).chisq()>20) continue;
    HepLorentzVector p4pi0( (*pi0_it).px(), (*pi0_it).py(),(*pi0_it).pz(),(*pi0_it).energy());
    Particle Pi0(*pi0_it);
    Particle &child0 = Pi0.relation().child(0);
    Particle &child1 = Pi0.relation().child(1);
    HepLorentzVector ch0_fit = child0.p();
    HepLorentzVector ch1_fit = child1.p();
    //double egam1=ch0_fit.e();
    //double egam2=ch1_fit.e();
   //cout<<"pi0 mass "<<p4pi0.m()<<endl;
    if(p4pi0.m()<pi0mass_lo || p4pi0.m()>pi0mass_up) continue;
   //cout<<"pi0 mass "<<p4pi0.m()<<endl;
   Vpi0.push_back(p4pi0);
   vPi0.push_back(Pi0);
   npi0++;
 }//end of loop pi0

//cout<<" pi0 size "<<Vpi0.size()<<endl;

//Ks list
int nKs=0;
Mdst_vee2_Manager& Mdst_vee2_Mgr = Mdst_vee2_Manager::get_manager();
FindKs findKs; int nLcan=0;
    Vint trkIndex_used; trkIndex_used.clear();
    for(Mdst_vee2_Manager::iterator vee_it=Mdst_vee2_Mgr.begin();vee_it!=Mdst_vee2_Mgr.end();vee_it++)
      {
    findKs.candidates((*vee_it),ip_position);
    int goodKs= findKs.goodKs();
    if(goodKs!=1) continue;
    HepLorentzVector p4( (*vee_it).px(), (*vee_it).py(),(*vee_it).pz(),(*vee_it).energy());
    VKs.push_back(p4);
    nKs++;
    Particle Ks(*vee_it);
    vKs.push_back(Ks);
    //Particle &child0 = KS.relation().child(0);
    //Particle &child1 = KS.relation().child(1);

    }

//cout<<" Ks size "<<VKs.size()<<endl;
//cout<<"FindD0Dp::size of pi "<<Vpion.size()<<endl;
//cout<<"FindD0Dp::size of ka "<<Vkaon.size()<<endl;
//cout<<"FindD0Dp::size of pi0 "<<Vpi0.size()<<endl;
//cout<<"FindD0Dp::size of Ks "<<VKs.size()<<endl;
//combind k/pi list
//D0--K pi 
for(unsigned int ipi=0;ipi<Vpion.size();ipi++ )
{
for(unsigned int ika=0;ika<Vkaon.size();ika++ )
{
if(chpion[ipi]*chkaon[ika]>0) continue;
HepLorentzVector Dp4 =  Vpion[ipi] + Vkaon[ika];
//cout<<" Dp4.m() "<< Dp4.m() <<" D0mass_lo  "<<D0mass_lo<<", D0mass_up "<<D0mass_up<<endl;
if(Dp4.m()<D0mass_lo||Dp4.m()>D0mass_up) continue;
D_cout++;

 Particle PD0(Dp4, Ptype("D0"));
 PD0.relation().append(vPion[ipi]);
 PD0.relation().append(vKaon[ika]);
 m_Dmeson.push_back(PD0);
 m_Dmode.push_back(0);
 m_charge.push_back(0);
}
}

//D0 --> K pi pi0
for(unsigned int ipi=0;ipi<Vpion.size();ipi++ )
{
for(unsigned int ika=0;ika<Vkaon.size();ika++ )
{
if(chpion[ipi]*chkaon[ika]>0) continue;
for(unsigned int ipi0=0;ipi0<Vpi0.size();ipi0++)
{
HepLorentzVector Dp4 =  Vpion[ipi] + Vkaon[ika] + Vpi0[ipi0];
if(Dp4.m()<D0mass_lo||Dp4.m()>D0mass_up) continue;
D_cout++;

 Particle PD0(Dp4, Ptype("D0"));
 PD0.relation().append(vPion[ipi]);
 PD0.relation().append(vKaon[ika]);
 PD0.relation().append(vPi0[ipi0]);
 m_Dmeson.push_back(PD0);
 m_Dmode.push_back(1);
 m_charge.push_back(0);
}
}
}

//D0 --> K 3pi
for(unsigned int ipi=0;ipi<Vpion.size();ipi++ )
{
for(unsigned int jpi=ipi+1;jpi<Vpion.size();jpi++ )
{
for(unsigned int tpi=jpi+1;tpi<Vpion.size();tpi++ )
{
if(chpion[ipi]==chpion[jpi]&&chpion[ipi]==chpion[tpi]) continue;
int sumCharge= chpion[ipi] + chpion[jpi] + chpion[tpi];
if(abs(sumCharge)!=1) cout<<"warning! charge of 3pi "<<sumCharge<<endl;
for(unsigned int ika=0;ika<Vkaon.size();ika++ )
{
if(sumCharge*chkaon[ika]>0) continue;
HepLorentzVector Dp4 =  Vpion[ipi] + Vpion[jpi] + Vpion[tpi] + Vkaon[ika];
if(Dp4.m()<D0mass_lo||Dp4.m()>D0mass_up) continue;
D_cout++;

 Particle PD0(Dp4, Ptype("D0"));
 PD0.relation().append(vPion[ipi]);
 PD0.relation().append(vPion[jpi]);
 PD0.relation().append(vPion[tpi]);
 PD0.relation().append(vKaon[ika]);
 m_Dmeson.push_back(PD0);
 m_Dmode.push_back(2);
 m_charge.push_back(0);
}
}
}
}


//Dp K+ pi- p-i
for(unsigned int ipi=0;ipi<Vpion.size();ipi++ )
{
for(unsigned int jpi=ipi+1;jpi<Vpion.size();jpi++ )
{
if(chpion[ipi]!=chpion[jpi]) continue;
//if(abs(sumCharge)!=1) cout<<"warning! charge of 3pi "<<sumCharge<<endl;
for(unsigned int ika=0;ika<Vkaon.size();ika++ )
{
if(chpion[ipi]*chkaon[ika]>0) continue;
HepLorentzVector Dp4 =  Vpion[ipi] + Vpion[jpi] + Vkaon[ika];
if(Dp4.m()<Dpmass_lo||Dp4.m()>Dpmass_up) continue;
int charge = chpion[ipi] + chpion[jpi] + chkaon[ika];
D_cout++;

 Particle PDp(Dp4, Ptype(charge>0 ? "D+":"D-"));
 PDp.relation().append(vPion[ipi]);
 PDp.relation().append(vPion[jpi]);
 PDp.relation().append(vKaon[ika]);
 m_Dmeson.push_back(PDp);
 m_Dmode.push_back(3);
 m_charge.push_back(charge);
}
}
}

//Dp K pi pi pi0
for(unsigned int ipi=0;ipi<Vpion.size();ipi++ )
{
for(unsigned int jpi=ipi+1;jpi<Vpion.size();jpi++ )
{
if(chpion[ipi]!=chpion[jpi]) continue;
//if(abs(sumCharge)!=1) cout<<"warning! charge of 3pi "<<sumCharge<<endl;
for(unsigned int ika=0;ika<Vkaon.size();ika++ )
{
if(chpion[ipi]*chkaon[ika]>0) continue;
for(unsigned int ipi0=0;ipi0<Vpi0.size();ipi0++)
{

HepLorentzVector Dp4 =  Vpion[ipi] + Vpion[jpi] + Vkaon[ika] + Vpi0[ipi0];
if(Dp4.m()<Dpmass_lo||Dp4.m()>Dpmass_up) continue;
int charge = chpion[ipi] + chpion[jpi] + chkaon[ika];
D_cout++;

 Particle PDp(Dp4, Ptype(charge>0 ? "D+":"D-"));
 PDp.relation().append(vPion[ipi]);
 PDp.relation().append(vPion[jpi]);
 PDp.relation().append(vKaon[ika]);
 PDp.relation().append(vPi0[ipi0]);
 m_Dmeson.push_back(PDp);
 m_Dmode.push_back(4);
 m_charge.push_back(charge);
}
}
}
}

//Dp --> Ks pi
for(unsigned int iks=0;iks<VKs.size();iks++)
{
for(unsigned int ipi=0;ipi<Vpion.size();ipi++ )
{

HepLorentzVector Dp4 =  Vpion[ipi] + VKs[iks];;
if(Dp4.m()<Dpmass_lo||Dp4.m()>Dpmass_up) continue;
D_cout++;
int charge = chpion[ipi];

 Particle PDp(Dp4, Ptype(charge>0 ? "D+":"D-"));
 PDp.relation().append(vPion[ipi]);
 PDp.relation().append(vKs[iks]);
 m_Dmeson.push_back(PDp);
 m_Dmode.push_back(5);
 m_charge.push_back(charge);
}
}

//Dp --> Ks pi pi0
for(unsigned int iks=0;iks<VKs.size();iks++)
{
for(unsigned int ipi=0;ipi<Vpion.size();ipi++ )
{
for(unsigned int ipi0=0;ipi0<Vpi0.size();ipi0++ )
{

HepLorentzVector Dp4 =  Vpion[ipi] + VKs[iks] + Vpi0[ipi0];
if(Dp4.m()<Dpmass_lo||Dp4.m()>Dpmass_up) continue;
D_cout++;
int charge = chpion[ipi];

 Particle PDp(Dp4, Ptype(charge>0 ? "D+":"D-"));
 PDp.relation().append(vPion[ipi]);
 PDp.relation().append(vKs[iks]);
 PDp.relation().append(vPi0[ipi0]);
 m_Dmeson.push_back(PDp);
 m_Dmode.push_back(6);
 m_charge.push_back(charge);
}
}
}

//Dp --> Ks 3pi 
for(unsigned int ipi=0;ipi<Vpion.size();ipi++ )
{
for(unsigned int jpi=ipi+1;jpi<Vpion.size();jpi++ )
{
for(unsigned int tpi=jpi+1;tpi<Vpion.size();tpi++ )
{
if(chpion[ipi]==chpion[jpi]&&chpion[ipi]==chpion[tpi]) continue;
int sumCharge= chpion[ipi] + chpion[jpi] + chpion[tpi];
if(abs(sumCharge)!=1) cout<<"warning! charge of 3pi "<<sumCharge<<endl;
for(unsigned int iks=0;iks<VKs.size();iks++ )
{
HepLorentzVector Dp4 =  Vpion[ipi] + Vpion[jpi] + Vpion[tpi] + VKs[iks];
int charge = chpion[ipi] + chpion[jpi] + chpion[tpi];
if(Dp4.m()<Dpmass_lo||Dp4.m()>Dpmass_up) continue;
D_cout++;

 Particle PDp(Dp4, Ptype(charge>0 ? "D+":"D-"));
 PDp.relation().append(vPion[ipi]);
 PDp.relation().append(vPion[jpi]);
 PDp.relation().append(vPion[tpi]);
 PDp.relation().append(vKs[iks]);
 m_Dmeson.push_back(PDp);
 m_Dmode.push_back(7);
 m_charge.push_back(charge);
}
}
}
}

//Dp --> K K pi 
for(unsigned int ika=0;ika<Vkaon.size();ika++ )
{
for(unsigned int jka=0;jka<Vkaon.size();jka++ )
{
if(chkaon[ika]==chkaon[jka]) continue;

for(unsigned int ipi=0;ipi<Vpion.size();ipi++ )
{
HepLorentzVector Dp4 =  Vkaon[ika] + Vkaon[jka] + Vpion[ipi];
if(Dp4.m()<Dpmass_lo||Dp4.m()>Dpmass_up) continue;
int charge = chkaon[ika] + chkaon[jka] + chpion[ipi];
D_cout++;

 Particle PDp(Dp4, Ptype(charge>0 ? "D+":"D-"));
 PDp.relation().append(vKaon[ika]);
 PDp.relation().append(vKaon[jka]);
 PDp.relation().append(vPion[ipi]);
 m_Dmeson.push_back(PDp);
 m_Dmode.push_back(8);
 m_charge.push_back(charge);

}
}
}

return D_cout;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif


