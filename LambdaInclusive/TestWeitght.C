void TestWeight() 
{

  TRandom r1(100);
  TRandom r2(200);

  TH1F *h1 = new TH1F("h1","",100, -1, 1);
  TH1F *h2 = new TH1F("h2","",100, -1, 1);

    double decayPar = 0.6;
    TRandom r1;
    for(int ii=0;ii<10000;ii++)
   {
    double costheta_mc = r1.Uniform(-1, 1);
    h1->Fill(costheta_mc);
    double random = r2.Uniform(0, 1+0.05*decayPar);
    double var = 1+decayPar*0.05*costheta_mc;
    //cout<<"random: "<<random<<" var "<<var<<endl;
    if(random>var) { /*cout<<"random: "<<random<<" var "<<var<<endl;*/ continue; }
    h2->Fill(costheta_mc);
   }


TCanvas *c1 = new TCanvas("c1","",1000,800);
c1->Divide(2);
c1->cd(1);
h1->Draw();
c1->cd(2);
h2->Draw();

}
