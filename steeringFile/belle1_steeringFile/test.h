//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep 19 11:14:50 2025 by ROOT version 6.26/14
// from TTree event/event
// found on file: basf_test.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class test {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           evtNo;
   Int_t           runNo;
   Int_t           expNo;
   Float_t         Q;
   Float_t         Energy_cms;
   Float_t         Evis_cms;
   Float_t         BalancePz_cms;
   Int_t           nGood;
   Int_t           nPip;
   Int_t           nPim;
   Int_t           nCluster;
   Float_t         thrust;
   Float_t         thrust_theta;
   Float_t         thrust_phi;
   Float_t         z;
   Float_t         pt;
   Float_t         cms_vecP[4];
   Int_t           nPhoton;
   Float_t         photon_p[19];   //[nPhoton]
   Float_t         photon_theta[19];   //[nPhoton]
   Float_t         photon_phi[19];   //[nPhoton]
   Int_t           nPip;
   Int_t           nPim;
   Float_t         pip_vecP[1][4];   //[nPip]
   Float_t         pim_vecP[1][4];   //[nPim]
   Float_t         pip_theta[1];   //[nPip]
   Float_t         pim_theta[1];   //[nPim]
   Float_t         pip_phi[1];   //[nPip]
   Float_t         pim_phi[1];   //[nPim]

   // List of branches
   TBranch        *b_evtNo;   //!
   TBranch        *b_runNo;   //!
   TBranch        *b_expNo;   //!
   TBranch        *b_Q;   //!
   TBranch        *b_Energy_cms;   //!
   TBranch        *b_Evis_cms;   //!
   TBranch        *b_BalancePz_cms;   //!
   TBranch        *b_nGood;   //!
   TBranch        *b_nPip;   //!
   TBranch        *b_nPim;   //!
   TBranch        *b_nCluster;   //!
   TBranch        *b_thrust;   //!
   TBranch        *b_thrust_theta;   //!
   TBranch        *b_thrust_phi;   //!
   TBranch        *b_z;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_cms_vecP;   //!
   TBranch        *b_nPhoton;   //!
   TBranch        *b_photon_p;   //!
   TBranch        *b_photon_theta;   //!
   TBranch        *b_photon_phi;   //!
   TBranch        *b_nPip;   //!
   TBranch        *b_nPim;   //!
   TBranch        *b_pip_vecP;   //!
   TBranch        *b_pim_vecP;   //!
   TBranch        *b_pip_theta;   //!
   TBranch        *b_pim_theta;   //!
   TBranch        *b_pip_phi;   //!
   TBranch        *b_pim_phi;   //!

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("basf_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("basf_test.root");
      }
      f->GetObject("event",tree);

   }
   Init(tree);
}

test::~test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void test::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evtNo", &evtNo, &b_evtNo);
   fChain->SetBranchAddress("runNo", &runNo, &b_runNo);
   fChain->SetBranchAddress("expNo", &expNo, &b_expNo);
   fChain->SetBranchAddress("Q", &Q, &b_Q);
   fChain->SetBranchAddress("Energy_cms", &Energy_cms, &b_Energy_cms);
   fChain->SetBranchAddress("Evis_cms", &Evis_cms, &b_Evis_cms);
   fChain->SetBranchAddress("BalancePz_cms", &BalancePz_cms, &b_BalancePz_cms);
   fChain->SetBranchAddress("nGood", &nGood, &b_nGood);
   fChain->SetBranchAddress("nPip", &nPip, &b_nPip);
   fChain->SetBranchAddress("nPim", &nPim, &b_nPim);
   fChain->SetBranchAddress("nCluster", &nCluster, &b_nCluster);
   fChain->SetBranchAddress("thrust", &thrust, &b_thrust);
   fChain->SetBranchAddress("thrust_theta", &thrust_theta, &b_thrust_theta);
   fChain->SetBranchAddress("thrust_phi", &thrust_phi, &b_thrust_phi);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("cms_vecP", cms_vecP, &b_cms_vecP);
   fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   fChain->SetBranchAddress("photon_p", photon_p, &b_photon_p);
   fChain->SetBranchAddress("photon_theta", photon_theta, &b_photon_theta);
   fChain->SetBranchAddress("photon_phi", photon_phi, &b_photon_phi);
//    fChain->SetBranchAddress("nPip", &nPip, &b_nPip);
//    fChain->SetBranchAddress("nPim", &nPim, &b_nPim);
   fChain->SetBranchAddress("pip_vecP", &pip_vecP, &b_pip_vecP);
   fChain->SetBranchAddress("pim_vecP", &pim_vecP, &b_pim_vecP);
   fChain->SetBranchAddress("pip_theta", &pip_theta, &b_pip_theta);
   fChain->SetBranchAddress("pim_theta", &pim_theta, &b_pim_theta);
   fChain->SetBranchAddress("pip_phi", &pip_phi, &b_pip_phi);
   fChain->SetBranchAddress("pim_phi", &pim_phi, &b_pim_phi);
   Notify();
}

Bool_t test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef test_cxx
