#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TVector.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TVector3.h>
#include "ROOT/RVec.hxx"  // 需要包含这个头文件
using ROOT::VecOps::RVec;

void MCtruth_nsig() {
    //TChain* fChain0 = new TChain("MCtruth");
    TChain* fChain0 = new TChain("truth");
    
    //fChain0->Add("/besfs5/groups/tauqcd/luruihua/cross_check/root/MC/root/349000/*.root");//MC root
    //fChain0->Add("/besfs5/groups/tauqcd/luruihua/cross_check/root/MC/root/350800/*.root");//MC root
    //fChain0->Add("/besfs5/groups/tauqcd/luruihua/cross_check/root/MC/root/350967/*.root");//MC root
    //fChain0->Add("/besfs5/groups/tauqcd/luruihua/cross_check/root/MC/root/351039/*.root");//MC root
    //fChain0->Add("/besfs5/groups/tauqcd/luruihua/cross_check/root/MC/root/351458/*.root");//MC root
    
    fChain0->Add("/besfs5/groups/tauqcd/luruihua/truth_Ks.root");//MC root
    /*
    double m_chisq[40];
    int m_MC_n_ks;
    double m_MC_ks_px[40],m_MC_ks_py[40],m_MC_ks_pz[40],m_MC_ks_E[40];
    double m_MC_pip_px[40],m_MC_pip_py[40],m_MC_pip_pz[40],m_MC_pip_E[40];

    fChain0->SetBranchAddress("nksMC1", &m_MC_n_ks);
    fChain0->SetBranchAddress("kspxMC1", &m_MC_ks_px[0]);
    fChain0->SetBranchAddress("kspyMC1", &m_MC_ks_py[0]);
    fChain0->SetBranchAddress("kspzMC1", &m_MC_ks_pz[0]);
    fChain0->SetBranchAddress("ksEnMC1", &m_MC_ks_E[0]);
    fChain0->SetBranchAddress("pippx_ksMC1", &m_MC_pip_px[0]);
    fChain0->SetBranchAddress("pippy_ksMC1", &m_MC_pip_py[0]);
    fChain0->SetBranchAddress("pippz_ksMC1", &m_MC_pip_pz[0]);
    fChain0->SetBranchAddress("pipEn_ksMC1", &m_MC_pip_E[0]);
    */

    RVec<double>* pip_px = nullptr;
    RVec<double>* pip_py = nullptr;
    RVec<double>* pip_pz = nullptr;
    RVec<double>* pip_E = nullptr;
    RVec<double>* pim_px = nullptr;
    RVec<double>* pim_py = nullptr;
    RVec<double>* pim_pz = nullptr;
    RVec<double>* pim_E = nullptr;
    
    // 设置分支地址
    fChain0->SetBranchAddress("pip_px_cms_truth", &pip_px);
    fChain0->SetBranchAddress("pip_py_cms_truth", &pip_py);
    fChain0->SetBranchAddress("pip_pz_cms_truth", &pip_pz);
    fChain0->SetBranchAddress("pip_E_cms_truth", &pip_E);
    fChain0->SetBranchAddress("pim_px_cms_truth", &pim_px);
    fChain0->SetBranchAddress("pim_py_cms_truth", &pim_py);
    fChain0->SetBranchAddress("pim_pz_cms_truth", &pim_pz);
    fChain0->SetBranchAddress("pim_E_cms_truth", &pim_E);


    double nsig_truth[40][10] = {0}; // 初始化为0的二维数组
    TLorentzVector ks_mc, pip_mc,pim_mc;
    double P_KS, costheta;
    int m_MC_n_ks=1;

    // 遍历树中的事件
    for (int i = 0; i < fChain0->GetEntries(); i++) {
    //for (int i = 0; i < 100000; i++) {
        fChain0->GetEntry(i);
        if(i % 10000 == 0) {
            std::cout << "Processing event: " << i << std::endl;
        }
        if(m_MC_n_ks <= 0) continue;
        for (int j = 0; j < m_MC_n_ks; j++) {
            //ks_mc.SetPxPyPzE(m_MC_ks_px[j], m_MC_ks_py[j], m_MC_ks_pz[j], m_MC_ks_E[j]);
            //pip_mc.SetPxPyPzE(m_MC_pip_px[j], m_MC_pip_py[j], m_MC_pip_pz[j], m_MC_pip_E[j]);
            pip_mc.SetPxPyPzE((*pip_px)[j], (*pip_py)[j], (*pip_pz)[j], (*pip_E)[j]);
            pim_mc.SetPxPyPzE((*pim_px)[j], (*pim_py)[j], (*pim_pz)[j], (*pim_E)[j]);
            ks_mc = pip_mc + pim_mc;
            P_KS = ks_mc.P();
            TVector3 boostVector = -ks_mc.BoostVector();
            pip_mc.Boost(boostVector);
            costheta = pip_mc.Vect().Dot(ks_mc.Vect())/(pip_mc.Vect().Mag()*ks_mc.Vect().Mag());
            int bin_p = floor(P_KS / 0.1);
            int bin_cos = floor((costheta + 1) / 0.2);
            if(bin_p < 0 || bin_p >= 16) continue;
            nsig_truth[bin_p][bin_cos]++;
            }
        }

    // 打开输出文本文件
    std::ofstream outfile("nsig_truth_wangz.txt");
    if (!outfile.is_open()) {
        std::cerr << "Error: Cannot open output file nsig_truth_masscut.txt" << std::endl;
        return;
    }

    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 10; j++) {
            double value = nsig_truth[i][j];
            double error = sqrt(value); // 假设误差为统计误差（sqrt(N））
            // 保存到文本文件
            outfile << i << " " << j << " " << value << " " << error << std::endl;
        }
    }
    // 关闭文本文件
    outfile.close();

    std::cout << "values saved to nsig_truth.txt" << std::endl;
}
