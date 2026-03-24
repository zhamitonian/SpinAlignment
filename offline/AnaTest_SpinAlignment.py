#!/usr/bin/env python3

import ROOT as R
from  OFFLINE_PROCESS import RDF_process, gMC_topoana, find_decay_indices
from DRAW import style_draw, HistStyle, Brush, graph_draw
from PHY_CALCULATOR import PhysicsCalculator
from math import sqrt,pi
from typing import Optional
import os
from typing import Tuple
from array import array
from functools import partial
from FIT import get_effCurve, perform_chisq_fit, QUICK_FIT, perform_resonance_fit
import re
from collections import defaultdict
from  phiSA_processor import SA
import pandas as pd

class Anatest_SpinAlignment:
    def __init__(self):
        self.tools = RDF_process()

    def check(self, rootFile = "temp.root"):
        file_b2bii = R.TFile(rootFile, 'READ')
        #tree= file_b2bii.Get("truth")
        tree= file_b2bii.Get("event")
        for i in range(5):
            tree.GetEntry(i)
            print("*", 25*"-", f"the {i}th event", 25*"-")
            for branch in tree.GetListOfBranches():
                name = branch.GetName()
                value = getattr(tree, name)
                if hasattr(value, '__len__') and not isinstance(value, str):
                    print(f"{name}: {[v for v in value]}")
                else:
                    print(f"{name}: {value}")

    def get_the_entry_index(self):
        file_b2bii = R.TFile("../steeringFile/b2bii_test.root", 'READ')
        tree = file_b2bii.Get('event_trk')
        for i in range(tree.GetEntries()):
            tree.GetEntry(i)
            if getattr(tree, "__event__") in [137, 573333, 4485]: 
                print(i)

    def little_test(self):
        path_old = "../steeringFile/belle1_steeringFile/exp71_rs2249_re2348_evtgen-uds_0_tree.root"
        path_new = "../steeringFile/belle1_steeringFile/test.root"
        df_old = R.RDataFrame("event", path_old)
        df_new = R.RDataFrame("event", path_new)
        h_old = df_old.Define("temp", "thrust[1]").Histo1D(("", ";cos#theta_{thrust};[]", 20, -1, 1), "temp").GetValue()
        h_new = df_new.Define("temp", "thrust[1]").Histo1D(("", ";cos#theta_{thrust};[]", 20, -1, 1), "temp").GetValue()

        def plot(h_b2bii, h_basf, label, var):
            
            c_all = R.TCanvas("c_all", "Combined", 1600, 1080)

            pad1 = R.TPad("pad1", "pad1", 0, 0.3, 1, 1)  
            pad1.SetBottomMargin(0.01)  
            pad1.SetTopMargin(0.1)
            pad1.SetLeftMargin(0.15)
            pad1.SetRightMargin(0.05)
            pad1.Draw()
            
            pad2 = R.TPad("pad2", "pad2", 0, 0, 1, 0.3)  
            pad2.SetTopMargin(0.02)     
            pad2.SetBottomMargin(0.3)
            pad2.SetLeftMargin(0.15)
            pad2.SetRightMargin(0.05)
            
            pad2.Draw()

            c_all.Update()

            style_draw([h_b2bii, h_basf], "", ["original","rotated"],styles = [HistStyle.line_hist(4, 1, 2), HistStyle.error_bars(1)], pad = pad1, save=False)

            h_ratio = h_basf.Clone("h_ratio")
            h_ratio.Sumw2()
            h_ratio.Divide(h_basf, h_b2bii, 1, 1, "B")
            h_ratio.GetYaxis().SetTitle("ratio")
            
            # 调整 pad2 (ratio plot) 的坐标轴标签大小和字体
            h_ratio.GetXaxis().SetLabelSize(0.10)   # X轴标签大小，因为pad2较小，需要更大
            h_ratio.GetYaxis().SetLabelSize(0.10)   # Y轴标签大小
            h_ratio.GetXaxis().SetTitleSize(0.16)   # X轴标题大小
            h_ratio.GetYaxis().SetTitleSize(0.14)   # Y轴标题大小
            h_ratio.GetXaxis().SetLabelFont(22)     # 字体 (42=Helvetica)
            h_ratio.GetYaxis().SetLabelFont(22)
            h_ratio.GetXaxis().SetTitleFont(22)
            h_ratio.GetYaxis().SetTitleFont(22)
            h_ratio.GetXaxis().SetTitleOffset(0.8)  # 标题偏移
            h_ratio.GetYaxis().SetTitleOffset(0.4)  # Y轴标题偏移

            #if var2 is not None:
                #var = var2.replace("[","_").replace("]","_") 
            style_draw([h_ratio], "./draw/b2bii_basf_comparison/bin_by_bin/{}.png".format(var), styles = [HistStyle.error_bars(1)], y_min= 0, y_max=2, use_user_y_range= True,pad= pad2)

            # diff
            h_diff = h_b2bii.Clone("h_diff")
            h_diff.Add(h_basf, -1)  # Subtract h1 from h0
            h_diff.SetTitle(label.split("/")[-1])
            style_draw([h_diff], "./draw/b2bii_basf_comparison/bin_by_bin/{}_diff.png".format(var),styles = [HistStyle.error_bars(1)])
        plot(h_old, h_new, ";cos#theta_{thrust};[]", "thrust[1]")

    def data_mc_comparsion_per_exp(self, output_dir, expNo):
        # 获取目录路径
        data_dir = "/gpfs/group/belle2/users2022/wangz/data_gMC_belle1/2025-10-28_SpinAlignmnet_Data/continuum/"
        gMC_dir = "/group/belle2/users2022/luruihua/for_wangz/data_gMC_belle1/2025-10-28_SpinAlignmnet_gMC/continuum/"

        exp_dir = os.path.join(output_dir, f"exp{expNo}")
        os.makedirs(exp_dir, exist_ok=True)
        
        def get_exp_df(path, expNo):
            path_list = []
            for file in os.listdir(path):
                if file.endswith(".root") and f"_e{expNo}_" in file:
                    full_path = os.path.join(path, file)
                    path_list.append(full_path)
            df = R.RDataFrame("event", path_list)
            df = SA(df)
            print(f"Number of files processed: {len(path_list)}")
            return df

        df_data = get_exp_df(data_dir, expNo)
        df_gMC = get_exp_df(gMC_dir, expNo) 

        def brush(h_mc, h_data, label, var):

            c_all = R.TCanvas("c_all", "Combined", 1600, 1080)

            pad1 = R.TPad("pad1", "pad1", 0, 0.3, 1, 1)  
            pad1.SetBottomMargin(0.01)  
            pad1.SetTopMargin(0.1)
            pad1.SetLeftMargin(0.15)
            pad1.SetRightMargin(0.05)
            pad1.Draw()
            
            pad2 = R.TPad("pad2", "pad2", 0, 0, 1, 0.3)  
            pad2.SetTopMargin(0.02)     
            pad2.SetBottomMargin(0.3)
            pad2.SetLeftMargin(0.15)
            pad2.SetRightMargin(0.05)
            
            pad2.Draw()

            c_all.Update()

            style_draw([h_mc, h_data], "", ["generic MC","data"],styles = [HistStyle.line_hist(4, 1, 2), HistStyle.error_bars(1)], pad = pad1, save=False)

            h_ratio = h_mc.Clone("h_ratio")
            h_ratio.Sumw2()
            h_ratio.Divide(h_mc, h_data, 1, 1, "B")
            h_ratio.GetYaxis().SetTitle("ratio")
            
            # 调整 pad2 (ratio plot) 的坐标轴标签大小和字体
            h_ratio.GetXaxis().SetLabelSize(0.10)   # X轴标签大小，因为pad2较小，需要更大
            h_ratio.GetYaxis().SetLabelSize(0.10)   # Y轴标签大小
            h_ratio.GetXaxis().SetTitleSize(0.16)   # X轴标题大小
            h_ratio.GetYaxis().SetTitleSize(0.14)   # Y轴标题大小
            h_ratio.GetXaxis().SetLabelFont(22)     # 字体 (42=Helvetica)
            h_ratio.GetYaxis().SetLabelFont(22)
            h_ratio.GetXaxis().SetTitleFont(22)
            h_ratio.GetYaxis().SetTitleFont(22)
            h_ratio.GetXaxis().SetTitleOffset(0.8)  # 标题偏移
            h_ratio.GetYaxis().SetTitleOffset(0.4)  # Y轴标题偏移

            png_name = f"{var.replace('[','_').replace(']','_')}"
            png_file = os.path.join(exp_dir, f"{png_name}.png")
            style_draw([h_ratio], png_file, styles = [HistStyle.error_bars(1)], y_min= 0, y_max=2, use_user_y_range= True,pad= pad2)

        def plot(label, nbin, xmin, xmax, var):
            def get_hist(df, var):
                if re.findall(r'\[.*?\]', var):
                    return df.Define("temp", var).Histo1D(("", label, nbin, xmin, xmax), "temp")
                return df.Histo1D(("", label, nbin, xmin, xmax), var)

            h_data = get_hist(df_data, var).GetValue()
            h_gMC = get_hist(df_gMC, var).GetValue()
            brush(h_gMC, h_data, label, var)

        plot(";N_{track};[]", 16, 2, 18, "nGood")
        plot(";N_{ECL};[]", 25, 0, 25,  "nECL") 
        plot(";N_{photon};[]", 25, 0, 25,  "nPhoton")
        plot(";E^{*}_{vis};[MeV]", 24, 1, 13, "Evis_cms")
        plot(";E^{*}_{sum};[MeV]", 24, 1, 13, "Esum_cms")
        plot(";#sum #vec P_{z}^{*};[MeV]", 24, -6, 6, "PzSum_cms")
        plot(";E_{vis};[MeV]", 24, 1, 13, "Evis")
        plot(";E_{sum};[MeV]", 24, 1, 9, "Esum")
        plot(";#sum #vec P_{z};[MeV]", 24, -6, 6, "PzSum")

        plot(";M_{heavy jet};[MeV]", 25, 0, 8, "HeavyJetMass")
        plot(";E_{heavy jet};[MeV]", 25, 0, 8, "HeavyJetEnergy")
        plot(";sphericity;[]", 25, 0, 1, "sphericity" )
        plot(";aplanarity;[]", 25, 0, 0.5, "aplanarity")
        plot(";foxWolframeR2;[]", 25, 0, 1, "foxWolfram[2]")
        plot(";cos#theta_{thrust};[]", 20, -1, 1, "thrust[1]")
        plot(";|#vec{T}|;[]", 51, 0.5, 1.01, "thrust[0]")
        plot(";#phi_{thrust};[rad]", 50, -pi, pi, "thrust[2]")

        # p theta phi for trk, cls, pho
        ## lab       
        plot(";p_{track};[MeV]", 50, 0, 5, "trk_p")
        plot(";#theta_{track};[rad]", 50, 0, pi, "trk_theta")
        plot(";#phi_{track};[rad]", 50, -pi, pi, "trk_phi")
        plot(";p_{cluster};[MeV]", 50, 0, 5, "cls_p")
        plot(";#theta_{cluster};[rad]", 50, 0, pi, "cls_theta")
        plot(";#phi_{cluster};[rad]", 50, -pi, pi, "cls_phi")
        plot(";p_{photon};[MeV]", 50, 0, 5, "pho_p")
        plot(";#theta_{photon};[rad]", 50, 0, pi, "pho_theta")
        plot(";#phi_{photon};[rad]", 50, -pi, pi, "pho_phi")

        ## CMS
        plot(";p_{track}^{*};[MeV]", 50, 0, 5, "trk_p_CMS")
        plot(";#theta_{track}^{*};[rad]", 50, 0, pi, "trk_theta_CMS")
        plot(";#phi_{track}^{*};[rad]", 50, -pi, pi, "trk_phi_CMS")
        plot(";p_{cluster}^{*};[MeV]", 50, 0, 5, "cls_p_CMS")
        plot(";#theta_{cluster}^{*};[rad]", 50, 0, pi, "cls_theta_CMS")
        plot(";#phi_{cluster}^{*};[rad]", 50, -pi, pi, "cls_phi_CMS")
        plot(";p_{photon}^{*};[MeV]", 50, 0, 5, "pho_p_CMS")
        plot(";#theta_{photon}^{*};[rad]", 50, 0, pi, "pho_theta_CMS")
        plot(";#phi_{photon}^{*};[rad]", 50, -pi, pi, "pho_phi_CMS")

        plot(";dr;[cm]", 20, 0, 1, "trk_dr") 
        plot(";|dz|;[cm]", 40, 0, 4, "trk_dz")


    def plot_nsig(self, txt_file, output_file):
        df_nsig = pd.read_csv(txt_file, sep=r"\s+",
                    names=["z_center", "z_width", "helicity_center", 
                            "helicity_width", "nsig", "nsig_err", "nsig_err2"],
                    skiprows=1)
    
        df_nsig = df_nsig[(df_nsig['nsig'] > 0)]

        hist_nsig = R.TH2D("hist_nsig", "nsig in cos#theta vs phi z;cos#theta*;phi z", 10, -1, 1, 20, 0, 1)

        for index, row in df_nsig.iterrows():
            binx = hist_nsig.GetXaxis().FindBin(row['helicity_center'])
            biny = hist_nsig.GetYaxis().FindBin(row['z_center'])

            hist_nsig.SetBinContent(binx, biny, row['nsig'])
            hist_nsig.SetBinError(binx, biny, row['nsig_err'])
        
        R.gStyle.SetOptStat(0)
        c = R.TCanvas("c_eff_nsig", "c_eff_nsig", 1600, 1080)
        R.gStyle.SetPaintTextFormat(".1f") 
        hist_nsig.Draw("COLZ TEXT")
        c.SaveAs(os.path.join("test_images", output_file))

    def test_qqbar_thrust(self):
        df_truth = R.RDataFrame("truth", "../steeringFile/belle1_steeringFile/test.root")
        h_thrust_cosTheta = df_truth.Define("temp", "abs(thrust_truth[1])").Histo1D(("", ";cos#theta_{thrust};[]", 20, 0, 1), "temp").GetValue()
        h_qqbar_cosTheta = df_truth.Define("temp", "abs(qqbar_axis[0])").Histo1D(("", ";cos#theta_{qqbar};[]", 20, 0, 1), "temp").GetValue()
        style_draw([h_thrust_cosTheta, h_qqbar_cosTheta], "./test_images/qqbarAxis_thrust_cosTheta.png", ["thrust axis","qqbar axis"], styles = [HistStyle.line_hist(4, 1, 2), HistStyle.error_bars(1)])

        h_thrust_phi = df_truth.Define("temp", "thrust_truth[2]").Histo1D(("", ";#phi_{thrust};[rad]", 20, -pi, pi), "temp").GetValue()
        h_qqbar_phi = df_truth.Define("temp", "qqbar_axis[1]").Histo1D(("", ";#phi_{qqbar};[rad]", 20, -pi, pi), "temp").GetValue()
        style_draw([h_thrust_phi, h_qqbar_phi], "./test_images/qqbarAxis_thrust_phi.png", ["thrust axis","qqbar axis"], styles = [HistStyle.line_hist(4, 1, 2), HistStyle.error_bars(1)])


        df_reco = R.RDataFrame("event", "../steeringFile/belle1_steeringFile/test.root")
        h_thrust_reco = df_reco.Define("temp", "thrust[0]").Histo1D(("", ";|#vec{T}|;[]", 51, 0.5, 1.01), "temp").GetValue()
        h_thrust_gen = df_reco.Define("temp", "thrust_truth[0]").Histo1D(("", ";|#vec{T}|;[]", 51, 0.5, 1.01), "temp").GetValue()
        style_draw([h_thrust_reco, h_thrust_gen], "./test_images/thrust_reco_gen.png", ["reco thrust","gen thrust"], styles = [HistStyle.line_hist(4, 1, 2), HistStyle.error_bars(1)])

    
    def check_xp_dist(self):
        truth_rootFile =  "/gpfs/group/belle/users/wangz/data_gMC_belle1/2025-12-03_SpinAlignment_qqbarMC/continuum_truth_processed_cutPt.root"
        #truth_rootFile = "continuum_truth_check.root"
        reco_truth_matched_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/2025-12-03_SpinAlignment_qqbarMC/continuum_reco_truth_matched.root"
        reco_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/2025-12-03_SpinAlignment_qqbarMC/continuum_reco_processed.root"
        df_reco = R.RDataFrame("event", reco_rootFile)
        df_truth = R.RDataFrame("truth", truth_rootFile)
        df_reco_truth_matched = R.RDataFrame("event", reco_truth_matched_rootFile)

        h_xp_reco_unmatch = df_reco.Histo1D(("", ";x_{p};[]", 50, 0, 1), "phi_xp").GetValue()
        h_xp_truth = df_truth.Histo1D(("", ";x_{p};[]", 50, 0, 1), "phi_xp").GetValue()
        h_xp_reco = df_reco_truth_matched.Histo1D(("", ";x_{p};[]", 50, 0, 1), "phi_xp").GetValue()
        h_xp_reco_gen = df_reco_truth_matched.Histo1D(("", ";x_{p};[]", 50, 0, 1), "phi_gen_xp").GetValue()

        h_xp_reco_unmatch.Scale(1/h_xp_reco_unmatch.Integral())
        h_xp_reco_gen.Scale(1/h_xp_reco_gen.Integral())
        h_xp_truth.Scale(1/h_xp_truth.Integral())
        h_xp_reco.Scale(1/h_xp_reco.Integral())
        style_draw([h_xp_reco, h_xp_reco_gen, h_xp_truth, h_xp_reco_unmatch], "./test_images/phi_xp_reco_gen_truth_matched.png", ["reco phi x_p","gen truth matched phi x_p","truth phi x_p", "reco_umatch"], styles = [HistStyle.line_hist(4, 1, 2), HistStyle.line_hist(2, 1, 2), HistStyle.error_bars(1)])
        style_draw([h_xp_reco, h_xp_reco_gen, h_xp_truth], "./test_images/phi_xp_reco_gen_truth.png", ["reco phi x_p","gen truth matched phi x_p","truth phi x_p"], styles = [HistStyle.line_hist(4, 1, 2), HistStyle.line_hist(2, 1, 2), HistStyle.error_bars(1)])
        style_draw([h_xp_truth, h_xp_reco], "./test_images/phi_xp_reco_truth_matched.png", ["truth phi x_p","reco truth matched phi x_p"], styles = [HistStyle.line_hist(4, 1, 2), HistStyle.error_bars(1)])
    
    def check_truth_calculation(self):
        truth_rootFile =  "/gpfs/group/belle/users/wangz/data_gMC_belle1/2025-12-03_SpinAlignment_qqbarMC/continuum_truth_processed_cutPt.root"
        truth_rootFile_x = "continuum_truth_check.root"
        df_truth = R.RDataFrame("truth", truth_rootFile)
        df_truth_x = R.RDataFrame("truth", truth_rootFile_x)

        h = df_truth.Histo1D(("", ";M_{#phi};[MeV]", 100, 0.99, 1.15), "phi_M").GetValue()
        h.Scale(1/h.Integral())
        h_x = df_truth_x.Histo1D(("", ";M_{#phi};[MeV]", 100, 0.99, 1.15), "phi_M").GetValue()
        h_x.Scale(1/h_x.Integral())
        style_draw([h, h_x], "./test_images/phi_M_truth_comparsion.png", ["not x","x"], 
                   styles = [HistStyle.filled_line_hist(6, 3019), HistStyle.filled_line_hist(7,3022)])
    
    def plot_nsig(self, txt_file, output_file="nsig.png"):
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        df = pd.read_csv(txt_file, sep='\s+', comment='#', 
                     names=['center', 'width', 'nsig', 'nsig_err', 'nsig_err2'])

        plt.figure(figsize=(10, 6))
        plt.errorbar(df['center'], df['nsig'], xerr=df['width'], yerr=df['nsig_err'], fmt='o', markersize=3, capsize=2)#, color='black')
        plt.yscale('log')
        plt.xlabel(r'$p_t$', fontsize=20)
        plt.ylabel(r'$N_{\mathrm{signal}}$', fontsize=20)
        plt.grid(True, alpha=0.3)
        plt.savefig(output_file)

        print(df)
        df.describe()

        """

        n_bins = len(df)
        bin_edges = np.zeros(n_bins + 1)

        for i in range(n_bins):
            bin_center = df.iloc[i]['center']
            bin_width = df.iloc[i]['width']
            bin_edges[i] = bin_center - bin_width
        
        bin_edges[n_bins] = df.iloc[n_bins-1]['center'] + df.iloc[n_bins-1]['width']

        hist = R.TH1D("", ";x_{p};[MeV]", n_bins, bin_edges)

        for i in range(n_bins):
            bin_idx = i + 1
            nsig = df.iloc[i]['nsig']
            nsig_err = df.iloc[i]['nsig_err']

            hist.SetBinContent(bin_idx, nsig)
            hist.SetBinError(bin_idx, nsig_err)
        
        style_draw([hist], output_file, styles = [HistStyle.error_bars(1)], log_y= True)
        """

    def Kshort_mass_distribution(self):
        rootFile = "/gpfs/home/belle2/wangz/Work/SpinAlignment/steeringFile/belle1_steeringFile/Ks_null_test/processed_test.root"
        df = R.RDataFrame("event", rootFile)
        model = R.RDF.TH1DModel("Ks_mass", ";M_{K_{S}^{0}};[MeV]", 100, 0.4, 0.6)
        h_mass_combine = df.Histo1D(model, "ks_m_combine")
        h_mass_read = df.Histo1D(model, "Ks_M")
        #h_mass_read = df.Histo1D(model, "ks_m_read")
        style_draw([h_mass_combine, h_mass_read], "./test_images/test.png", ["combine", "read"], styles=[HistStyle.line_hist(4, 1), HistStyle.line_hist(2,5)])

if __name__ == "__main__":
    ana = Anatest_SpinAlignment()
    ana.Kshort_mass_distribution()
    #ana.check()

    nsig_txt1 = "/gpfs/home/belle2/wangz/Work/SpinAlignment/offline/draw/Extract_rho00/images_test/fit_phi/nsig_results.txt"
    nsig_txt2 = "/gpfs/home/belle2/wangz/Work/SpinAlignment/offline/draw/Extract_rho00/images/fit_phi_test/nsig_results.txt"
    #ana.plot_nsig(nsig_txt1, "nsig_test1.png")
    #ana.plot_nsig(nsig_txt2, "nsig_test2.png")
    #ana.test_qqbar_thrust()
    #ana.check("../steeringFile/belle1_steeringFile/test_truth_processed.root")
    #ana.check("../steeringFile/belle1_steeringFile/test_reco_truth_matched.root")
    #ana.check_xp_dist()
    #ana.check_truth_calculation()
    #ana.plot_nsig("./test_images/nsig_results_binning_v2.txt", "nsig_binning_v2.png")
    #ana.plot_nsig("./test_images/nsig_results_pt_uniform_bining.txt", "pt_uniform_binning.png")

