#!/usr/python3
import ROOT as R
from DRAW import HistStyle, style_draw
from OFFLINE_PROCESS import RDF_process

import os
import re
from math import pi
from collections import defaultdict

def bin_by_bin(b2bii_path, basf_path, output_dir = "."):
    """
    Compare histograms from b2bii and basf ROOT files bin by bin.
    """
    df_b2bii = R.RDataFrame("event_trk", b2bii_path)
    df_basf = R.RDataFrame("event", basf_path)

    ### define |cos(theta)_thrust| and pt relative to thrust axis
    #df_basf = df_basf.Define("abs_cosTheta_thrust", "abs(thrust[1])")
    #df_b2bii = df_b2bii.Define("abs_cosTheta_thrust", "abs(thrustAxisCosTheta)")

    count_b2bii = int(df_b2bii.Count().GetValue())
    count_basf = int(df_basf.Count().GetValue())
    print(f"event entries,\n b2bii: {count_b2bii}\nbasf: {count_basf}")

    def brush(h_b2bii, h_basf, label, var):

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

        style_draw([h_b2bii, h_basf], "", ["b2bii","basf"],styles = [HistStyle.line_hist(4, 1, 2), HistStyle.error_bars(1)], pad = pad1, save=False)

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

        png_file = os.path.join(output_dir, "bin_by_bin", f"{var}.png")
        style_draw([h_ratio], png_file, styles = [HistStyle.error_bars(1)], y_min= 0.9, y_max=1.1, use_user_y_range= True,pad= pad2)
        #style_draw([h_ratio], png_file, styles = [HistStyle.error_bars(1)], y_min= 0, y_max=2, use_user_y_range= True,pad= pad2)

        """
        # diff
        diff_png_file = os.path.join(output_dir, "bin_by_bin", f"{var}_diff.png")
        h_diff = h_b2bii.Clone("h_diff")
        h_diff.Add(h_basf, -1)  # Subtract h1 from h0
        h_diff.SetTitle(label.split("/")[-1])
        style_draw([h_diff], diff_png_file, ["b2bii - basf"],styles = [HistStyle.error_bars(1)])
        """

    def plot_bin_by_bin(label, nbin, xmin, xmax, var, var2=None):
        def get_hist(df, var):
            if re.findall(r'\[.*?\]', var):
                return df.Define("temp", var).Histo1D(("", label, nbin, xmin, xmax), "temp")
            return df.Histo1D(("", label, nbin, xmin, xmax), var)

        h_b2bii = get_hist(df_b2bii, var).GetValue()
        h_basf = get_hist(df_basf, var2 if var2 else var).GetValue()
        brush(h_b2bii, h_basf, label, var)

   #plot_bin_by_bin(";N_{track};[]", 16, 2, 18, "Ntrk")
   #plot_bin_by_bin(";N_{Cluster};[]", 25, 0, 25, "Ncls") 
   #plot_bin_by_bin(";E_{vis};[MeV]", 24, 1, 13, "Evis")
   #plot_bin_by_bin(";#sum #vec P_{z};[MeV]", 24, -6, 6, "Pz")
   #plot_bin_by_bin(";|#vec{T}|;[]", 24, 0.5, 1.1, "thrust", "thrust[0]")
   #plot_bin_by_bin(";M_{heavy jet};[MeV]", 25, 0, 8, "heavyJetMass","HeavyJetMass")
   ##plot_bin_by_bin(";foxWolframeR2;[]", 25, 0, 1, "foxWolframR2", "foxWolfram[")
   #plot_bin_by_bin(";cos#theta_{thrust};[]", 20, -1, 1, "thrustAxisCosTheta", "thrust[1]")
   #plot_bin_by_bin(";|cos#theta_{thrust}|;[1]", 25, 0, 1, "abs_cosTheta_thrust")


    ## trk , cls and pho 's p theta phi
    def plot_hist_list(tree1, tree2, var1, var2, label, nbin, xmin, xmax):
        h1 = R.TH1F("h1", label, nbin, xmin, xmax)
        h2 = R.TH1F("h2", label, nbin, xmin, xmax)
        tree1.Draw(f"{var1}>>h1", "", "goff")
        tree2.Draw(f"{var2}>>h2", "", "goff")
        #tree1.Draw(f"{var1}>>h1", "p > 0.105", "goff")
        #tree2.Draw(f"{var2}>>h2", "gam_p > 0.105", "goff")
        
        brush(h1, h2, label, var2)

        h1.Delete()
        h2.Delete()

    file_basf = R.TFile(basf_path)
    tree_basf = file_basf.Get("event")
    file_b2bii = R.TFile(b2bii_path)
    #tree_b2bii_trk = file_b2bii.Get("track")
    #tree_b2bii_gam = file_b2bii.Get("photon")
    tree_b2bii_gam = file_b2bii.Get("event_cls")
    #plot_hist_list(tree_b2bii_trk, tree_basf, "p", "trk_p", ";p_{track};[MeV]", 50, 0, 5)    
    #plot_hist_list(tree_b2bii_trk, tree_basf, "cosTheta", "trk_costheta", ";cos#theta_{track};[1]", 50, -1, 1)    
    #plot_hist_list(tree_b2bii_trk, tree_basf, "phi", "trk_phi", ";#phi_{track};[rad]", 50, -pi, pi)    
    #plot_hist_list(tree_b2bii_trk, tree_basf, "p_CMS", "trk_p_CMS", ";p_{track};[MeV]", 50, 0, 5)    
    #plot_hist_list(tree_b2bii_trk, tree_basf, "theta_CMS", "trk_costheta_CMS", ";cos#theta_{track};[1]", 50, -1, 1)    
    #plot_hist_list(tree_b2bii_trk, tree_basf, "phi_CMS", "trk_phi_CMS", ";#phi_{track};[rad]", 50, -pi, pi)    

    #plot_hist_list(tree_b2bii_gam, tree_basf, "p", "gam_p", ";p_{cluster};[MeV]", 50, 0, 5)    
    #plot_hist_list(tree_b2bii_gam, tree_basf, "cosTheta", "gam_costheta", ";cos#theta_{cluster};[1]", 50, -1, 1)    
    #plot_hist_list(tree_b2bii_gam, tree_basf, "phi", "gam_phi", ";#phi_{cluster};[rad]", 50, -pi, pi)
    #plot_hist_list(tree_b2bii_gam, tree_basf, "p_CMS", "gam_p_CMS", ";p_{cluster};[MeV]", 50, 0, 5)    
    #plot_hist_list(tree_b2bii_gam, tree_basf, "theta_CMS", "gam_costheta_CMS", ";cos#theta_{cluster};[1]", 50, -1, 1)    
    #plot_hist_list(tree_b2bii_gam, tree_basf, "phi_CMS", "gam_phi_CMS", ";#phi_{cluster};[rad]", 50, -pi, pi)
    plot_hist_list(tree_b2bii_gam, tree_basf, "phi_CMS", "gam_phi_CMS", ";#phi_{cluster};[rad]", 50, -pi, pi)


def evt_by_evt(b2bii_path, basf_path, output_dir = "."):
    """
    Compare variables from b2bii and basf ROOT files event by event.
    """
    file = R.TFile(basf_path)
    tree = file.Get("event")
    tree.AddFriend("track", b2bii_path)
    df = R.RDataFrame(tree)

    def evt_by_evt_diff(label, nbin, xmin, xmax, var, log_y=True):
        h = df.Define("temp", var).Histo1D(("", label, nbin, xmin, xmax), "temp")
        xTitle =  h.GetXaxis().GetTitle()  + " difference"
        h.GetXaxis().SetTitle(xTitle)
        var_name = var.split(" - ")[0].strip().replace("[","_").replace("]","_").replace(".","_")
        output_png = os.path.join(output_dir, "evt_by_evt", f"evt_diff_{var_name}.png")
        style_draw([h], output_png, styles = [HistStyle.error_bars(1)], log_y= log_y)

    evt_by_evt_diff(";nGood;[]", 4, -2, 2, "track.Ntrk - Ntrk", True)
    evt_by_evt_diff(";nCluster;[]", 4, -2, 2, "track.Ncls- Ncls")
    evt_by_evt_diff(";Evis;[MeV]", 30, -3, 3, "track.Evis - Evis")
    evt_by_evt_diff(";#sum P_{z};[MeV]", 30, -3, 3, "track.Pz - Pz")
    evt_by_evt_diff(";Thrust;[]", 40, -2, 2, "event_trk.thrust - thrust[0]")
    evt_by_evt_diff(";cos#theta_{thrust};[]", 40, -2, 2, "thrustAxisCosTheta - thrust[1]")
    evt_by_evt_diff(";HeavyJetMass;[GeV]", 50, -2, 2, "heavyJetMass - HeavyJetMass")


if __name__ == "__main__":
    b2bii_path = "with_Evis_cut/hadronic_b2bii.root"
    basf_path_rotated = "with_Evis_cut/hadronic_belle.root"  ## after coordinates rotation
    #bin_by_bin(b2bii_path, basf_path_rotated, "with_Evis_cut")

    b2bii_path = "without_Evis_cut/hadronic_b2bii.root"
    basf_path_rotated = "without_Evis_cut/hadronic_belle.root"

    b2bii_temp = "/gpfs/home/belle2/wangz/Work/SpinAlignment/steeringFile/hadronic_b2bii.root"
    basf_temp = "/gpfs/home/belle2/wangz/Work/SpinAlignment/steeringFile/belle1_steeringFile/hadronic_selection_test/hadronic_belle.root"
    #bin_by_bin(b2bii_path, basf_path_rotated, "without_Evis_cut")

    b2bii_path = "/gpfs/group/belle2/users2022/wangz/other_rootFiles/SpinAlignment/basf2_hadronic_sel_validation/b2bii_test.root"
    bin_by_bin(b2bii_path, "without_Evis_cut/hadronic_belle.root", "without_Evis_cut")

