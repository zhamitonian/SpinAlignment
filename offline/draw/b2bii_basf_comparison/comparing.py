#!/usr/python3
import ROOT as R
from DRAW import HistStyle, style_draw

import os
import re
from math import pi
from collections import defaultdict

b2bii_path = "/gpfs/group/belle2/users2022/wangz/other_rootFiles/SpinAlignment/b2bii_basf_comparing_rootProd/b2bii.root"
basf_path = "/gpfs/group/belle2/users2022/wangz/other_rootFiles/SpinAlignment/b2bii_basf_comparing_rootProd/basf/exp71_rs2249_re2348_evtgen-uds_0_tree.root"
basf_path_rotated = "/gpfs/group/belle2/users2022/wangz/other_rootFiles/SpinAlignment/b2bii_basf_comparing_rootProd/basf/after_rotate.root"  ## after coordinates rotation

def bin_by_bin(b2bii_path, basf_path, output_dir = "."):
    """
    Compare histograms from b2bii and basf ROOT files bin by bin.
    """
    df_b2bii = R.RDataFrame("event_trk", b2bii_path)
    df_basf = R.RDataFrame("event", basf_path)
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
        style_draw([h_ratio], png_file, styles = [HistStyle.error_bars(1)], y_min= 0, y_max=2, use_user_y_range= True,pad= pad2)

        # diff
        diff_png_file = os.path.join(output_dir, "bin_by_bin", f"{var}_diff.png")
        h_diff = h_b2bii.Clone("h_diff")
        h_diff.Add(h_basf, -1)  # Subtract h1 from h0
        h_diff.SetTitle(label.split("/")[-1])
        style_draw([h_diff], diff_png_file, ["b2bii - basf"],styles = [HistStyle.error_bars(1)])


    def plot_bin_by_bin(label, nbin, xmin, xmax, var, var2=None):
        def get_hist(df, var):
            if re.findall(r'\[.*?\]', var):
                return df.Define("temp", var).Histo1D(("", label, nbin, xmin, xmax), "temp")
            return df.Histo1D(("", label, nbin, xmin, xmax), var)

        h_b2bii = get_hist(df_b2bii, var).GetValue()
        h_basf = get_hist(df_basf, var2 if var2 else var).GetValue()
        brush(h_b2bii, h_basf, label, var)

    plot_bin_by_bin(";N_{track};[]", 16, 2, 18, "NosTrack", "nGood")
    plot_bin_by_bin(";N_{Cluster};[]", 25, 0, 25, "NosCluster", "nCluster") 
    plot_bin_by_bin(";N_{photon};[]", 25, 0, 25, "NosPhoton", "nPhoton")
    plot_bin_by_bin(";E_{vis};[MeV]", 24, 1, 13, "EvisCMS", "Evis_cms")
    plot_bin_by_bin(";#sum #vec P_{z};[MeV]", 24, -6, 6, "BalancePzCMS", "BalancePz_cms")
    plot_bin_by_bin(";Eenergy_{cms};[MeV]", 24, 1, 13, "EnergyCMS", "Energy_cms")
    plot_bin_by_bin(";|#vec{T}|;[]", 24, 0.5, 1.1, "thrust", "thrust[0]")
    #plot_bin_by_bin(";E_{cms};[MeV]", 36, 2, 11, "Ecms")
    plot_bin_by_bin(";M_{heavy jet};[MeV]", 25, 0, 8, "HeavyJetMass")
    plot_bin_by_bin(";E_{heavy jet};[MeV]", 25, 0, 8, "HeavyJetEnergy")
    plot_bin_by_bin(";sphericity;[]", 25, 0, 1, "sphericity" )
    plot_bin_by_bin(";aplanarity;[]", 25, 0, 0.5, "aplanarity")
    plot_bin_by_bin(";foxWolframeR2;[]", 25, 0, 1, "foxWolframR2", "foxWolfram[2]")
    plot_bin_by_bin(";cos#theta_{thrust};[]", 20, -1, 1, "thrustAxisCosTheta", "thrust[1]")


    ## trk , cls and pho 's p theta phi
    def plot_hist_list(tree1, tree2, var1, var2, label, nbin, xmin, xmax):
        h1 = R.TH1F("h1", label, nbin, xmin, xmax)
        h2 = R.TH1F("h2", label, nbin, xmin, xmax)
        tree1.Draw(f"{var1}>>h1", "", "goff")
        tree2.Draw(f"{var2}>>h2", "", "goff")
        
        brush(h1, h2, label, var2)

        h1.Delete()
        h2.Delete()

    file_basf = R.TFile(basf_path)
    tree_basf = file_basf.Get("event")
    file_b2bii = R.TFile(b2bii_path)
    tree_b2bii_trk = file_b2bii.Get("event_trk")
    tree_b2bii_cls = file_b2bii.Get("event_cls")
    tree_b2bii_pho = file_b2bii.Get("event_pho")
    plot_hist_list(tree_b2bii_trk, tree_basf, "p", "trk_p", ";p_{track};[MeV]", 50, 0, 5)    
    plot_hist_list(tree_b2bii_trk, tree_basf, "theta", "trk_theta", ";#theta_{track};[rad]", 50, 0, pi)    
    plot_hist_list(tree_b2bii_trk, tree_basf, "phi", "trk_phi", ";#phi_{track};[rad]", 50, -pi, pi)    

    plot_hist_list(tree_b2bii_cls, tree_basf, "p", "cls_p", ";p_{cluster};[MeV]", 50, 0, 5)    
    plot_hist_list(tree_b2bii_cls, tree_basf, "theta", "cls_theta", ";#theta_{cluster};[rad]", 50, 0, pi)    
    plot_hist_list(tree_b2bii_cls, tree_basf, "phi", "cls_phi", ";#phi_{cluster};[rad]", 50, -pi, pi)

    plot_hist_list(tree_b2bii_pho, tree_basf, "p", "pho_p", ";p_{photon};[MeV]", 50, 0, 5)       
    plot_hist_list(tree_b2bii_pho, tree_basf, "theta", "pho_theta", ";#theta_{photon};[rad]", 50, 0, pi)    
    plot_hist_list(tree_b2bii_pho, tree_basf, "phi", "pho_phi", ";#phi_{photon};[rad]", 50, -pi, pi)

    # ------------- cms -------------------     
    plot_hist_list(tree_b2bii_trk, tree_basf, "p_CMS", "trk_p_CMS", ";p_{track}^{cms};[MeV]", 50, 0, 5)    
    plot_hist_list(tree_b2bii_trk, tree_basf, "theta_CMS", "trk_theta_CMS", ";#theta_{track}^{cms};[rad]", 50, 0, pi)    
    plot_hist_list(tree_b2bii_trk, tree_basf, "phi_CMS", "trk_phi_CMS", ";#phi_{track}^{cms};[rad]", 50, -pi, pi)    

    plot_hist_list(tree_b2bii_cls, tree_basf, "p_CMS", "cls_p_CMS", ";p_{cluster}^{cms};[MeV]", 50, 0, 5)    
    plot_hist_list(tree_b2bii_cls, tree_basf, "theta_CMS", "cls_theta_CMS", ";#theta_{cluster}^{cms};[rad]", 50, 0, pi)    
    plot_hist_list(tree_b2bii_cls, tree_basf, "phi_CMS", "cls_phi_CMS", ";#phi_{cluster}^{cms};[rad]", 50, -pi, pi)

    plot_hist_list(tree_b2bii_pho, tree_basf, "p_CMS", "pho_p_CMS", ";p_{photon}^{cms};[MeV]", 50, 0, 5)       
    plot_hist_list(tree_b2bii_pho, tree_basf, "theta_CMS", "pho_theta_CMS", ";#theta_{photon}^{cms};[rad]", 50, 0, pi)    
    plot_hist_list(tree_b2bii_pho, tree_basf, "phi_CMS", "pho_phi_CMS", ";#phi_{photon}^{cms};[rad]", 50, -pi, pi)


    def relative_difference_before_and_after_boost(tree1, tree2, var1, var2, label, nbin, x_min, x_max, y_min = 1e-4, y_max = 10):
        h1 = R.TH1F("h1", label, nbin, x_min, x_max)
        h1_cms = R.TH1F("h1_cms", label, nbin, x_min, x_max)
        h2 = R.TH1F("h2", label, nbin, x_min, x_max)
        h2_cms = R.TH1F("h2_cms", label, nbin, x_min, x_max)
        tree1.Draw(f"{var1}>>h1", "", "goff")
        tree2.Draw(f"{var2}>>h2", "", "goff")
        tree1.Draw(f"{var1}_CMS>>h1_cms", "", "goff")
        tree2.Draw(f"{var2}_CMS>>h2_cms", "", "goff")

        h_diff_lab = h1.Clone("lab_diff")
        h_diff_lab.Add(h2, -1)  
        for i in range(1, h_diff_lab.GetNbinsX() + 1):
            h_diff_lab.SetBinContent(i, abs(h_diff_lab.GetBinContent(i)))  # 取绝对值
        h_diff_lab.Divide(h1)  

        h_diff_cms = h1_cms.Clone("cms_diff")
        h_diff_cms.Add(h2_cms, -1)
        for i in range(1, h_diff_cms.GetNbinsX() + 1):
            h_diff_cms.SetBinContent(i, abs(h_diff_cms.GetBinContent(i)))
        h_diff_cms.Divide(h1_cms)

        style_draw([h_diff_lab, h_diff_cms], f"{output_dir}/bin_by_bin/relative_{var2}_boost.png", ["before boost","after boost"],styles= [HistStyle.line_hist(4,1,2), HistStyle.line_hist(2,1,2)],log_y=True, y_min = y_min, y_max = y_max, use_user_y_range=True)

    relative_difference_before_and_after_boost(tree_b2bii_pho, tree_basf, "theta", "pho_theta", ";|h_{b2bii} - h_{basf}|/h_{b2bii} (#theta);[]", 50, 0, pi, 1e-4, 10)
    relative_difference_before_and_after_boost(tree_b2bii_pho, tree_basf, "phi", "pho_phi", ";|h_{b2bii} - h_{basf}|/h_{b2bii} (#phi);[]", 50, -pi, pi, 1e-5, 0.1)
    relative_difference_before_and_after_boost(tree_b2bii_pho, tree_basf, "p", "pho_p", ";|h_{b2bii} - h_{basf}|/h_{b2bii} (p);[]", 50, 0, 5, 1e-7, 0.1)

    relative_difference_before_and_after_boost(tree_b2bii_trk, tree_basf, "theta", "trk_theta", ";|h_{b2bii} - h_{basf}|/h_{b2bii} (#theta);[]", 50, 0, pi, 1e-4, 10)
    relative_difference_before_and_after_boost(tree_b2bii_trk, tree_basf, "phi", "trk_phi", ";|h_{b2bii} - h_{basf}|/h_{b2bii} (#phi);[]", 50, -pi, pi, 1e-5, 0.1)
    relative_difference_before_and_after_boost(tree_b2bii_trk, tree_basf, "p", "trk_p", ";|h_{b2bii} - h_{basf}|/h_{b2bii} (p);[]", 50, 0, 5, 1e-7, 0.1)


def evt_by_evt(b2bii_path, basf_path, output_dir = "."):
    """
    Compare variables from b2bii and basf ROOT files event by event.
    """
    file = R.TFile(basf_path)
    tree = file.Get("event")
    tree.AddFriend("event_trk", b2bii_path)
    df = R.RDataFrame(tree)

    def evt_by_evt_diff(label, nbin, xmin, xmax, var, log_y=True):
        h = df.Define("temp", var).Histo1D(("", label, nbin, xmin, xmax), "temp")
        xTitle =  h.GetXaxis().GetTitle()  + " difference"
        h.GetXaxis().SetTitle(xTitle)
        var_name = var.split(" - ")[0].strip().replace("[","_").replace("]","_").replace(".","_")
        output_png = os.path.join(output_dir, "evt_by_evt", f"evt_diff_{var_name}.png")
        style_draw([h], output_png, styles = [HistStyle.error_bars(1)], log_y= log_y)

    evt_by_evt_diff(";nGood;[]", 4, -2, 2, "NosTrack - nGood", True)
    evt_by_evt_diff(";nCluster;[]", 4, -2, 2, "NosCluster - nCluster")
    evt_by_evt_diff(";nPhoton;[]", 4, -2, 2, "NosPhoton - nPhoton")
    evt_by_evt_diff(";Evis_cms;[MeV]", 30, -3, 3, "EvisCMS - Evis_cms")
    evt_by_evt_diff(";BalancePz_cms;[MeV]", 30, -3, 3, "BalancePzCMS - BalancePz_cms")
    evt_by_evt_diff(";Energy_cms;[MeV]", 30, -3, 3, "EnergyCMS - Energy_cms")
    evt_by_evt_diff(";Thrust;[]", 40, -2, 2, "event_trk.thrust - thrust[0]")
    evt_by_evt_diff(";cos#theta_{thrust};[]", 40, -2, 2, "thrustAxisCosTheta - thrust[1]")
    ###evt_by_evt_diff(";Ecms;[MeV]", 48, -12, 12, "event_trk.Ecms - Ecms")
    evt_by_evt_diff(";HeavyJetMass;[GeV]", 50, -2, 2, "event_trk.HeavyJetMass - HeavyJetMass")
    evt_by_evt_diff(";HeavyJetEnergy;[GeV]", 50, -2, 2, "event_trk.HeavyJetEnergy - HeavyJetEnergy")
    evt_by_evt_diff(";sphericity;[]", 25, -0.5, 0.5, "event_trk.sphericity - sphericity" )
    evt_by_evt_diff(";aplanarity;[]", 25, -0.5, 0.5, "event_trk.aplanarity - aplanarity")
    evt_by_evt_diff(";foxWolframeR2;[]", 25, -0.5, 0.5, "event_trk.foxWolframR2 - foxWolfram[2]")


def find_the_event(b2bii_path, basf_path, diff_var = ("NosTrack", "nGood")):
    """
    Find the event with differences between two root files.
    """
    dict_event = defaultdict(int)

    file_b2bii = R.TFile(b2bii_path)
    tree_b2bii = file_b2bii.Get("event_trk")
    #tree_b2bii = file_b2bii.Get("event_pho")

    file_basf = R.TFile(basf_path)
    tree_basf = file_basf.Get("event")

    def check_event(tree, exp_name, run_name, evt_name, var_name):
        diffent_entries =[]
        for entry in range(tree.GetEntries()):
            tree.GetEntry(entry)
            (exp, run, evt)  = (getattr(tree, exp_name), getattr(tree, run_name), getattr(tree, evt_name))
            var_value = getattr(tree, var_name)
            if dict_event[(exp, run, evt)] != 0 and dict_event[(exp, run, evt)] !=var_value :
                print(f"Event with different var_value found: exp {exp}, run {run}, evt {evt}, var_value in tree1 {dict_event[(exp, run, evt)]}, in tree2 {var_value }")
                print(f"coordinate entry: {entry}")
                diffent_entries.append(entry)
            dict_event[(exp, run, evt)] =var_value 
        print("total event with different var_value: ", len(diffent_entries))
        return diffent_entries

    # --------- trk ----------------
    check_event(tree_basf, "expNo", "runNo", "evtNo", diff_var[1] )
    entries =  check_event(tree_b2bii, "__experiment__", "__run__", "__event__", diff_var[0] )
    # entry 71759

    # ---------- pho ---------------
    #check_event(tree_basf, "expNo", "runNo", "evtNo", "nPhoton" )
    #check_event(tree_b2bii, "__experiment__", "__run__", "__event__", "NosPhoton")

    def print_var(tree, EntryNo):
        """
        print all variables of the event, also the array variables
        """
        tree.GetEntry(EntryNo)
        for branch in tree.GetListOfBranches():
            name = branch.GetName()
            value = getattr(tree, name)
            if hasattr(value, '__len__') and not isinstance(value, str):
                print(f"{name}: {[v for v in value]}")
            else:
                print(f"{name}: {value}")
                
    print_var(tree_b2bii, entries[0])
    print_var(tree_basf, entries[0])



if __name__ == "__main__":
    #bin_by_bin(b2bii_path, basf_path)
    #evt_by_evt(b2bii_path, basf_path)
    find_the_event(b2bii_path, basf_path, ("NosTrack", "nGood"))
    pass