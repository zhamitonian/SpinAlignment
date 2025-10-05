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

class Anatest_SpinAlignment:
    def __init__(self):
        self.tools = RDF_process()

    # bin by bin comparison
    def bin_by_bin_comparing(self):
        #df_b2bii = R.RDataFrame("event_trk", "../steeringFile/b2bii_test.root")
        df_b2bii = R.RDataFrame("event_trk", "../steeringFile/b2bii_test_thrust_cal.root")
        #df_basf =  R.RDataFrame("event", "../steeringFile/belle1_steeringFile/exp71_rs2249_re2348_evtgen-uds_0_tree.root")
        df_basf =  R.RDataFrame("event", "../steeringFile/belle1_steeringFile/basf_test_thrust_cal.root")
        print(df_b2bii.Count().GetValue(), df_basf.Count().GetValue())

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

            style_draw([h_b2bii, h_basf], "", ["b2bii","basf"],styles = [HistStyle.line_hist(4, 1, 2), HistStyle.error_bars(1)], pad = pad1, save=False)

            h_ratio = h_basf.Clone("h_ratio")
            h_ratio.Sumw2()
            h_ratio.Divide(h_basf, h_b2bii, 1, 1, "B")
            h_ratio.GetYaxis().SetTitle("ratio")
            
            # 调整 pad2 (ratio plot) 的坐标轴标签大小和字体
            h_ratio.GetXaxis().SetLabelSize(0.10)   # X轴标签大小，因为pad2较小，需要更大
            h_ratio.GetYaxis().SetLabelSize(0.10)   # Y轴标签大小
            h_ratio.GetXaxis().SetTitleSize(0.12)   # X轴标题大小
            h_ratio.GetYaxis().SetTitleSize(0.12)   # Y轴标题大小
            h_ratio.GetXaxis().SetLabelFont(22)     # 字体 (42=Helvetica)
            h_ratio.GetYaxis().SetLabelFont(22)
            h_ratio.GetXaxis().SetTitleFont(22)
            h_ratio.GetYaxis().SetTitleFont(22)
            h_ratio.GetXaxis().SetTitleOffset(1.0)  # 标题偏移
            h_ratio.GetYaxis().SetTitleOffset(0.4)  # Y轴标题偏移

            #if var2 is not None:
                #var = var2.replace("[","_").replace("]","_") 
            style_draw([h_ratio], "./draw/b2bii_basf_comparison/bin_by_bin/{}.png".format(var), styles = [HistStyle.error_bars(1)], y_min= 0, y_max=2, use_user_y_range= True,pad= pad2)

            # diff
            h_diff = h_b2bii.Clone("h_diff")
            h_diff.Add(h_basf, -1)  # Subtract h1 from h0
            h_diff.SetTitle(label.split("/")[-1])
            style_draw([h_diff], "./draw/b2bii_basf_comparison/bin_by_bin/{}_diff.png".format(var), ["b2bii - basf"],styles = [HistStyle.error_bars(1)])


        def plot_bin_by_bin(label, nbin, xmin, xmax, var, var2=None):
            def get_hist(df, var):
                if re.findall(r'\[.*?\]', var):
                    return df.Define("temp", var).Histo1D(("", label, nbin, xmin, xmax), "temp")
                return df.Histo1D(("", label, nbin, xmin, xmax), var)

            h_b2bii = get_hist(df_b2bii, var).GetValue()
            h_basf = get_hist(df_basf, var2 if var2 else var).GetValue()
            plot(h_b2bii, h_basf, label, var)

           
        #lot_bin_by_bin(";N_{track};[]", 16, 2, 18, "NosTrack", "nGood")
        #lot_bin_by_bin(";N_{Cluster};[]", 25, 0, 25, "NosCluster", "nCluster") 
        #lot_bin_by_bin(";N_{photon};[]", 25, 0, 25, "NosPhoton", "nPhoton")
        #lot_bin_by_bin(";E_{vis};[MeV]", 24, 1, 13, "EvisCMS", "Evis_cms")
        #lot_bin_by_bin(";#sum #vec P_{z};[MeV]", 24, -6, 6, "BalancePzCMS", "BalancePz_cms")
        #lot_bin_by_bin(";Eenergy_{cms};[MeV]", 24, 1, 13, "EnergyCMS", "Energy_cms")
        plot_bin_by_bin(";|#vec{T}|;[]", 24, 0.5, 1.1, "thrust", "thrust[0]")
        ###plot_bin_by_bin(";E_{cms};[MeV]", 36, 2, 11, "Ecms")
        #lot_bin_by_bin(";M_{heavy jet};[MeV]", 25, 0, 8, "HeavyJetMass")
        #lot_bin_by_bin(";M_{heavy jet};[MeV]", 25, 0, 8, "HeavyJetEnergy")
        #lot_bin_by_bin(";sphericity;[]", 25, 0, 1, "sphericity" )
        #lot_bin_by_bin(";aplanarity;[]", 25, 0, 0.5, "aplanarity")
        #lot_bin_by_bin(";foxWolframeR2;[]", 25, 0, 1, "foxWolframR2", "foxWolfram[2]")
        plot_bin_by_bin(";cos#theta_{thrust};[]", 20, -1, 1, "thrustAxisCosTheta", "thrust[1]")

    
        ## trk , cls and pho 's p theta phi
        def plot_hist_list(tree1, tree2, var1, var2, label, nbin, xmin, xmax):
            h1 = R.TH1F("h1", label, nbin, xmin, xmax)
            h2 = R.TH1F("h2", label, nbin, xmin, xmax)
            tree1.Draw(f"{var1}>>h1", "", "goff")
            tree2.Draw(f"{var2}>>h2", "", "goff")
            
            plot(h1, h2, label, var2)

            h1.Delete()
            h2.Delete()

        file_basf = R.TFile("../steeringFile/belle1_steeringFile/exp71_rs2249_re2348_evtgen-uds_0_tree.root")
        tree_basf = file_basf.Get("event")
        file_b2bii = R.TFile("../steeringFile/b2bii_test.root")
        tree_b2bii_trk = file_b2bii.Get("event_trk")
        tree_b2bii_cls = file_b2bii.Get("event_cls")
        tree_b2bii_pho = file_b2bii.Get("event_pho")
        #plot_hist_list(tree_b2bii_trk, tree_basf, "p", "trk_p", ";p_{track};[MeV]", 50, 0, 5)    
        #plot_hist_list(tree_b2bii_trk, tree_basf, "theta", "trk_theta", ";p_{track};[rad]", 50, 0, pi)    
        #plot_hist_list(tree_b2bii_trk, tree_basf, "phi", "trk_phi", ";p_{track};[rad]", 50, -pi, pi)    

        #plot_hist_list(tree_b2bii_cls, tree_basf, "p", "cls_p", ";p_{cluster};[MeV]", 50, 0, 5)    
        #plot_hist_list(tree_b2bii_cls, tree_basf, "theta", "cls_theta", ";theta_{cluster};[rad]", 50, 0, pi)    
        #plot_hist_list(tree_b2bii_cls, tree_basf, "phi", "cls_phi", ";phi_{cluster};[rad]", 50, -pi, pi)

        #plot_hist_list(tree_b2bii_pho, tree_basf, "p", "pho_p", ";p_{photon};[MeV]", 50, 0, 5)       
        #plot_hist_list(tree_b2bii_pho, tree_basf, "theta", "pho_theta", ";theta_{photon};[rad]", 50, 0, pi)    
        #plot_hist_list(tree_b2bii_pho, tree_basf, "phi", "pho_phi", ";phi_{photon};[rad]", 50, -pi, pi)
            
    
    def event_by_event_comparing(self):
        #file = R.TFile('../steeringFile/belle1_steeringFile/exp71_rs2249_re2348_evtgen-uds_0_tree.root', 'READ')
        file = R.TFile('../steeringFile/belle1_steeringFile/basf_test_thrust_cal.root', 'READ')
        tree= file.Get("event")
        #tree.AddFriend("event_trk", "../steeringFile/b2bii_test.root")
        tree.AddFriend("event_trk",  "../steeringFile/b2bii_test_thrust_cal.root")
        df = R.RDataFrame(tree)

        #print(df.GetColumnNames())

        def evt_by_evt_diff(label, nbin, xmin, xmax, var, log_y=True):
            h = df.Define("temp", var).Histo1D(("", label, nbin, xmin, xmax), "temp")
            var_name = var.split(" - ")[0].strip().replace("[","_").replace("]","_").replace(".","_")
            style_draw([h], f"./draw/b2bii_basf_comparison/evt_by_evt/evt_diff_{var_name}.png", [var], styles = [HistStyle.error_bars(1)], log_y= log_y)

        #evt_by_evt_diff(";nGood;[]", 4, -2, 2, "NosTrack - nGood", True)
        #evt_by_evt_diff(";nCluster;[]", 4, -2, 2, "NosCluster - nCluster")
        #evt_by_evt_diff(";nPhoton;[]", 4, -2, 2, "NosPhoton - nPhoton")
        #evt_by_evt_diff(";Evis_cms;[MeV]", 30, -3, 3, "EvisCMS - Evis_cms")
        #evt_by_evt_diff(";BalancePz_cms;[MeV]", 30, -3, 3, "BalancePzCMS - BalancePz_cms")
        #evt_by_evt_diff(";Energy_cms;[MeV]", 30, -3, 3, "EnergyCMS - Energy_cms")
        evt_by_evt_diff(";Thrust;[]", 40, -2, 2, "event_trk.thrust - thrust[0]")
        evt_by_evt_diff(";cos#theta_{thrust};[]", 40, -2, 2, "thrustAxisCosTheta - thrust[1]")
        ###evt_by_evt_diff(";Ecms;[MeV]", 48, -12, 12, "event_trk.Ecms - Ecms")
        #evt_by_evt_diff(";HeavyJetMass;[GeV]", 50, -2, 2, "event_trk.HeavyJetMass - HeavyJetMass")
        #evt_by_evt_diff(";HeavyJetEnergy;[GeV]", 50, -2, 2, "event_trk.HeavyJetEnergy - HeavyJetEnergy")
        #evt_by_evt_diff(";sphericity;[]", 25, -0.5, 0.5, "event_trk.sphericity - sphericity" )
        #evt_by_evt_diff(";aplanarity;[]", 25, -0.5, 0.5, "event_trk.aplanarity - aplanarity")
        #evt_by_evt_diff(";foxWolframeR2;[]", 25, -0.5, 0.5, "event_trk.foxWolframR2 - foxWolfram[2]")


    def test(self):
        
        df =  R.RDataFrame("event", "../steeringFile/belle1_steeringFile/basf_test.root")
        h0 = df.Define("temp", "thrust[0]").Histo1D(("h_0", ";;", 25, 0.5, 1.1), "temp")
        h1 = df.Histo1D(("h_0", ";;", 25, 0.5, 1.1), "Thrust_test")
        style_draw([h0, h1], "./test_images/n.png", ["0","1"],styles = [HistStyle.line_hist(4,1,2), HistStyle.error_bars(1)])

        h0 = df.Histo1D(("h_0", ";;", 25, 0, 25), "nCluster")
        h1 = df.Histo1D(("h_0", ";;", 25, 0, 25), "nCluster_test")
        style_draw([h0, h1], "./test_images/nCluster.png", ["0","1"],styles = [HistStyle.line_hist(4,1,2), HistStyle.error_bars(1)])

        h0 = df.Histo1D(("h_0", ";;", 25, 0, 25), "nGood")
        h1 = df.Histo1D(("h_0", ";;", 25, 0, 25), "nGood_test")
        style_draw([h0, h1], "./test_images/nGood.png", ["0","1"],styles = [HistStyle.line_hist(4,1,2), HistStyle.error_bars(1)])

        h0 = df.Histo1D(("h_0", ";;", 25, 0, 8), "HeavyJetMass")
        h1 = df.Histo1D(("h_0", ";;", 25, 0, 8), "HeavyJetMass_test")
        style_draw([h0, h1], "./test_images/heavyJetMass.png", ["0","1"],styles = [HistStyle.line_hist(4,1,2), HistStyle.error_bars(1)])    
    
    def test2(self):
        df_b2bii = R.RDataFrame("event", "../steeringFile/b2bii_test.root")
        #df_b2bii = R.RDataFrame("event", "../steeringFile/test.root")
        h0 = df_b2bii.Histo1D(("h_0", ";;", 24, 1, 12), "EvisCMS")
        h1 = df_b2bii.Histo1D(("h_1", ";;", 24, 1, 12), "EnergyCMS")
        style_draw([h0, h1], "./test_images/Evis_Energy.png", ["EvisCMS","EnergyCMS"],styles = [HistStyle.line_hist(4,1,2), HistStyle.error_bars(1)])

        h0 = df_b2bii.Histo1D(("h_0", ";;", 25, 0, 25), "NosCluster").GetValue()
        h1 = df_b2bii.Histo1D(("h_1", ";;", 25, 0, 25), "NosPhoton").GetValue()
        style_draw([h0, h1], "./test_images/Ncls_Nphoton.png", ["NosCluster","NosPhoton"],styles = [HistStyle.line_hist(4,1,2), HistStyle.error_bars(1)])

        h0.Add(h1, -1)  # Subtract h1 from h0
        style_draw([h0], "./test_images/Ncls_minus_Nphoton.png", ["NosCluster - NosPhoton"],styles = [HistStyle.error_bars(1)])

        df_basf = R.RDataFrame("event", "../steeringFile/belle1_steeringFile/exp71_rs2249_re2348_evtgen-uds_0_tree.root")
        #df_basf = R.RDataFrame("event", "../steeringFile/belle1_steeringFile/basf_test.root")
        h0 = df_basf.Histo1D(("h_0", ";;", 24, 1, 12), "Evis_cms")
        h1 = df_basf.Histo1D(("h_1", ";;", 24, 1, 12), "Energy_cms")
        style_draw([h0, h1], "./test_images/Evis_Energy_basf.png", ["Evis_cms","Energy_cms"],styles = [HistStyle.line_hist(4,1,2), HistStyle.error_bars(1)])

        h0 = df_basf.Histo1D(("h_0", ";;", 25, 0, 25), "nCluster").GetValue()
        h1 = df_basf.Histo1D(("h_1", ";;", 25, 0, 25), "nPhoton").GetValue()
        #h0.Add(h1, -1)  # Subtract h1 from h0
        #style_draw([h0], "./test_images/Ncls_minus_Nphoton_basf.png", ["nCluster - nPhoton"],styles = [HistStyle.error_bars(1)])

        h2 = df_b2bii.Histo1D(("h_1", ";;", 25, 0, 25), "NosCluster").GetValue()
        h2.Add(h0, -1)  # Subtract h1 from h0
        style_draw([h2], "./test_images/Ncls_minus_Nphoton_basf.png", ["nCluster - nPhoton"],styles = [HistStyle.error_bars(1)])

    def test3(self):
        # test VariablesToEventBasedTree, okay
        df_1 = R.RDataFrame("event", "../steeringFile/test.root")
        df_2 = R.RDataFrame("event", "../steeringFile/test1.root")
        def plot(label, nbin, xmin, xmax, var, var2=None):
            def get_hist(df, var):
                if re.findall(r'\[.*?\]', var):
                    return df.Define("temp", var).Histo1D(("", label, nbin, xmin, xmax), "temp")
                return df.Histo1D(("", label, nbin, xmin, xmax), var)

            h_b2bii = get_hist(df_1, var).GetValue()
            h_basf = get_hist(df_2, var2 if var2 else var).GetValue()

            if var2 is not None:
                var = var2.replace("[","_").replace("]","_") 
            style_draw([h_b2bii, h_basf], "./test_images/{}.png".format(var), ["b2bii","basf"],styles = [HistStyle.line_hist(4,1,2), HistStyle.error_bars(1)])

            h_b2bii.Add(h_basf, -1)  # Subtract h1 from h0
            h_b2bii.SetTitle(label.split("/")[-1])
            style_draw([h_b2bii], "./test_images/{}_diff.png".format(var), ["b2bii - basf"],styles = [HistStyle.error_bars(1)])

        #plot(";N_{track};", 16, 2, 18, "NosTrack")
        #plot(";N_{Cluster};", 25, 0, 25, "NosCluster") 
        plot(";N_{photon};", 25, 0, 25, "NosPhoton" )
        plot(";E_{vis};MeV", 24, 1, 13, "EvisCMS")
        plot(";#sum #vec P_{z};MeV", 24, -6, 6, "BalancePzCMS")
        #plot(";Eenergy_{cms};MeV", 24, 1, 13, "EnergyCMS")
        #plot(";|#vec{T}|;", 24, 0.5, 1.1, "thrust")
        #plot(";E_{cms};MeV", 36, 2, 11, "Ecms")
        #plot(";M_{heavy jet};MeV", 25, 0, 8, "HeavyJetMass")
        #plot(";M_{heavy jet};MeV", 25, 0, 8, "HeavyJetEnergy")
        #plot(";sphericity;", 25, 0, 1, "sphericity" )
        #plot(";aplanarity;", 25, 0, 0.5, "aplanarity")
        #plot(";foxWolframeR2;", 25, 0, 1, "foxWolframR2")
        #plot(";cos#theta_{thrust};", 20, -1, 1, "thrustAxisCosTheta")

    
    def find_the_event(self):
        dict_event = defaultdict(int)

        file_b2bii = R.TFile("../steeringFile/b2bii_test.root")
        tree_b2bii = file_b2bii.Get("event")

        file_basf = R.TFile("../steeringFile/belle1_steeringFile/exp71_rs2249_re2348_evtgen-uds_0_tree_02.root")
        tree_basf = file_basf.Get("event")

        def check_event(tree, exp_name, run_name, evt_name, nGood_name):
            for entry in range(tree.GetEntries()):
                tree.GetEntry(entry)
                (exp, run, evt)  = (getattr(tree, exp_name), getattr(tree, run_name), getattr(tree, evt_name))
                nGood = getattr(tree, nGood_name)
                if dict_event[(exp, run, evt)] != 0 and dict_event[(exp, run, evt)] != nGood:
                    print(f"Event with different nGood found: exp {exp}, run {run}, evt {evt}, nGood in tree1 {dict_event[(exp, run, evt)]}, in tree2 {nGood}")
                    print(entry)
                dict_event[(exp, run, evt)] = nGood

        #check_event(tree_basf, "expNo", "runNo", "evtNo", "nGood" )
        #check_event(tree_b2bii, "__experiment__", "__run__", "__event__", "NosTrack")
        # entry 71759

        def print_var(tree):
            tree.GetEntry(71759)
            for branch in tree.GetListOfBranches():
                name = branch.GetName()
                value = getattr(tree, name)
                if hasattr(value, '__len__') and not isinstance(value, str):
                    print(f"{name}: {[v for v in value]}")
                else:
                    print(f"{name}: {value}")
                    
        print_var(tree_b2bii)
        print_var(tree_basf)


if __name__ == "__main__":
    ana = Anatest_SpinAlignment()
    ana.bin_by_bin_comparing()
    ana.event_by_event_comparing()
    #ana.test()
    #ana.test2()
    #ana.find_the_event()
    #ana.test3()
