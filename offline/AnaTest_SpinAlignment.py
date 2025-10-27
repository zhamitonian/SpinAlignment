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

    def check(self):
        file_b2bii = R.TFile("../steeringFile/b2bii_test.root", 'READ')
        tree= file_b2bii.Get("event")
        tree.GetEntry(0)
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




if __name__ == "__main__":
    ana = Anatest_SpinAlignment()
   # ana.check()
    #ana.get_the_entry_index()
    ana.little_test()
