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

class Anatest_SpinAlignment:
    def __init__(self):
        self.tools = RDF_process()

    def comparing_basf_b2bii(self):
        #df_b2bii = R.RDataFrame("event", "../steeringFile/b2bii_test02.root")
        df_b2bii = R.RDataFrame("event", "../steeringFile/b2bii_test.root")
        df_basf =  R.RDataFrame("event", "../steeringFile/belle1_steeringFile/exp71_rs2249_re2348_evtgen-uds_0_tree_02.root")
        print(df_b2bii.Count().GetValue(), df_basf.Count().GetValue())
        #df_basf = R.RDataFrame("event", "../steeringFile/belle1_steeringFile/basf_test.root")

        #df_basf = df_basf.Define("thrustAxisCosTheta", "cos(thrust[1])")
        
        def plot(label, nbin, xmin, xmax, var, var2=None):
            def get_hist(df, var):
                if re.findall(r'\[.*?\]', var):
                    return df.Define("temp", var).Histo1D(("", label, nbin, xmin, xmax), "temp")
                return df.Histo1D(("", label, nbin, xmin, xmax), var)

            h_b2bii = get_hist(df_b2bii, var).GetValue()
            h_basf = get_hist(df_basf, var2 if var2 else var).GetValue()

            if var2 is not None:
                var = var2.replace("[","_").replace("]","_") 
            style_draw([h_b2bii, h_basf], "./draw/b2bii_basf_comparison/{}.png".format(var), ["b2bii","basf"],styles = [HistStyle.line_hist(4,1,2), HistStyle.error_bars(1)])

            h_b2bii.Add(h_basf, -1)  # Subtract h1 from h0
            h_b2bii.SetTitle(label.split("/")[-1])
            style_draw([h_b2bii], "./draw/b2bii_basf_comparison/{}_diff.png".format(var), ["b2bii - basf"],styles = [HistStyle.error_bars(1)])

        plot(";N_{track};", 16, 2, 18, "NosTrack", "nGood")
        plot(";N_{Cluster};", 25, 0, 25, "NosCluster", "nCluster") 
        plot(";N_{photon};", 25, 0, 25, "NosPhoton", "nPhoton")
        plot(";E_{vis};MeV", 24, 1, 13, "EvisCMS", "Evis_cms")
        plot(";#sum #vec P_{z};MeV", 24, -6, 6, "BalancePzCMS", "BalancePz_cms")
        plot(";Eenergy_{cms};MeV", 24, 1, 13, "EnergyCMS", "Energy_cms")
        plot(";|#vec{T}|;", 24, 0.5, 1.1, "thrust", "thrust[0]")
        plot(";E_{cms};MeV", 36, 2, 11, "Ecms")
        plot(";M_{heavy jet};MeV", 25, 0, 8, "HeavyJetMass")
        plot(";M_{heavy jet};MeV", 25, 0, 8, "HeavyJetEnergy")
        plot(";sphericity;", 25, 0, 1, "sphericity" )
        plot(";aplanarity;", 25, 0, 0.5, "aplanarity")
        plot(";foxWolframeR2;", 25, 0, 1, "foxWolframR2", "foxWolfram[2]")
        plot(";cos#theta_{thrust};", 20, -1, 1, "thrustAxisCosTheta", "thrust[1]")

        '''
        plot(";p_{track}; MeV", 25, 0, 6, "", "pip_p")
        plot(";#theta_{track};rad", 20, 0, pi, "", "pip_theta")
        plot(";#phi_{track};rad", 20, -pi, pi, "", "pip_phi")
        plot(";p_{#gamma}",21, 0, 7, "photon_p")
        plot(";#theta_{#gamma};rad", 20, 0, pi, "photon_theta")
        plot(";#phi_{#gamma};rad", 20, -pi, pi, "photon_phi")
        '''

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



if __name__ == "__main__":
    ana = Anatest_SpinAlignment()
    ana.comparing_basf_b2bii()
    #ana.test()
    #ana.test2()
