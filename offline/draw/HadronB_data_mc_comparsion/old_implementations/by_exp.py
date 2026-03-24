#!/usr/bin/env python3

import ROOT as R
from DRAW import style_draw, HistStyle
from math import pi
import os
import re
import sys
import pandas as pd
from collections import defaultdict

sys.path.append("/gpfs/home/belle2/wangz/Work/SpinAlignment/offline")
from Process import SA


class comapring_by_exp:
    """
    comparing data with qqbar mc for a given expNo , and save the histograms into root files
    """
    def __init__(self, output_dir, data_dir, qqbar_dir, expNo):
        exp_dir = os.path.join(output_dir, "images", f"exp{expNo}")
        os.makedirs(exp_dir, exist_ok=True)
        self.exp_dir = exp_dir

        root_dir = os.path.join(output_dir, "hist_roots")
        os.makedirs(root_dir, exist_ok=True)
        self.rootFilePath = os.path.join(root_dir, f"exp{expNo}.root")

        self.data_dir = data_dir
        self.qqbar_dir = qqbar_dir
        self.expNo = expNo
        self.correct_dz = False

    def read_dz_shift(self):
        """
        read the dz shift for different exp from file
        """
        if os.path.exists('images/primeVz_fitting/center_primeVz_byExp.txt'):
            pd_rf = pd.read_csv('images/primeVz_fitting/center_primeVz_byExp.txt', 
                        sep='\s+', 
                        comment='#',
                        names=['exp_id', 'mean', 'sigma'])
            mean_dict = pd_rf.set_index('exp_id')['mean'].to_dict()
            self.correct_dz = True
        else: 
            mean_dict = defaultdict(lambda: 0)
        
        self.z_shift_dict = mean_dict

    def _get_data_df(self):
        path_list = []
        for file in os.listdir(self.data_dir):
            if file.endswith(".root") and f"_e{self.expNo}_" in file:
                full_path = os.path.join(self.data_dir, file)
                path_list.append(full_path)
        df = R.RDataFrame("event", path_list)
        df = SA(df)
        if self.correct_dz:
            df = df.Redefine("PrimeVz", f"PrimeVz - {self.z_shift_dict[self.expNo]}")

        # tempolarily fix
        df = df.Define("trk_cosTheta", "cos(trk_theta)") 
        df = df.Define("pho_cosTheta", "cos(pho_theta)")
        
        
        print(f"Number of files processed: {len(path_list)}")
        return df

    def _get_qqbar_df(self):
        path_list_uds = []
        path_list_ccbar = []
        for file in os.listdir(self.qqbar_dir):
            if file.endswith(".root") and f"_e{self.expNo}_" in file:
                full_path = os.path.join(self.qqbar_dir, file)
                if "evtgen-charm" in file:
                    path_list_ccbar.append(full_path)
                elif "evtgen-uds" in file:
                    path_list_uds.append(full_path)
        df_ccbar = R.RDataFrame("event", path_list_ccbar)
        df_uds = R.RDataFrame("event", path_list_uds)
        df_ccbar = SA(df_ccbar)
        df_uds = SA(df_uds)
        if self.correct_dz:
            df_ccbar = df_ccbar.Redefine("PrimeVz", f"PrimeVz - {self.z_shift_dict[self.expNo]}").Filter("PrimeVz > -1.5 && PrimeVz < 1.5")
            df_uds = df_uds.Redefine("PrimeVz", f"PrimeVz - {self.z_shift_dict[self.expNo]}").Filter("PrimeVz > -1.5 && PrimeVz < 1.5")
        
        # tempolarily fix
        df_ccbar = df_ccbar.Define("trk_cosTheta", "cos(trk_theta)") 
        df_ccbar = df_ccbar.Define("pho_cosTheta", "cos(pho_theta)")
        df_uds = df_uds.Define("trk_cosTheta", "cos(trk_theta)") 
        df_uds = df_uds.Define("pho_cosTheta", "cos(pho_theta)")
        print(f"Number of files processed: {len(path_list_ccbar) + len(path_list_uds)}")
        return df_ccbar, df_uds

    def _get_hist(self, df, var):

        hist_config = {
            "nGood": (";N_{track};[]", 16, 2, 18),
            "nECL": (";N_{ECL};[]", 25, 0, 25),
            "nPhoton": (";N_{photon};[]", 25, 0, 25),
            "Evis_cms": (";E^{*}_{vis};[MeV]", 24, 1, 13),
            "Esum_cms": (";E^{*}_{sum};[MeV]", 24, 1, 13),
            "PzSum_cms": (";#sum #vec P_{z}^{*};[MeV]", 24, -6, 6),
            "Evis": (";E_{vis};[MeV]", 24, 1, 13),
            "Esum": (";E_{sum};[MeV]", 24, 1, 9),
            "PzSum": (";#sum #vec P_{z};[MeV]", 24, -6, 6),
            
            "HeavyJetMass": (";M_{heavy jet};[MeV]", 25, 0, 8),
            "HeavyJetEnergy": (";E_{heavy jet};[MeV]", 25, 0, 8),
            "sphericity": (";sphericity;[]", 25, 0, 1),
            "aplanarity": (";aplanarity;[]", 25, 0, 0.5),
            "foxWolfram[2]": (";foxWolframeR2;[]", 25, 0, 1),
            "thrust[1]": (";cos#theta_{thrust};[]", 20, -1, 1),
            "thrust[0]": (";|#vec{T}|;[]", 51, 0.5, 1.01),
            "thrust[2]": (";#phi_{thrust};[rad]", 50, -pi, pi),
            
            # Track, cluster, photon
            "trk_p": (";p_{track};[MeV]", 50, 0, 5),
            "trk_theta": (";#theta_{track};[rad]", 50, 0, pi),
            "trk_phi": (";#phi_{track};[rad]", 50, -pi, pi),
            "cls_p": (";p_{cluster};[MeV]", 50, 0, 5),
            "cls_theta": (";#theta_{cluster};[rad]", 50, 0, pi),
            "cls_phi": (";#phi_{cluster};[rad]", 50, -pi, pi),
            "pho_p": (";p_{photon};[MeV]", 50, 0, 5),
            "pho_theta": (";#theta_{photon};[rad]", 50, 0, pi),
            "pho_phi": (";#phi_{photon};[rad]", 50, -pi, pi),
            
            "trk_dr": (";dr;[cm]", 20, 0, 1),
            "trk_dz": (";|dz|;[cm]", 40, 0, 4),
            
            "PrimeVz": (";PrimeVz;[cm]", 70, -3.5, 3.5),
            "PrimeVr": (";PrimeVr;[cm]", 30, 0, 1.5),

            "trk_cosTheta": (";cos#theta_{track};[]", 40, -1, 1),
            "pho_cosTheta": (";cos#theta_{photon};[]", 40, -1, 1),
        }
        
        if var not in hist_config:
            print(f"Warning: Variable {var} not found in hist_config")
            return None
        
        label, nbin, xmin, xmax = hist_config[var]
        
        # Create histogram from RDataFrame
        def get_hist(df, var):
            if re.findall(r'\[.*?\]', var):
                return df.Define("temp", var).Histo1D(("", label, nbin, xmin, xmax), "temp").GetValue()
            return df.Histo1D(("", label, nbin, xmin, xmax), var).GetValue() 
        
        hist = get_hist(df, var)

        return hist

    def brush(self, h_data, h_uds, h_ccbar, var):
        h_mc = h_uds.Clone("h_mc")
        h_mc.Add(h_ccbar)

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

        style_draw([h_uds, h_ccbar, h_data], "", ["uds", "c#bar{c}", "data"],
                   styles = [HistStyle.filled_hist(6), HistStyle.filled_hist(7), HistStyle.error_bars(1)], 
                   pad = pad1, save=False)

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
        png_file = os.path.join(self.exp_dir, f"{png_name}.png")
        style_draw([h_ratio], png_file, styles = [HistStyle.error_bars(1)], y_min= 0, y_max=2, use_user_y_range= True,pad= pad2)


    def start_plotting(self):
        self.read_dz_shift()
        df_data = self._get_data_df()
        if self.correct_dz:
            df_data = df_data.Filter("PrimeVz > -1.5 && PrimeVz < 1.5")
            df_sideband = df_data.Filter(" (PrimeVz >=2.0 && PrimeVz <=3) || (PrimeVz <= -2.0 && PrimeVz >= -3) ")
        
        df_ccbar, df_uds = self._get_qqbar_df()

        for var in ['nGood', 'nECL', 'nPhoton', 'Evis_cms', 'Esum_cms', 'PzSum_cms', 'Evis', 'Esum', 'PzSum',
                'HeavyJetMass', 'HeavyJetEnergy', 'sphericity', 'aplanarity', 'foxWolfram[2]',
                'thrust[0]', 'thrust[1]', 'thrust[2]', 'trk_p', 'trk_theta', 'trk_phi',
                'cls_p', 'cls_theta', 'cls_phi', 'pho_p', 'pho_theta', 'pho_phi', 'trk_dr', 
                'trk_dz', 'PrimeVr', 'PrimeVz', 'trk_cosTheta', 'pho_cosTheta']:
            h_data = self._get_hist(df_data, var)
            h_uds = self._get_hist(df_uds, var)
            h_ccbar = self._get_hist(df_ccbar, var)
            if self.correct_dz:
                h_sideband = self._get_hist(df_sideband, var)

            self.brush( h_data, h_uds, h_ccbar, var )

            rootFile = R.TFile.Open(self.rootFilePath, "UPDATE")
            rootFile.cd()
            h_data.Write(f"{var}_data")
            h_uds.Write(f"{var}_uds")
            h_ccbar.Write(f"{var}_ccbar")
            if self.correct_dz:
                h_sideband.Write(f"{var}_data_sideband")
            rootFile.Close()

def main():
    data_dir = "/gpfs/group/belle2/users2022/wangz/data_gMC_belle1/2025-10-30_SpinAlignmnet_Data/continuum/"
    gMC_dir = "/group/belle2/users2022/luruihua/for_wangz/data_gMC_belle1/2025-10-30_SpinAlignmnet_gMC/continuum/"
    output_dir = sys.argv[1]
    expNo = int(sys.argv[2])

    comparer = comapring_by_exp(output_dir, data_dir, gMC_dir, expNo)
    comparer.start_plotting()


if __name__ == "__main__":
    main()