#!/usr/bin/env python3

import ROOT as R
from DRAW import style_draw, HistStyle
from math import pi
import os
import re
import sys

class concise_comparing:
    """
    comparing data with qqbar qed mc 
    """
    def __init__(self, output_dir, data_dir, qqbar_dir, qed_dir):
        output_dir = os.path.join(output_dir, "concise_comparing")
        os.makedirs(output_dir, exist_ok=True)
        self.output_dir = output_dir

        self.data_dir = data_dir
        self.qqbar_dir = qqbar_dir
        self.qed_dir = qed_dir
        pass

    def _get_data_df(self):
        path_list = []
        
        for file in os.listdir(self.data_dir):
            if file.endswith(".root"):
                full_path = os.path.join(self.data_dir, file)
                path_list.append(full_path)
        df = R.RDataFrame("event", path_list)

        print(f"Number of files processed: {len(path_list)}")
        return df

    def _get_qqbar_df(self):
        path_list_uds = []
        path_list_ccbar = []

        def get_df_paths(type):
            paths = []
            for file in os.listdir(self.qqbar_dir):
                if file.endswith(".root") and type in file:
                    full_path = os.path.join(self.qqbar_dir, file)
                    
                    try:
                        test_file = R.TFile.Open(full_path, "READ")
                        if not test_file.IsZombie():
                            tree = test_file.Get("event")
                            if tree and tree.GetEntries() > 0:
                                paths.append(full_path)
                            else:
                                print(f"Corrupted file skipped: {full_path}")
                        else :
                            print("Skipping invalid file:", full_path)
                        if test_file:
                            test_file.Close()
                    except Exception as e:
                        print(f"Error opening file {full_path}: {e}")
                        continue
            return paths

        types = ['evtgen-charm', 'evtgen-uds']
        path_list_ccbar = get_df_paths(types[0])
        path_list_uds = get_df_paths(types[1])

        df_ccbar = R.RDataFrame("event", path_list_ccbar)
        df_uds = R.RDataFrame("event", path_list_uds)
        
        print(f"Number of files processed: {len(path_list_ccbar) + len(path_list_uds)}")
        return df_ccbar, df_uds

    def _get_qed_df(self):
        dfs = []

        def get_df_paths(type):
            paths = []
            for file in os.listdir(self.qed_dir):
                if file.endswith("_0.root") and file.startswith(type):
                    full_path = os.path.join(self.qed_dir, file)

                    try:
                        test_file = R.TFile.Open(full_path, "READ")
                        if not test_file.IsZombie():
                            tree = test_file.Get("event")
                            if tree and tree.GetEntries() > 0:
                                paths.append(full_path)
                            else:
                                print(f"Corrupted file skipped: {full_path}")
                        else :
                            print("Skipping invalid file:", full_path)
                        if test_file:
                            test_file.Close()
                    except Exception as e:
                        print(f"Error opening file {full_path}: {e}")
                        continue
            return paths

        types = ['tautau'] 
        other_types = ['mumu' ,'bhabha' ,'eemm', 'eeee', 'eecc', 'eeuu', 'eess']

        for type in types:
            paths = get_df_paths(type)
            if paths:  # Only create RDataFrame if paths is not empty
                df = R.RDataFrame("event", paths)
                dfs.append(df)
            else:
                print(f"Warning: No valid files found for type {type}")

        other_qed_paths = []
        for type in other_types:
            paths = get_df_paths(type)
            other_qed_paths.extend(paths)
        
        if other_qed_paths:  # Only create RDataFrame if diPhoton_paths is not empty
            df_diPhoton = R.RDataFrame("event", other_qed_paths)
            dfs.append(df_diPhoton)
        else:
            print("Warning: No valid files found for diPhoton types")

        return dfs

    def _get_hist(self, df, var):

        hist_config = {
            "Ntrk": (";N_{track};[]", 16, 2, 18),
            "Ncls": (";N_{cls};[]", 25, 0, 25),
            "Evis": (";E^{*}_{vis};[MeV]", 14, 6, 13),
            "Esum": (";E^{*}_{sum};[MeV]", 24, 1, 13),
            "Pz": (";#sum #vec P_{z}^{*};[MeV]", 24, -6, 6),
            "Psum" :(";|#sum #vec P^{*}|;[MeV]", 18, 1, 10),
            
            "HeavyJetMass": (";M_{heavy jet};[MeV]", 25, 0, 8),
            "R2": (";foxWolframeR2;[]", 25, 0, 1),
            "thrust[1]": (";cos#theta_{thrust};[]", 20, -1, 1),
            "thrust[0]": (";|#vec{T}|;[]", 51, 0.5, 1.01),
            "thrust[2]": (";#phi_{thrust};[rad]", 50, -pi, pi),
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

    def _combine_hist_list(self, hist_list):
        """
        Combine all histograms in the list into a single histogram
        Args:
            hist_list: List of ROOT histograms
        Returns:
            Combined histogram
        """
        if not hist_list:
            return None
        
        # Clone first histogram to preserve binning and axis properties
        result = hist_list[0].Clone()
        
        # Add all other histograms
        for hist in hist_list[1:]:
            result.Add(hist)
        
        return result

    def brush(self, hist_list, output_dir, var, leg_position=2, normalize=False):
        h_data = hist_list[0]
        MC_hists = hist_list[1:]
        h_MC = self._combine_hist_list(MC_hists)
        scale_factor = h_data.Integral() / h_MC.Integral()
        if normalize:
            h_MC.Scale(scale_factor)
            for h in MC_hists:
                h.Scale(scale_factor)

        h_uds = MC_hists[0]
        h_ccbar = MC_hists[1]
        # wrong ?  udsc mc didn't scale ?
        h_tautau = MC_hists[2]
        h_other_qed = MC_hists[3]
        
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
        style_draw([h_uds, h_ccbar, h_tautau, h_other_qed, h_data], "", 
                   ["uds", "c#bar{c}",  "#tau^{+}#tau^{-}", "other qed mc", 'data'],
                   styles = [HistStyle.filled_hist(i+2) for i in range(4)] + [HistStyle.error_bars(1)], pad = pad1, save=False, legend_position= leg_position)
        h_ratio = h_MC.Clone("h_ratio")
        h_ratio.Sumw2()
        h_ratio.Divide(h_MC, h_data, 1, 1, "B")
        h_ratio.GetYaxis().SetTitle("ratio")

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
        png_file = os.path.join(output_dir, f"{png_name}.png")
        style_draw([h_ratio], png_file, styles = [HistStyle.error_bars(1)], y_min= 0, y_max=2, use_user_y_range= True,pad= pad2 )


    def start_plotting(self):
        df_list = []
        data_df = self._get_data_df()
        df_list.append(data_df)
        qqbar_ccbar_df, qqbar_uds_df = self._get_qqbar_df()
        df_list.append(qqbar_uds_df)
        df_list.append(qqbar_ccbar_df)
        qed_dfs = self._get_qed_df()
        df_list.extend(qed_dfs)

        var_list = [
            "Ntrk", "Ncls", "Evis", "Esum", "Pz", "Psum",
            "HeavyJetMass", "R2", "thrust[1]", "thrust[0]", "thrust[2]",
        ]

        for var in var_list:
            hist_list = []
            for df in df_list:
                hist = self._get_hist(df, var)
                hist_list.append(hist)
            #style_draw(hist_list, f"./temp/{var}.png")
            
            self.brush(hist_list, self.output_dir, var, leg_position=2, normalize= False)
            normalized_output_dir = self.output_dir + "_normalized"
            os.makedirs(normalized_output_dir, exist_ok=True)
            self.brush(hist_list, normalized_output_dir, var, leg_position=2, normalize= True)
            
    
if __name__ == "__main__":
    concise_comparing_instance = concise_comparing(
        output_dir = "./images",
        data_dir = "/home/belle2/wangz/disk2/data_gMC_belle1/2025-12-01_SpinAlignment_data/continuum",
        #qqbar_dir= "/home/belle2/wangz/disk2/data_gMC_belle1/2025-12-01_SpinAlignment_gMC/continuum",
        qqbar_dir= "/gpfs/group/belle/users/wangz/data_gMC_belle1/2025-12-02_SpinAlignment_qqbarMC/continuum",
        qed_dir=   "/home/belle2/wangz/disk2/data_gMC_belle1/2025-12-02_SpinAlignment_qedMC/rootFiles")
    concise_comparing_instance.start_plotting()

