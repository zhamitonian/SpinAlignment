#!/usr/bin/env python3

"""
with HadronB_dataMC_comp.cc 's output root files
"""

import ROOT as R
from DRAW import style_draw, HistStyle
from math import pi
import os
import re
import sys
from collections import defaultdict

class comparing:
    """
    comparing data with qqbar qed mc 
    """
    def __init__(self, output_dir, data_dir, qqbar_dir, qed_dir):
        output_dir = os.path.join(output_dir, "comparing")
        os.makedirs(output_dir, exist_ok=True)
        self.output_dir = output_dir

        self.data_dir = data_dir
        self.qqbar_dir = qqbar_dir
        self.qed_dir = qed_dir
        pass


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


    def _get_data_hist(self, var):
        data_hists = []

        hist_roots_dir = self.data_dir
        for rootFile_name in os.listdir(hist_roots_dir):
            rootFile_path = os.path.join(hist_roots_dir, rootFile_name)
            rootFile = R.TFile.Open(rootFile_path, "READ")

            if "e73" in rootFile_name: # no coordinate e73 MC
                continue
            
            # Check if file is zombie or invalid
            if not rootFile or rootFile.IsZombie():
                print(f"Warning: Skipping zombie or invalid file: {rootFile_name}")
                if rootFile:
                    rootFile.Close()
                continue

            h_data = rootFile.Get(var)
            if not h_data:
                print(f"Warning: Histogram '{var}' not found in {rootFile_name}")
                rootFile.Close()
                continue
                
            h_data = h_data.Clone()
            h_data.SetDirectory(0)  
            data_hists.append(h_data)
            
            rootFile.Close()
        
        combined_hist = self._combine_hist_list(data_hists)
        return combined_hist
            

    def _get_qqbar_hist(self, var):
        uds_hists = []
        charm_hists = []

        hist_roots_dir = self.qqbar_dir
        for rootFile_name in os.listdir(hist_roots_dir):
            rootFile_path = os.path.join(hist_roots_dir, rootFile_name)
            rootFile = R.TFile.Open(rootFile_path, "READ")
            
            # Check if file is zombie or invalid
            if not rootFile or rootFile.IsZombie():
                print(f"Warning: Skipping zombie or invalid file: {rootFile_name}")
                if rootFile:
                    rootFile.Close()
                continue

            h_qqbar = rootFile.Get(var)
            if not h_qqbar:
                print(f"Warning: Histogram '{var}' not found in {rootFile_name}")
                rootFile.Close()
                continue
                
            h_qqbar = h_qqbar.Clone()
            h_qqbar.SetDirectory(0)

            if "evtgen-uds" in rootFile_name:
                uds_hists.append(h_qqbar)
            elif "evtgen-charm" in rootFile_name:
                charm_hists.append(h_qqbar)
            
            rootFile.Close()    
        
        combined_uds_hist = self._combine_hist_list(uds_hists)
        combined_charm_hist = self._combine_hist_list(charm_hists)

        return combined_uds_hist, combined_charm_hist

    def _get_qed_hist1(self, var): 
        tautau_hists = []
        other_qed_hists = []

        hist_roots_dir = self.qed_dir
        for rootFile_name in os.listdir(hist_roots_dir):
            if rootFile_name.endswith("_0.root"):

                rootFile_path = os.path.join(hist_roots_dir, rootFile_name)
                rootFile = R.TFile.Open(rootFile_path, "READ")
            
                # Check if file is zombie or invalid
                if not rootFile or rootFile.IsZombie():
                    print(f"Warning: Skipping zombie or invalid file: {rootFile_name}")
                    if rootFile:
                        rootFile.Close()
                    continue

                h_qed = rootFile.Get(var)
                if not h_qed:
                    print(f"Warning: Histogram '{var}' not found in {rootFile_name}")
                    rootFile.Close()
                    continue
                
                h_qed = h_qed.Clone()
                h_qed.SetDirectory(0)

                if "tautau" in rootFile_name:
                    tautau_hists.append(h_qed)
                else:
                    other_qed_hists.append(h_qed)
            
                rootFile.Close()

        combined_tautau_hist = self._combine_hist_list(tautau_hists)
        combined_other_qed_hist = self._combine_hist_list(other_qed_hists)

        return combined_tautau_hist, combined_other_qed_hist

    def _get_qed_hist(self, var): 
        tautau_hists = []
        other_qed_hists = []
        hist_groups = defaultdict(list)
        hist_roots_dir = self.qed_dir
        for rootFile_name in os.listdir(hist_roots_dir):
            match = re.match(r'^([a-z]+_e\d+)_\d+\.root$', rootFile_name)
            if match:
                prefix = match.group(1)

                rootFile_path = os.path.join(hist_roots_dir, rootFile_name)
                rootFile = R.TFile.Open(rootFile_path, "READ")

                # Check if file is zombie or invalid
                if not rootFile or rootFile.IsZombie():
                    print(f"Warning: Skipping zombie or invalid file: {rootFile_name}")
                    if rootFile:
                        rootFile.Close()
                    continue

                h_qed = rootFile.Get(var)
                if not h_qed:
                    print(f"Warning: Histogram '{var}' not found in {rootFile_name}")
                    rootFile.Close()
                    continue
                
                h_qed = h_qed.Clone()
                h_qed.SetDirectory(0)

                hist_groups[prefix].append(h_qed)
            
                rootFile.Close()
        
        for prefix, hist_list in hist_groups.items():
            hist = self._combine_hist_list(hist_list)
            hist.Scale(1.0/len(hist_list))
            print(f"Combined {len(hist_list)} files for prefix {prefix}")
            if "tautau" in prefix:
                tautau_hists.append(hist)
            else:
                other_qed_hists.append(hist)

        combined_tautau_hist = self._combine_hist_list(tautau_hists)
        combined_other_qed_hist = self._combine_hist_list(other_qed_hists)

        return combined_tautau_hist, combined_other_qed_hist

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

        x_label_dict = {
            "h_photon_costheta": "cos#theta of cluster",
            "h_photon_p": "energy of cluster (GeV)",
            "h_trk_costheta": "cos#theta of track",
            "h_trk_p": "momentum of track (GeV/c)",
        }
        

        style_draw([h_other_qed, h_tautau, h_ccbar, h_uds, h_data], "", 
                   ["other qed mc", "#tau^{+}#tau^{-}", 'c#bar{c}' , 'uds','data'],
                   styles = [HistStyle.filled_hist(6, 3006), 
                             HistStyle.filled_hist(4, 3007),
                             HistStyle.filled_hist(3, 3011),
                             HistStyle.filled_hist(2, 0),
                             HistStyle.error_bars(1)],
                   pad = pad1, save=False, legend_position= leg_position)
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

        h_ratio.GetXaxis().CenterTitle()
        h_ratio.GetXaxis().SetTitle(x_label_dict.get(var, h_ratio.GetXaxis().GetTitle()))

        png_name = f"{var.replace('[','_').replace(']','_')}"
        png_file = os.path.join(output_dir, f"{png_name}.png")
        style_draw([h_ratio], png_file, styles = [HistStyle.error_bars(1)], y_min= 0, y_max=2, use_user_y_range= True,pad= pad2 )


    def start_plotting(self):
        var_list = [
            "h_Aplanarity",
            "h_Esum",
            "h_Evis",
            "h_HeavyJetEnergy",
            "h_HeavyJetMass",
            "h_HeavyJetMass_over_Evis",
            "h_Ncls",
            "h_Ntrk",
            "h_photon_costheta",
            "h_photon_p",
            "h_photon_phi",
            "h_PrimeR",
            "h_PrimeZ",
            "h_Psum",
            "h_Pz",
            "h_R2",
            "h_Sphericity",
            "h_thrust",
            "h_thrust_costheta",
            "h_thrust_phi",
            "h_trk_costheta",
            "h_trk_dr",
            "h_trk_dz",
            "h_trk_p",
            "h_trk_phi",
        ]


        #for var in var_list:
        #for var in ["h_trk_p", "h_photon_p", "h_trk_costheta", "h_photon_costheta"]:
        for var in ["h_photon_costheta"]:
            data_hist = self._get_data_hist(var)
            uds_hist, ccbar_hist = self._get_qqbar_hist(var)
            uds_hist.Scale(1.0/4)
            ccbar_hist.Scale(1.0/4)
            tautau_hist, other_qed_hist = self._get_qed_hist(var)

            hist_list = [data_hist, uds_hist, ccbar_hist, tautau_hist, other_qed_hist]

            leg_position = 2
            if var in ['h_HeavyJetEnergy', 'h_thrust_costheta', 'h_thrust', 'h_trk_costheta']:
                leg_position = 0
            
            self.brush(hist_list, self.output_dir, var, leg_position=leg_position, normalize= False)
            normalized_output_dir = self.output_dir + "_normalized"
            os.makedirs(normalized_output_dir, exist_ok=True)
            self.brush(hist_list, normalized_output_dir, var, leg_position=leg_position, normalize= True)
            

    def test_qedMC_read(self):
        for var in ["h_Esum", "h_Evis"]:
            tautau_hist, other_qed_hist = self._get_qed_hist(var)
            tautau_hist1, other_qed_hist1 = self._get_qed_hist1(var)

            style_draw([other_qed_hist ,other_qed_hist1], f"{self.output_dir}/comparison_other_qed_{var}.png",
                       ["other qed mc: method 1", "other qed mc: method 2"],
                       [HistStyle.line_hist(4), 
                        HistStyle.error_bars(1)]) 
            style_draw([tautau_hist ,tautau_hist1], f"{self.output_dir}/comparison_tautau_{var}.png",
                       ["#tau^{+}#tau^{-} mc: method 1", "#tau^{+}#tau^{-} mc: method 2"],
                       [HistStyle.line_hist(4), 
                        HistStyle.error_bars(1)]) 
    
if __name__ == "__main__":
    dataMC_dir = "/home/belle2/wangz/disk2/data_gMC_belle1/HadronB_DataMC_comp_v1.0.2"

    comparing_instance = comparing(
        output_dir = "./images_v1.0.2",
        data_dir = f"{dataMC_dir}_data/4S_offres",
        qqbar_dir= f"{dataMC_dir}_qqbar/4S_offres",
        qed_dir =  f"{dataMC_dir}_qed/rootFiles")
    comparing_instance.start_plotting()

