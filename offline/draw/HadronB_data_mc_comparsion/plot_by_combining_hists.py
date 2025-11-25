import ROOT as R
from DRAW import style_draw, HistStyle
import os
import sys
import re
from math import pi
sys.path.append("/gpfs/home/belle2/wangz/Work/SpinAlignment/offline")
from Process import SA

class combine_plot:
    def __init__(self):
        pass

    def combine_hist_list(self, hist_list):
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

    def _get_data_gMC_hists(self, var):
        data_hists = []
        ccbar_hists = []
        uds_hists = []
        data_sideband_hists = []

        for rootFile_name in os.listdir("./hist_roots"):
            rootFile_path = os.path.join("./hist_roots", rootFile_name)
            rootFile = R.TFile.Open(rootFile_path, "READ")
                    
            h_data = rootFile.Get(f"{var}_data").Clone()
            h_uds = rootFile.Get(f"{var}_uds").Clone()
            h_ccbar = rootFile.Get(f"{var}_ccbar").Clone()
            h_data_sideband = rootFile.Get(f"{var}_data_sideband").Clone()
            h_data.SetDirectory(0)  # Detach from file
            h_uds.SetDirectory(0)   
            h_ccbar.SetDirectory(0)
            h_data_sideband.SetDirectory(0) 

            if h_data and h_uds and h_ccbar:  # Check if histograms exist
                data_hists.append(h_data)  # Clone to keep after file closes
                ccbar_hists.append(h_ccbar)
                uds_hists.append(h_uds)
                data_sideband_hists.append(h_data_sideband)
                    
            rootFile.Close()
                
        # Combine histograms
        combined_data = self.combine_hist_list(data_hists)
        combined_ccbar = self.combine_hist_list(ccbar_hists)
        combined_uds = self.combine_hist_list(uds_hists)
        combined_data_sideband = self.combine_hist_list(data_sideband_hists)
        
        return combined_data, combined_uds, combined_ccbar, combined_data_sideband

    def include_QED_mc(self):
        """
        Include QED MC histograms into the combined plots.
        """
        QED_mc_dir = "/group/belle2/users2022/wangz/other_rootFiles/SpinAlignment/hadronB_data_mc_comparing_QED/rootFiles"

        dfs = []

        def get_df_paths(type):
            paths = []
            for file in os.listdir(QED_mc_dir):
                if file.endswith("_0.root") and file.startswith(type):
                    full_path = os.path.join(QED_mc_dir, file)

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

        types = ['bhabha', 'mumu', 'tautau', 'eeee', 'eemm'] 
        diPhoton_types = ['eecc', 'eeuu', 'eess']

        for type in types:
            paths = get_df_paths(type) 
            df = R.RDataFrame("event", paths)
            df = SA(df)
            df = df.Filter("PrimeVz <1.5 && PrimeVz > -1.5")
            dfs.append(df)

        diPhoton_paths = []
        for type in diPhoton_types:
            paths = get_df_paths(type)
            diPhoton_paths.extend(paths)
        df_diPhoton = R.RDataFrame("event", diPhoton_paths)
        df_diPhoton = SA(df_diPhoton)
        df_diPhoton = df.Filter("PrimeVz <1.5 && PrimeVz > -1.5")
        dfs.append(df_diPhoton)

        self.qed_dfs = dfs

        return paths
    
    def _get_qed_hist(self, var):
        if not hasattr(self, 'qed_df'):
            self.include_QED_mc()

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
        
        hists = []
        for df in self.qed_dfs:
            hist = get_hist(df, var)
            hists.append(hist)

        return hists

    def brush(self, MC_hists, h_data, var , output_dir, leg_position =2, normalize = False):

        h_MC = self.combine_hist_list(MC_hists)
        scale_factor = h_data.Integral() / h_MC.Integral()
        if normalize:
            h_MC.Scale(scale_factor)
            for h in MC_hists:
                h.Scale(scale_factor)

        h_uds = MC_hists[0]
        h_ccbar = MC_hists[1]
        h_bhabha = MC_hists[2]
        h_mumu = MC_hists[3]
        h_tautau = MC_hists[4]
        h_eeee = MC_hists[5]
        h_eemm = MC_hists[6]
        h_diPhoton = MC_hists[7]
        h_beamwall = MC_hists[8]
        
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
        style_draw([h_uds, h_ccbar, h_bhabha, h_mumu,h_tautau, h_eeee, h_eemm, h_diPhoton,h_beamwall, h_data], "", 
                   ["uds", "c#bar{c}", "bhabha", "#mu^{+}#mu^{-}", "#tau^{+}#tau^{-}", "e^{+}e^{-}e^{+}e^{-}", "e^{+}e^{-}#mu^{+}#mu^{-}", "e^{+}e^{-}#gamma#gamma","beam wall", 'data'],
                   styles = [HistStyle.filled_hist(i+2) for i in range(8)] + [HistStyle.filled_hist(11), HistStyle.error_bars(1)], pad = pad1, save=False, legend_position= leg_position)
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


    def start_plotting(self, output_dir, normalize = False):

        var_list = ['nGood', 'nECL', 'nPhoton', 'Evis_cms', 'Esum_cms', 'PzSum_cms', 'Evis', 'Esum', 'PzSum',
                'HeavyJetMass', 'HeavyJetEnergy', 'sphericity', 'aplanarity', 'foxWolfram[2]', 
                'thrust[0]', 'thrust[1]', 'thrust[2]', 'trk_p', 'trk_theta', 'trk_phi',
                'cls_p', 'cls_theta', 'cls_phi', 'pho_p', 'pho_theta', 'pho_phi', 'trk_dr', 'trk_dz', 'PrimeVr', 'PrimeVz']

        for var in var_list:
            h_data, h_uds, h_ccbar, h_data_sideband = self._get_data_gMC_hists(var)
            h_data_sideband.Scale(1.5)
            h_qed_mcs = self._get_qed_hist(var)
            MC_hists = [h_uds, h_ccbar] +  h_qed_mcs +  [h_data_sideband]

            leg_position = 2
            if var in ["Evis_cms" , "HeavyJetEnergy", 'thrust[0]', 'thrust[1]', 'PzSum_cms']:
                leg_position = 0

            self.brush(MC_hists ,h_data, var , output_dir, leg_position= leg_position, normalize= normalize)
            self.brush(MC_hists ,h_data, var ,output_dir= "./images/combined_plots_normalized", leg_position= leg_position, normalize= True)


if __name__ == "__main__":
    combiner = combine_plot()
    output_dir = "./images/combined_plots"
    os.makedirs(output_dir, exist_ok=True)
    combiner.start_plotting(output_dir)
    






