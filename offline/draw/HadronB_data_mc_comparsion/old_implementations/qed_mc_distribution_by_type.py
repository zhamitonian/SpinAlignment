import ROOT as R
from DRAW import style_draw, HistStyle
import os
import sys
import re
from math import pi
sys.path.append("/gpfs/home/belle2/wangz/Work/SpinAlignment/offline")
from Process import SA

def include_QED_mc():
    """
    Include QED MC histograms into the combined plots.
    """
    QED_mc_dir = "/group/belle2/users2022/wangz/other_rootFiles/SpinAlignment/hadronB_data_mc_comparing_QED/rootFiles"
    all_paths = []

    def get_df(type):
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
                 
        df = R.RDataFrame("event", paths)
        df = SA(df)
        return df, paths
    
    types = ['bhabha', 'mumu', 'tautau', 'eeee', 'eemm', 'eecc', 'eeuu', 'eess']
    dfs = []
    for type in types:
        df, paths = get_df(type) 
        dfs.append(df)
        all_paths.extend(paths)

    return dfs, paths

def qed_mc_stack(output_dir):
    types = ['bhabha', 'mumu', 'tautau', 'eeee', 'eemm', 'eecc', 'eeuu', 'eess']
    dfs, _ = include_QED_mc()

    def plot_qed_mc_stack(label, nbin, xmin, xmax, var):
        hists = []

        def get_hist(df, var):
            if re.findall(r'\[.*?\]', var):
                return df.Define("temp", var).Histo1D(("", label, nbin, xmin, xmax), "temp")
            return df.Histo1D(("", label, nbin, xmin, xmax), var)

        for df in dfs:
            hist = get_hist(df, var)
            hists.append(hist)

        png_name = f"{var.replace('[','_').replace(']','_')}"
        png_file = os.path.join(output_dir, f"{png_name}.png")
        style_draw(hists, png_file, styles = [HistStyle.filled_hist(i + 1) for i in range(len(hists))], leg_texts= types)

    plot_qed_mc_stack(";N_{track};[]", 16, 2, 18, "nGood")
    plot_qed_mc_stack(";N_{ECL};[]", 25, 0, 25,  "nECL") 
    plot_qed_mc_stack(";N_{photon};[]", 25, 0, 25,  "nPhoton")
    plot_qed_mc_stack(";E^{*}_{vis};[MeV]", 24, 1, 13, "Evis_cms")
    plot_qed_mc_stack(";E^{*}_{sum};[MeV]", 24, 1, 13, "Esum_cms")
    plot_qed_mc_stack(";#sum #vec P_{z}^{*};[MeV]", 24, -6, 6, "PzSum_cms")
    plot_qed_mc_stack(";E_{vis};[MeV]", 24, 1, 13, "Evis")
    plot_qed_mc_stack(";E_{sum};[MeV]", 24, 1, 9, "Esum")
    plot_qed_mc_stack(";#sum #vec P_{z};[MeV]", 24, -6, 6, "PzSum")

    plot_qed_mc_stack(";M_{heavy jet};[MeV]", 25, 0, 8, "HeavyJetMass")
    plot_qed_mc_stack(";E_{heavy jet};[MeV]", 25, 0, 8, "HeavyJetEnergy")
    plot_qed_mc_stack(";sphericity;[]", 25, 0, 1, "sphericity" )
    plot_qed_mc_stack(";aplanarity;[]", 25, 0, 0.5, "aplanarity")
    plot_qed_mc_stack(";foxWolframeR2;[]", 25, 0, 1, "foxWolfram[2]")
    plot_qed_mc_stack(";cos#theta_{thrust};[]", 20, -1, 1, "thrust[1]")
    plot_qed_mc_stack(";|#vec{T}|;[]", 51, 0.5, 1.01, "thrust[0]")
    plot_qed_mc_stack(";#phi_{thrust};[rad]", 50, -pi, pi, "thrust[2]")

    # p theta phi for trk, cls, pho
    ## lab       
    plot_qed_mc_stack(";p_{track};[MeV]", 50, 0, 5, "trk_p")
    plot_qed_mc_stack(";#theta_{track};[rad]", 50, 0, pi, "trk_theta")
    plot_qed_mc_stack(";#phi_{track};[rad]", 50, -pi, pi, "trk_phi")
    plot_qed_mc_stack(";p_{cluster};[MeV]", 50, 0, 5, "cls_p")
    plot_qed_mc_stack(";#theta_{cluster};[rad]", 50, 0, pi, "cls_theta")
    plot_qed_mc_stack(";#phi_{cluster};[rad]", 50, -pi, pi, "cls_phi")
    plot_qed_mc_stack(";p_{photon};[MeV]", 50, 0, 5, "pho_p")
    plot_qed_mc_stack(";#theta_{photon};[rad]", 50, 0, pi, "pho_theta")
    plot_qed_mc_stack(";#phi_{photon};[rad]", 50, -pi, pi, "pho_phi")

    plot_qed_mc_stack(";dr;[cm]", 20, 0, 1, "trk_dr") 
    plot_qed_mc_stack(";|dz|;[cm]", 40, 0, 4, "trk_dz")

    plot_qed_mc_stack(";PrimeVz;[cm]", 70, -3.5, 3.5, "PrimeVz")
    plot_qed_mc_stack(";PrimeVr;[cm]", 30, 0, 1.5, "PrimeVr")

    
if __name__ == "__main__":
    output_dir = "./test/"
    os.makedirs(output_dir, exist_ok=True)
    #draw_comparison(output_dir)

    output_dir = "qed_mc_stacked_plots"
    os.makedirs(output_dir, exist_ok=True)
    qed_mc_stack(output_dir)






