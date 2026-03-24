#!/usr/bin/env python3

import ROOT as R
from DRAW import style_draw, HistStyle
import os

def bkg_dist(MC_rootFile, output_dir = "./"):
    os.makedirs(output_dir, exist_ok=True)

    rdf = R.RDataFrame("event", MC_rootFile)
    condition = "pip_isSignal == false && pim_isSignal == false"
    for var in ["Ks_M", "Ks_z", "Ks_helicity_angle"]:
        rdf = rdf.Redefine(var, f"{var}[{condition}]") 
    
    rdf = rdf.Filter("Ks_M.size() > 0", "has Ks after truth match")
    rdf = rdf.Filter("n_ks_truth == 1", "no truth Ks") 

    rdf.Report().Print()

    # -------  1d distributions
    hist_M = rdf.Histo1D(("", ";M_{K_{s}};[MeV]", 100, 0.47, 0.525), "Ks_M")
    style_draw([hist_M], os.path.join(output_dir, "bkg_M_Ks_dist.png"), styles=[HistStyle.filled_line_hist(4,3011)])

    hist_z = rdf.Histo1D(("", ";z;[1]", 20, 0, 1), "Ks_z")
    style_draw([hist_z], os.path.join(output_dir, "bkg_z_Ks_dist.png"), styles=[HistStyle.filled_line_hist(4,3011)])

    hist_helicity_angle = rdf.Histo1D(("", ";cos#theta*;[1]", 10, -1, 1), "Ks_helicity_angle")
    style_draw([hist_helicity_angle], os.path.join(output_dir, "bkg_helicity_angle_Ks_dist.png"), styles=[HistStyle.filled_line_hist(4,3011)])

    # -------  2d distribution
    hist_model = R.RDF.TH2DModel("Ks_z_vs_helicity_angle_bkg", "Ks_z vs helicity_angle for bkg;z;cos#theta^{#star}", 20, 0, 1, 10, -1, 1)
    bkg_dist = rdf.Histo2D(hist_model, "Ks_z", "Ks_helicity_angle")
        
    R.gStyle.SetOptStat(0)
    R.gStyle.SetPaintTextFormat(".0f")

    c = R.TCanvas("c0", "c0", 1600, 1080)
    bkg_dist.Draw("COLZ TEXT")
    c.SaveAs(os.path.join(output_dir, "Ks_costheta_vs_z_reco_truth_matched.png"))

    rdf.Snapshot("event", os.path.join(output_dir, "processed_bkg.root"))

def bkg_dist_check_p_cut(MC_rootFile, MC_rootFile2, output_dir = "./"):
    os.makedirs(output_dir, exist_ok=True)

    rdf = R.RDataFrame("event", MC_rootFile)
    rdf2 = R.RDataFrame("event", MC_rootFile2)

    condition = "pip_isSignal == false && pim_isSignal == false"
    for var in ["Ks_M", "Ks_z", "Ks_helicity_angle"]:
        rdf = rdf.Redefine(var, f"{var}[{condition}]") 
        rdf2 = rdf2.Redefine(var, f"{var}[{condition}]")
    
    rdf = rdf.Filter("Ks_M.size() > 0", "has Ks after truth match")
    rdf2 = rdf2.Filter("Ks_M.size() > 0", "has Ks after truth match")

    rdf.Report().Print()
    rdf2.Report().Print()

    # -------  1d distributions
    hist_M = rdf.Histo1D(("", ";M_{K_{s}};[MeV]", 110, 0.47, 0.525), "Ks_M")
    hist_M2 = rdf2.Histo1D(("", ";M_{K_{s}};[MeV]", 110, 0.47, 0.525), "Ks_M")
    style_draw([hist_M, hist_M2], os.path.join(output_dir, "momemtum_cut_M_ks_dist.png"), 
               leg_texts= ["no momentum cut", "p_{lab} > 0.5 GeV/c"],
               styles=[HistStyle.filled_line_hist(4,3011)])

if __name__ == "__main__":
    MC_rootFile_v210 = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.1.0_qqbar/exp55_reco_processed.root"
    MC_rootFile_v220 = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.2.0_qqbar/exp55_reco_processed.root"
    MC_rootFile_v210_p = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.1.0_qqbar/exp55_reco_processed_test.root"
    MC_rootFile_v230 = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.0_qqbar/exp55_reco_processed.root"
    MC_rootFile_v231 = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.1_qqbar/exp55_reco_processed.root"

    #bkg_dist_check_p_cut(MC_rootFile, MC_rootFile2)
    #bkg_dist(MC_rootFile_v220, "v2.2.0")
    #bkg_dist(MC_rootFile_v210, "v2.1.0")
    #bkg_dist(MC_rootFile_v210_p, "v2.1.0_with_p_cut")
    
    #bkg_dist(MC_rootFile_v230, "v2.3.0_n_ks_truth_eq_0")
    bkg_dist(MC_rootFile_v231, "test")
    

