#!/usr/bin/env python3

import ROOT as R
from rho00_data_mc_comparsion import process_single_bin

import numpy as np
import sys

binning_config = {"Ks_z": np.linspace(0., 1, 21).tolist(),
                      "Ks_helicity_angle": np.linspace(-1.0, 1.0, 11).tolist()
                     }
plot_labels = {"Ks_p":r"$p(K_s)$",
               "Ks_costheta":r"$\cos\theta(K_s)$",
               "Ks_phi":r"$\varphi(K_s)$",
               "pip_p":r"$p(\pi^+)$",
               "pim_p":r"$p(\pi^-)$",
               "pip_costheta":r"$\cos\theta(\pi^+)$",
               "pim_costheta":r"$\cos\theta(\pi^-)$",
               "pip_phi":r"$\varphi(\pi^+)$",
               "pim_phi":r"$\varphi(\pi^-)$"
              }

var_labels = {
    "Ks_z" : r"$z$",
    "Ks_helicity_angle" : r"$cos\theta^{\star}$"
}
title = r"Data/MC Comparison for $K_s$ Spin Alignment Analysis"

def wrapper(flat_bin_idx):

    data_rootFile_dir = "/gpfs/home/belle2/wangz/Work/SpinAlignment/offline/draw/Ks_null_test/images/fit_dt_e55_v2.1.0"
    splot_rootFile_dir = data_rootFile_dir 
    MC_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.1.0_qqbar/exp55_reco_truth_processed.root" 
    output_dir = "../images/data_mc_comparing_e55_v2.1.0" 
    var_config = {"Ks_M": (100, 0.47, 0.52, ";M_{K_{s}};[MeV]"), 
                    "Ks_p": (100, 0, 5, ";p(K_{s});[GeV/c]"),
                    "Ks_costheta": (20, -1, 1, ";cos#theta(K_{s});[]"),
                    "Ks_phi" : (20, -3.14, 3.14, ";#varphi(K_{s});[rad]"),
                    "pip_p": (90, 0, 4.5, ";p(#pi^{+});[MeV]"),
                    "pim_p": (90, 0, 4.5, ";p(#pi^{-});[MeV]"),
                    "pip_costheta": (20, -1, 1, ";cos#theta(#pi^{+});[]"),
                    "pim_costheta": (20, -1, 1, ";cos#theta(#pi^{-});[]"),
                    "pip_phi": (20, -3.14, 3.14, ";#phi(#pi^{+});[rad]"),
                    "pim_phi": (20, -3.14, 3.14, ";#phi(#pi^{-});[rad]"),
                    }

    process_single_bin(flat_bin_idx, binning_config, var_config, data_rootFile_dir, splot_rootFile_dir, MC_rootFile, output_dir)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "--bin" and len(sys.argv) > 2:
            # Process single bin (called by LSF job)
            bin_idx = int(sys.argv[2])
            wrapper(bin_idx)
        if sys.argv[1] == "--gen_summary":
            from generate_latex_summary import generate_latex_document
            generate_latex_document(binning_config, plot_labels, var_labels, title,
                                    #fit_plot_dir= "../images/fit_dt_e55_v2.1.0",
                                    comparison_plot_dir= "../images/data_mc_comparing_e55_v2.1.0/",
                                    valid_bins= [i+10 for i in range(140)])
        else:
            print("Usage:")
            print("  python3 start_comparing.py --bin <idx>   # Process single bin")
