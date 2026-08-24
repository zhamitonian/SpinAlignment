#!/usr/bin/env python3

import ROOT as R
import sys
sys.path.append("/gpfs/group/belle/users/wangz/Work/SpinAlignment/offline/draw/")
from common.physics.comparer import process_single_bin
from common.config.comparison import DEFAULT_KS_WEIGHT_CONFIG
from common.config.binning import DEFAULT_KS_BINNING
from common.utils.latex_summary import generate_latex_document

comp_plots = ["Ks_p", "Ks_costheta", "Ks_phi", "pip_p", "pim_p", "pip_costheta", "pim_costheta", "pip_phi", "pim_phi"]

title = r"Data/MC Comparison for $K_s$ Spin Alignment Analysis"

vec_branches = ["Ks_M", "Ks_p", "Ks_costheta", "Ks_phi", "Ks_weight", "Ks_z", "Ks_helicity_angle"]
vec_branches += ["pip_p", "pim_p", "pip_costheta", "pim_costheta", "pip_phi", "pim_phi", "pip_pid_weight", "pim_pid_weight"]

def wrapper(flat_bin_idx):
    data_rootFile_dir = "../fit_results/data_fit"
    splot_rootFile_dir = data_rootFile_dir 
    MC_rootFile = "../rootFiles/sig_isSignal_v2.3.3.root"
    output_dir = "../images/data_mc_comparing/" 

    process_single_bin(flat_bin_idx, DEFAULT_KS_BINNING, comp_plots, 
                       data_rootFile_dir, splot_rootFile_dir, MC_rootFile, 
                       output_dir, vec_branches, DEFAULT_KS_WEIGHT_CONFIG, reweight_var="pip_costheta")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "--bin" and len(sys.argv) > 2:
            bin_idx = int(sys.argv[2])
            wrapper(bin_idx)
        if sys.argv[1] == "--gen_summary":
            generate_latex_document(DEFAULT_KS_BINNING, comp_plots, title,
                                    #fit_plot_dir= "../images/fit_dt_e55_v2.1.0",
                                    comparison_plot_dir =  "../images/data_mc_comparing/" ,
                                    valid_bins= list(range(150, 160)))
        else:
            print("Usage:")
            print("  python3 start_comparing.py --bin <idx>   # Process single bin")
