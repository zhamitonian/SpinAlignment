#!/usr/bin/env python3
from FIT import QUICK_FIT
import ROOT as R
import os
from math import pi

from FIT.generic_fit import *

def wrapper_func(tree, output_dir, log_file, bin_fit_range, binned_fit):

    ranges = bin_fit_range.split(';')
    range_txt = []
    for range_str in ranges:
        range_str = range_str.strip()
        range_txt.append(range_str)
    bin_index = int(log_file.split("_")[-1].split(".")[0])

    # data fit
    n_entry = tree.GetEntries()
    #----------------------------------------------
    mean , width = 1.0195, 0.004249

    pdf_config = FitDefinition([Variable("phi_M", 1, 1.06, 120)],
                [PDFSpec("sig", "phi_M", "bw_gauss",
                            {"mass": (mean, mean-0.005, mean+0.05), "width": width, "reso_mean":0, 
                             "resolution" : (0.001, 0.0002, 0.01),}),
                PDFSpec("bkg", "phi_M", "chebychev", {"order":2, "coef1": (0, -10, 10), "coef2":(0, -1, 1)}),],
                model = f"SUM(nsig[{0.5 * n_entry}, 0, {n_entry}] *sig, nbkg[{0.5 * n_entry},0, {n_entry}]*bkg)")

    plot_config = PlotConfiguration(plot_config={
        "xlabel" : {"phi_M" : "M_{K^{+}K^{-}} (GeV/c^{2})"},
        "components" : {"model" : {"label" : "Total Fit", "color" : 4},
                        "sig": {"label" : "Signal", "color" : 2, "style" : 4, "width" : 3},
                        "bkg": {"label" : "Background", "color" : R.kAzure + 1, "style" : 3, "width" : 3},},
        "legend" : {"extra_text" : range_txt,}, 
        "show_pull" :False, "logy":False})
    dataset_config = DatasetConfig(binned_fit= binned_fit, target_branch=["phi_M"], perform_splot=False)
    fit_config = FitterConfig(two_step_fit=True, use_minos=False)

    fitter = GenericFit(tree, output_dir, log_file=log_file,fit_definition=pdf_config,
                    dataset_config=dataset_config, plot_config=plot_config, fitter_config=fit_config)
    fitter.run()
    result_data, yield_results = fitter.run()
    nsig, nsig_err = yield_results["nsig"], yield_results["nsig_err"]

    return result_data , nsig, nsig_err 


def fit(): 
    quick_fit = QUICK_FIT(wrapper_func, {"phi_z": (20, 0, 1), "phi_helicity_phi": (20, 0, pi/2)})
    quick_fit.parse_arguments()

if __name__ == "__main__":
    fit()