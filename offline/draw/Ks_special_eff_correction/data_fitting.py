#!/usr/bin/env python3
from FIT import QUICK_FIT
import ROOT as R
import os

from FIT.generic_fit import *

def wrapper_func(tree, output_dir, log_file, bin_fit_range, binned_fit):

    ranges = bin_fit_range.split(';')
    range_txt = []
    for range_str in ranges:
        range_str = range_str.strip()
        range_txt.append(range_str)
    bin_index = int(log_file.split("_")[-1].split(".")[0])

    n_entry = tree.GetEntries()
    bkg_range = 1

    pdf_config = FitDefinition([Variable("Ks_M", 0.47, 0.520, 200)],
                [PDFSpec("DSCB", "Ks_M", "crystal_ball",
                        {"mean":(0.497,0.496,0.498,"mean"),
                            "sigma":(f"sigma[0.01,0.001,0.03]*(1-a[0, -1, 1])"),
                            "alpha":(1.5,0.01,5, "alpha_l"),
                            "n":(2.0,0.01,100, "nl"),
                            "n_right":(2.0,0.01,100, "nr"),
                            "alpha_right":(1.5,0.01,5, "alpha_r"),
                            "sigma_right":"sigma*(1+a)"}),
                PDFSpec("gauss", "Ks_M", "gaussian",
                        {"mean":"mean", "sigma":"sigma * k[1, 0.001, 100]"}),
                PDFSpec("bkg", "Ks_M", "chebychev", {"order":1, "coef1": (0, -bkg_range, bkg_range)}),],
                model = f"SUM(nsig[{n_entry}, {0.5 * n_entry}, {1.5 * n_entry}]*SUM(frac[0.3,0,0.8]*DSCB, gauss), nbkg[{0.05 * n_entry},0, {0.1 * n_entry}]*bkg)")

    plot_config = PlotConfiguration(plot_config={
        "xlabel" : {"Ks_M" : "M_{#pi^{+}#pi^{-}} (GeV/c^{2})"},
        "components" : {"model" : {"label" : "Total Fit", "color" : 4},
                        "DSCB": {"label" : "DSCB", "color" : 2, "style" : 4, "width" : 3},
                        "gauss": {"label" : "Gaussian", "color" : R.kGreen + 2, "style" : 7, "width" : 3},
                        "bkg": {"label" : "Background", "color" : R.kAzure + 1, "style" : 3, "width" : 3},},
        "legend" : {"extra_text" : range_txt,}, 
        "show_pull" :False, "logy":True})
    dataset_config = DatasetConfig(binned_fit= binned_fit,  target_branch=["Ks_M"], weight_branch="Ks_eff_correction",perform_splot=False)
    fit_config = FitterConfig(two_step_fit=True, use_minos=False)

    fitter_data = GenericFit(tree, output_dir, log_file=log_file,fit_definition=pdf_config,
                    dataset_config=dataset_config, plot_config=plot_config, fitter_config=fit_config)
    fitter_data.run()
    result_data, yield_results = fitter_data.run()
    nsig, nsig_err = yield_results["nsig"], yield_results["nsig_err"]

    return result_data , nsig, nsig_err 


def fit():
    quick_fit = QUICK_FIT(wrapper_func, {"Ks_z": (20, 0, 1), "Ks_helicity_angle": (10, -1, 1)})
    quick_fit.parse_arguments()

if __name__ == "__main__":
    fit()