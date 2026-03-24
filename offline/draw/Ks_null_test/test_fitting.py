#!/usr/bin/env python3
from FIT import QUICK_FIT
import ROOT as R

from FIT.generic_fit import *

def wrapper_func(tree, output_dir, log_file, bin_fit_range, branches_name, binned_fit):

    ranges = bin_fit_range.split(';')
    range_txt = []
    for range_str in ranges:
        range_str = range_str.strip()
        range_txt.append(range_str)

    plot_config = PlotConfiguration(plot_config={
        "xlabel" : {"Ks_M" : "M_{#pi^{+}#pi^{-}} (GeV/c^{2})"},
        "components" : {"model" : {"label" : "Total Fit", "color" : 4},
                        "DSCB": {"label" : "DSCB", "color" : 2, "style" : 4, "width" : 3},
                        "gauss": {"label" : "Gaussian", "color" : R.kGreen + 2, "style" : 7, "width" : 3},},
        "legend" : {"extra_text" : range_txt,}})
    dataset_config = DatasetConfig(binned_fit= binned_fit, branches_name=branches_name, perform_splot=False,weight_branch="Ks_weight")
    fit_config = FitterConfig(two_step_fit=True, use_minos=False)

    pdf_config = FitDefinition([Variable("Ks_M", 0.47, 0.52, 200)],
            [PDFSpec("DSCB", "Ks_M", "crystal_ball",
                        {"mean":(0.497,0.496,0.498,"mean"),
                        "sigma":(f"sigma[0.006,0.001,0.03]*(1-a[0, -1, 1])"),
                        "alpha":(1.5,0.01,5),
                        "n":(2,0.01,100),
                        "n_right":(2.0,0.01,100),
                        "alpha_right":(1.5,0.01,5),
                        "sigma_right":"sigma*(1+a)"}),
            PDFSpec("gauss", "Ks_M", "gaussian",
                    {"mean":"mean", 
                        "sigma":"sigma * k[1, 0.001, 100]"}),
            PDFSpec("gauss2", "Ks_M", "gaussian",
                    {"mean":(0.001, -0.005, 0.005),
                        "sigma":(0.0006, 0.0005, 0.005)}),
            PDFSpec("bkg", "Ks_M", "chebychev", {"order":1})],
            model = "SUM(nsig[20000,0,4000000] * SUM(frac[0.3,0,0.8]*DSCB, gauss), nbkg[5000,0,10000] * bkg)")

    fitter = GenericFit(tree, output_dir, log_file=log_file, fit_definition=pdf_config,
                    dataset_config=dataset_config, plot_config=plot_config, fitter_config=fit_config)
    result, fit_results = fitter.run() 
    nsig, nsig_err = fit_results["nsig"], fit_results["nsig_err"]
    print(f"Fit result: nsig = {nsig} ± {nsig_err}")

    return result , nsig, nsig_err


def fit():
    quick_fit = QUICK_FIT(wrapper_func, [("Ks_z", 0, 1, 20), ("Ks_helicity_angle", -1, 1, 10)])
    quick_fit.parse_arguments()

if __name__ == "__main__":
    fit()