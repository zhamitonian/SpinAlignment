#!/usr/bin/env python3
from FIT import QUICK_FIT
import ROOT as R
import os

from FIT.generic_fit import *
from math import pi

def wrapper_func(tree, output_dir, log_file, var_name):

    #----------------------------------------------
    if var_name.startswith("phi"):
        mass , width = 1.019461, 0.004266
        xmin, xmax, nbins = 0.99, 1.06, 140
        target_branch = ["phi_M"]
    elif var_name.startswith("Kstar"):
        mass , width = 0.89555, 0.0473
        xmin, xmax, nbins = 0.69, 1.09, 200
        target_branch = ["Kstar_M", "Kstar_bar_M"]
    elif var_name.startswith("D0"):
        mass , width = 1.8648, 0.001
        xmin, xmax, nbins = 1.80, 1.92, 120
        xmin, xmax, nbins = 1.80, 1.92, 120
        target_branch = ["D0_M", "D0_bar_M"]
    else:
        raise ValueError(f"Unknown variable {var_name} for fitting")

    var_name = "a"
    pdf_config = FitDefinition([Variable(var_name, xmin, xmax, nbins)],
                    [PDFSpec("sig", var_name, "bw_gauss", 
                             {"mass":mass, "width":width , 
                              "resolution":(0.001, 0.0001, 0.01), 
                              "reso_mean": (0, -0.005, 0.005)}),
                    PDFSpec("bkg", var_name, "chebychev", {"order":1, "coef1":(0, -10, 10)}),
                    ],
                    model = "SUM(nsig[10000,0,1e6]*sig, nbkg1[50, 0, 1e6]*bkg)")
    plot_config = PlotConfiguration(plot_config={
        "xlabel" : {"Kstar_M" : "M_{K^{#pm}#pi^{#mp}} (GeV/c^{2})" ,"phi_M" : "M_{K^{+}K^{-}} (GeV/c^{2})", "D0_M": "M_{K^{#mp}#pi^{#pm}} (GeV/c^{2})"},
        "components" : {"model" : {"label" : "Total Fit", "color" : 4},
                        "sig": {"label" : "Signal", "color" : 2, "style" : 4, "width" : 3},
                        "bkg": {"label" : "Bkg", "color" : 7, "style" : 2, "width" : 3},},
        "legend" : {"extra_text" : ["Belle data @ 9.43 GeV", "L_{int} = 1.816 fb^{-1}"]},
        "show_pull" :False, "logy":False})
    dataset_config = DatasetConfig(binned_fit= True, target_branch=target_branch, perform_splot=False)
    fit_config = FitterConfig(two_step_fit=True, use_minos=False)
    fitter_data = GenericFit(tree, output_dir, log_file=log_file,fit_definition=pdf_config,
                    dataset_config=dataset_config, plot_config=plot_config, fitter_config=fit_config)
    
    result_data, yield_results = fitter_data.run()
    nsig, nsig_err = yield_results["nsig"], yield_results["nsig_err"]

    return result_data , nsig, nsig_err 


def fit():
    rootFile = "./reco.root"
    file = R.TFile(rootFile)
    tree = file.Get("event")
    for par in ["D0"]:
        wrapper_func(tree, "./fit_results/",f"./fit_results/{par}_fit.log", var_name=f"{par}_M")

if __name__ == "__main__":
    fit()