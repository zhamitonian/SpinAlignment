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

    #---------------------------------------------
    # sig MC fit
    sigma_initial_vaule = 0.006
    bin_cate2 = [91, 92, 93, 97, 98, 101, 105, 123, 127, 132, 145, 150, 151, 152, 154, 157]
    if bin_index in bin_cate2:
        sigma_initial_vaule = 0.009
    elif bin_index in [80, 135]:
        sigma_initial_vaule = 0.005
    if bin_index in [132, 140, 153,159]:
        sigma_initial_vaule = 0.008
    if bin_index in [155, 156]:
        sigma_initial_vaule = 0.007
    print(f"Bin {bin_index}: initial sigma value = {sigma_initial_vaule}")

    MC_output_dir = f"./fit_results/MCsig_fit/bin_{bin_index:03d}"
    MC_log_file = MC_output_dir + ".log"
    MC_rootFile = R.TFile(f"../Ks_null_test/fit_results/MCsig_fit/temp_bin_{bin_index}.root", "READ")
    MC_tree = MC_rootFile.Get("event")

    pdf_config = FitDefinition([Variable("Ks_M", 0.47, 0.520, 200)],
                    [PDFSpec("DSCB", "Ks_M", "crystal_ball",
                                {"mean":(0.497,0.496,0.498,"mean"),
                                "sigma":(f"sigma[{sigma_initial_vaule},0.001,0.03]*(1-a[0, -1, 1])"),
                                "alpha":(1.5,0.01,5, "alpha_l"),
                                "n":(2.0,0.01,100, "nl"),
                                "n_right":(2.0,0.01,100, "nr"),
                                "alpha_right":(1.5,0.01,5, "alpha_r"),
                                "sigma_right":"sigma*(1+a)"}),
                    PDFSpec("gauss", "Ks_M", "gaussian",
                            {"mean":"mean", "sigma":"sigma * k[1, 0.001, 100]"}),],
                    #model = "SUM(frac[0.3,0,0.8]*DSCB, gauss)")
                    model = "SUM(frac[0.3,0,0.5]*DSCB, gauss)")
    plot_config = PlotConfiguration(plot_config={
        "xlabel" : {"Ks_M" : "M_{#pi^{+}#pi^{-}} (GeV/c^{2})"},
        "components" : {"model" : {"label" : "Total Fit", "color" : 4},
                        "DSCB": {"label" : "DSCB", "color" : 2, "style" : 4, "width" : 3},
                        "gauss": {"label" : "Gaussian", "color" : R.kGreen + 2, "style" : 7, "width" : 3},},
        "legend" : {"extra_text" : range_txt,}, "logy" : True, "show_pull" : False})

    dataset_config = DatasetConfig(binned_fit=True, target_branch=["Ks_M"], perform_splot=False,weight_branch="Ks_weight")
    fit_config = FitterConfig(two_step_fit=True, use_minos=False)

    fitter_mc = GenericFit(MC_tree, MC_output_dir, log_file=MC_log_file, fit_definition=pdf_config, 
                           dataset_config=dataset_config, plot_config=plot_config, fitter_config=fit_config)
    result_mc, _ = fitter_mc.run() 
    
    #----------------------------------------------
    # data fit
    n_entry = tree.GetEntries()
    alpha_devi = 2.0
    n_devi = 3.0
    bkg_range = 1

    pdf_config = FitDefinition([Variable("Ks_M", 0.47, 0.520, 200)],
    #pdf_config = FitDefinition([Variable("Ks_M", 0.47, 0.528, 232)],
                [PDFSpec("DSCB", "Ks_M", "crystal_ball",
                            {"mean":("mean[0.497,0.496,0.498] + diff1[0, -0.01, 0.01]"),
                            "sigma":("sigma[0.006,0.001,0.03] * A[1,0.7,1.3] *(1-a[0, -1, 1])"),
                            "alpha":(1.5,0.01,5, "alpha_l"),
                            "n":(2.0,0.01,100, "nl"),
                            "n_right":(2.0,0.01,100, "nr"),
                            #"n":(f"nr[2.0,0.01,100] * nr_div[1, {1/n_devi}, {n_devi}]"),
                            #"n_right":(f"nl[2.0,0.01,100] * nl_div[1, {1/n_devi}, {n_devi}]"),
                            "alpha_right":(1.5,0.01,5, "alpha_r"),
                            "sigma_right":"sigma*A*(1+a)"}),
                PDFSpec("gauss", "Ks_M", "gaussian",
                        {"mean":"mean + diff2[0, -0.01, 0.01]", "sigma":"sigma * A * k[1, 0.001, 100]"}),
                PDFSpec("bkg", "Ks_M", "chebychev", {"order":1, "coef1": (0, -bkg_range, bkg_range)}),],
                model = f"SUM(nsig[{n_entry}, {0.5 * n_entry}, {10 * n_entry}]*SUM(frac[0.3,0,0.8]*DSCB, gauss), nbkg[{0.05 * n_entry},0, {3 * n_entry}]*bkg)")

    plot_config = PlotConfiguration(plot_config={
        "xlabel" : {"Ks_M" : "M_{#pi^{+}#pi^{-}} (GeV/c^{2})"},
        "components" : {"model" : {"label" : "Total Fit", "color" : 4},
                        "DSCB": {"label" : "DSCB", "color" : 2, "style" : 4, "width" : 3},
                        "gauss": {"label" : "Gaussian", "color" : R.kGreen + 2, "style" : 7, "width" : 3},
                        "bkg": {"label" : "Background", "color" : R.kAzure + 1, "style" : 3, "width" : 3},},
        "legend" : {"extra_text" : range_txt,}, 
        "show_pull" :False, "logy":True})
    mc_constrins = (result_mc, ["a", "alpha_r", "nl", "nr", "alpha_l", "k", "frac" , "sigma", "mean"])
    dataset_config = DatasetConfig(binned_fit= binned_fit,  target_branch=["Ks_M"],
                                   weight_branch="Ks_eff_correction", perform_splot=False)
                                    #perform_splot=False)
    fit_config = FitterConfig(two_step_fit=True, use_minos=False)

    fitter_data = GenericFit(tree, output_dir, log_file=log_file,fit_definition=pdf_config,
                    dataset_config=dataset_config, plot_config=plot_config, fitter_config=fit_config, mc_constrains=mc_constrins)
    fitter_data.run()
    result_data, yield_results = fitter_data.run()
    nsig, nsig_err = yield_results["nsig"], yield_results["nsig_err"]

    """
    chi2_var = R.RooChi2Var("chi2", "chi2", fitter_data.model, fitter_data.dataset)
    chi2 = chi2_var.getVal()
    n_float = result_data.floatParsFinal().getSize()
    ndf = fitter_data.dataset.numEntries() - n_float  
    chi2_ndf = chi2 / ndf
    chi2_ndf_err = (2.0 / ndf) ** 0.5  
    print(f"Fit result for bin {bin_index}: nsig = {nsig:.2f} ± {nsig_err:.2f}, chi2/ndf = {chi2:.1f}/{ndf} = {chi2_ndf:.2f} ")
    """

    return result_data , nsig, nsig_err 


def fit():
    quick_fit = QUICK_FIT(wrapper_func, {"Ks_z": (20, 0, 1), "Ks_helicity_angle": (10, -1, 1)})
    quick_fit.parse_arguments()

if __name__ == "__main__":
    fit()