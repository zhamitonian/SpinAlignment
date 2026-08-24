#!/usr/bin/env python3
from FIT.utils import QUICK_FIT
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

    plot_config = PlotConfiguration(plot_config={
        "xlabel" : {"Ks_M" : "M_{#pi^{+}#pi^{-}} (GeV/c^{2})"},
        "components" : {"model" : {"label" : "Total Fit", "color" : 4},
                        "DSCB": {"label" : "DSCB", "color" : 2, "style" : 4, "width" : 3},
                        "gauss": {"label" : "Gaussian", "color" : R.kGreen + 2, "style" : 7, "width" : 3},},
                        #"sig": {"label" : "Gaussian", "color" : R.kGreen + 2, "style" : 7, "width" : 3},},
        "legend" : {"extra_text" : range_txt,}, "logy" :True, "show_pull" : False})
    dataset_config = DatasetConfig(binned_fit= binned_fit, perform_splot=False,weight_branch="Ks_weight")
    fit_config = FitterConfig(two_step_fit=True, use_minos=False)

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

    #hist = R.TH1F("", ";sigma(GeV/c^{2});#chi^{2}/ndof", 10, 0, 10) 
    hist = R.TH1F("", ";sigma(GeV/c^{2});#chi^{2}/ndof", 1, 0, 1) 
    pad_width = len(str(hist.GetNbinsX()))
    for i in range(hist.GetNbinsX()):
        pdf_config = FitDefinition([Variable("Ks_M", 0.47, 0.52, 200)],
                [PDFSpec("DSCB", "Ks_M", "crystal_ball",
                            {"mean":(0.497,0.496,0.498,"mean"),
                            #"sigma":(f"sigma[{i*0.001 + 0.001},0.001,0.03]*(1-a[0, -1, 1])"),
                            "sigma":(f"sigma[{sigma_initial_vaule},0.001,0.03]*(1-a[0, -1, 1])"),
                            "alpha":(1.5,0.01,5),
                            "n":(2,0.01,100),
                            "n_right":(2.0,0.01,100),
                            "alpha_right":(1.5,0.01,5),
                            "sigma_right":"sigma*(1+a)"}),
                PDFSpec("gauss", "Ks_M", "gaussian",
                        {"mean":"mean", "sigma":"sigma * k[1, 0.001, 100]"}),
                PDFSpec("gauss2", "Ks_M", "gaussian",
                        {"mean":"mean", "sigma":"sigma * k2[1, 0.001, 100]"}),
                PDFSpec("bifur_gauss", "Ks_M", "bifur_gauss",
                        {"mean":"mean", "sigma_left":"sigma *kl[1, 0.001, 100]", "sigma_right":"sigma * kr[1, 0.001, 100]"}),
                PDFSpec("bkg", "Ks_M", "chebychev", {"order":1}),
                PDFSpec("sig", "Ks_M", "composite", {"formula" : "SUM(frac[0.3,0,0.8]*DSCB, gauss)"})],
                model = "SUM(frac[0.3,0,0.8]*DSCB, gauss)")
                #model = "sig")
                #model = "SUM(frac[0.3,0,0.8]*DSCB, frac2[0.5,0.3, 0.9]*gauss, gauss2)")

        log_file_i = os.path.splitext(log_file)[0] + f"_{i:0{pad_width}d}.log"
        fitter = GenericFit(tree, output_dir + f"_{i:0{pad_width}d}", log_file=log_file_i, fit_definition=pdf_config,
                    dataset_config=dataset_config, plot_config=plot_config, fitter_config=fit_config)
        result, fit_results = fitter.run() 
        """
        chi2_var = R.RooChi2Var("chi2", "chi2", fitter.model, fitter.dataset)
        chi2 = chi2_var.getVal()
        n_float = result.floatParsFinal().getSize()
        ndf = fitter.dataset.numEntries() - n_float  
        chi2_ndf = chi2 / ndf
        chi2_ndf_err = (2.0 / ndf) ** 0.5  

        hist.SetBinContent(i+1, chi2_ndf)
        hist.SetBinError(i+1, chi2_ndf_err)
        """

    bin_index = log_file.split("_")[-1].split(".")[0]
    from DRAW import style_draw,HistStyle
    style_draw([hist], f"test/bin{bin_index}_chi2.png", styles= [HistStyle.error_bars(1)], y_min=0, y_max=5, use_user_y_range=True)
    
    return result , 0, 0


def fit():
    quick_fit = QUICK_FIT(wrapper_func, {"Ks_z": (20, 0, 1), "Ks_helicity_angle": (10, -1, 1)})
    quick_fit.parse_arguments()

if __name__ == "__main__":
    fit()