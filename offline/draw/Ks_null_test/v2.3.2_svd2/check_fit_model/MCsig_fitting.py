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
                        #"bkg": {"label" : "Background", "color" : R.kGray + 2, "style" : 3, "width" : 3},},
        "legend" : {"extra_text" : range_txt,}})
    dataset_config = DatasetConfig(binned_fit= binned_fit, branches_name=branches_name, perform_splot=False,weight_branch="Ks_weight")
    fit_config = FitterConfig(two_step_fit=True, use_minos=False)

    #hist = R.TH1F("", ";sigma(GeV/c^{2});#chi^{2}/ndof", 14, 0.001, 0.015) 
    hist = R.TH1F("", ";sigma(GeV/c^{2});#chi^{2}/ndof", 8, 0, 8) 
    pad_width = len(str(hist.GetNbinsX()))
    for i in range(hist.GetNbinsX()):
        pdf_config = FitDefinition([Variable("Ks_M", 0.47, 0.52, 200)],
                [PDFSpec("DSCB", "Ks_M", "crystal_ball",
                            {"mean":(0.497,0.496,0.498,"mean"),
                            #"sigma":(f"sigma[{i*0.001 + 0.001},0.001,0.03]*(1-a[0, -1, 1])"),
                            "sigma":("sigma[0.006,0.001,0.03]*(1-a[0, -1, 1])"),
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
                model = "SUM(frac[0.3,0,0.8]*DSCB, gauss)")

        log_file_i = log_file.split(".")[0] + f"_{i:0{pad_width}d}.log"
        print(log_file_i)
        fitter = GenericFit(tree, output_dir + f"_{i:0{pad_width}d}", log_file=log_file_i, fit_definition=pdf_config,
                    dataset_config=dataset_config, plot_config=plot_config, fitter_config=fit_config)
        result, fit_results = fitter.run() 

        chi2_var = R.RooChi2Var("chi2", "chi2", fitter.model, fitter.dataset)
        chi2 = chi2_var.getVal()
        n_float = result.floatParsFinal().getSize()
        ndf = fitter.dataset.numEntries() - n_float  
        chi2_ndf = chi2 / ndf
        chi2_ndf_err = (2.0 / ndf) ** 0.5  

        hist.SetBinContent(i+1, chi2_ndf)
        hist.SetBinError(i+1, chi2_ndf_err)

    bin_index = log_file.split("_")[-1].split(".")[0]
    from DRAW import style_draw,HistStyle
    style_draw([hist], f"test/bin{bin_index}_chi2.png", styles= [HistStyle.error_bars(1)], y_min=0, y_max=5, use_user_y_range=True)
    
    return result , 0, 0


def fit():
    quick_fit = QUICK_FIT(wrapper_func, [("Ks_z", 0, 1, 20), ("Ks_helicity_angle", -1, 1, 10)])
    quick_fit.parse_arguments()

if __name__ == "__main__":
    fit()