#!/usr/bin/env python3
from FIT import QUICK_FIT
import ROOT as R

import sys
sys.path.append("/home/belle2/wangz/Work/ANA_TOOLS/FIT/")
from FIT.generic_fit import * 


def wrapper_func(tree, output_dir, log_file, bin_fit_range, branches_name, binned_fit):

    ranges = bin_fit_range.split(';')
    range_txt = []
    for range_str in ranges:
        range_str = range_str.strip()
        range_txt.append(range_str)

    double_gauss_config = FitDefinition([Variable("Ks_M", 0.475, 0.52, 90)],
                                        [PDFSpec("sig", "Ks_M", "double_gauss", 
                                                 {"same_mean":True,
                                                  "mean1": (0.49761, 0.49561, 0.49861),
                                                  "sigma1":(0.002, 0.0008, 0.004),  # narror gaussian
                                                  "sigma2":(0.004, 0.0025, 0.006),   # wide
                                                  "frac":(0.7, 0.5, 1)
                                                 }),
                                        PDFSpec("bkg", "Ks_M", "chebychev",
                                                {"order":3})],
                                        model = "nsig[10000,0,4000000]*sig + nbkg[5000,0,400000]* bkg")    

    plot_config = PlotConfiguration(plot_config={
        "xlabel" : {"Ks_M" : "M_{#pi^{+}#pi^{-}} (GeV/c^{2})"},
        "components" : {
            "model" : {"label" : "Total Fit", "color" : 4},
            "sig" : {"label" : "Signal", "color" : 2, "style" : 4, "width" : 3},
            "bkg" : {"label" : "Background", "color" : R.kGreen + 2, "style" : 7, "width" : 3},
        },
        "legend" : {"extra_text" : range_txt,
                   }
    }) 

    dataset_config = DatasetConfig(binned_fit = binned_fit,
                                    branches_name = branches_name,    
                                    if_save = False, # save the combined root file
                                    perform_splot = True)

    fitter = GenericFit(tree, output_dir, log_file = log_file,
                        fit_definition= double_gauss_config,
                        dataset_config= dataset_config,
                        plot_config= plot_config
                     )


    result, fit_results = fitter.run()
    nsig , nsig_err = fit_results["nsig"], fit_results["nsig_err"]

    return result, nsig, nsig_err

            

def fit():
    quick_fit = QUICK_FIT(wrapper_func, [("Ks_z", 0, 1, 20), ("Ks_helicity_angle", -1, 1, 10)])
    quick_fit.parse_arguments()

if __name__ == "__main__":
    fit()

"""
cmd example:
./start_fitting.py -i ./draw/Data_fitting/data_processed.root -od ./images/ -BrN "phi_M" --batch  
"""