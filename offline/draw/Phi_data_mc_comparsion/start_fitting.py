#!/usr/bin/env python3

from FIT import QUICK_FIT, perform_resonance_fit
import ROOT as R
import numpy as np
from functools import partial

def fit():
    #func = partial(perform_resonance_fit, which_bkg = 1, bkg_order = 1)#, if_save = True)
    #func = partial(perform_resonance_fit, which_bkg = 0, bkg_order = 2, fit_config = (1,1.04, 40))#, if_save = True)
    # bin [0, 39, 50]
    func = partial(perform_resonance_fit, which_bkg = 1, bkg_order = 2, fit_config = (1,1.03, 30))#, if_save = True)

    pt_exp_binning = [0.0, 0.125, 0.25, 0.375, 0.5, 0.6611, 0.8688, 1.1366, 1.4817, 1.9265, 2.5]
    quick_fit = QUICK_FIT(func, [("phi_z", 0.2, 1, 10), ("phi_thrust_pt", pt_exp_binning), ("phi_helicity_angle", -1, 1, 10)])

    quick_fit.parse_arguments()

if __name__ == "__main__":
    fit()

"""
cmd example:
./start_fitting.py -i ./draw/Data_fitting/data_processed.root -od ./images/ -BrN "phi_M" --batch  
"""