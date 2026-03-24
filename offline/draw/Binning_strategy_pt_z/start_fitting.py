#!/usr/bin/env python3

from FIT import QUICK_FIT, perform_resonance_fit
import ROOT as R
import numpy as np
from functools import partial

def fit():
    func = partial(perform_resonance_fit, which_bkg = 1, bkg_order = 1)#, if_save = True)

    pt_exp_binning = [0.0, 0.125, 0.25, 0.375, 0.5, 0.6611, 0.8688, 1.1366, 1.4817, 1.9265, 2.5]
    #quick_fit = QUICK_FIT(perform_resonance_fit, ("phi_thrust_pt", pt_exp_binning), binned_fit=True)
    #quick_fit = QUICK_FIT(perform_resonance_fit, ("phi_thrust_pt", 0, 2.5, 10), binned_fit=True)

    #z_exp_binning = [0.2, 0.275, 0.35, 0.425, 0.5, 0.5773, 0.6569, 0.7389, 0.8233, 0.9104, 1.0]
    #quick_fit = QUICK_FIT(func, ("phi_z",0.2, 1, 10))
    #quick_fit = QUICK_FIT(perform_resonance_fit, ("phi_z", z_exp_binning), binned_fit=True)

    #quick_fit = QUICK_FIT(perform_resonance_fit, [("phi_z", 0.2, 1, 10), ("phi_thrust_pt", pt_exp_binning), ("phi_helicity_angle", -1, 1, 10)])
    quick_fit = QUICK_FIT(func, [("phi_z", 0.2, 1, 10), ("phi_thrust_pt", pt_exp_binning), ("phi_helicity_angle", -1, 1, 10)])

    quick_fit.parse_arguments()

if __name__ == "__main__":
    fit()

"""
cmd example:
./start_fitting.py -i ./draw/Data_fitting/data_processed.root -od ./images/ -BrN "phi_M" --batch  
"""