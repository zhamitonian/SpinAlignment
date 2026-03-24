#!/usr/bin/env python3

from FIT import QUICK_FIT, perform_resonance_fit
import ROOT as R

def fit():
    quick_fit = QUICK_FIT(perform_resonance_fit, [("phi_z", 0, 1, 20), ("phi_helicity_angle", -1, 1, 10)] )
    quick_fit.parse_arguments()

if __name__ == "__main__":
    fit()

"""
cmd example:
./start_fitting.py -i ./draw/Data_fitting/data_processed.root -od ./images/ -BrN "phi_M" --batch  
"""