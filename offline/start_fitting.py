#!/usr/bin/env python3

from FIT import QUICK_FIT, perform_resonance_fit
import ROOT as R
import numpy as np

def fit():
    #quick_fit = QUICK_FIT(perform_resonance_fit, [("phi_xp", 0, 1, 20), ("phi_helicity_angle", -1, 1, 5), ("phi_thrust_pt", 0, 1, 5 )] , binned_fit=False)
    #quick_fit = QUICK_FIT(perform_resonance_fit, [("phi_xp", 0, 1, 20), ("phi_helicity_angle", -1, 1, 10)] )
    #quick_fit = QUICK_FIT(perform_resonance_fit, ("phi_xp", 0, 1, 40))
    
    #pt_bins = np.concatenate((np.linspace(0, 1, num=40, endpoint=False), 
    #                          np.linspace(1, 1.5, num=10, endpoint=False), 
    #                          np.linspace(1.5, 2.5, num=6)))
    #quick_fit = QUICK_FIT(perform_resonance_fit, ("phi_thrust_pt", list(pt_bins)))

    #bin_boundaries = [0.0, 0.0626, 0.125, 0.1876, 0.25, 0.3126, 0.375, 0.4376, 0.5, 0.5625, 0.6378, 0.7285, 0.8378, 0.9694, 1.1279, 1.319, 1.5491, 1.8263, 2.1602, 2.5625]
    #quick_fit = QUICK_FIT(perform_resonance_fit, ("phi_thrust_pt", bin_boundaries))

    quick_fit = QUICK_FIT(perform_resonance_fit, ("phi_thrust_pt", 0, 2.5, 10), binned_fit=True)

    quick_fit.parse_arguments()

def test():
    file = R.TFile("./temp_bin_25.root", "READ")
    tree = file.Get("event")
    perform_resonance_fit(tree, output_dir="./test_images/", binned_fit=True ,branches_name=["phi_M"])

if __name__ == "__main__":
    fit()
    #test()

"""
cmd example:
./start_fitting.py -i ./draw/Data_fitting/data_processed.root -od ./images/ -BrN "phi_M" --batch  
"""