#!/usr/bin/env python3

from FIT import QUICK_FIT, perform_resonance_fit
import ROOT as R
from functools import partial

#Ks_FIT_WINDOW = (0.4, 0.6, 1000)
#Ks_FIT_WINDOW = (0.467, 0.527, 80)
Ks_FIT_WINDOW = (0.470, 0.525, 110)
Ks_PARTICLE_CONFIG = ("Ks", 0.49761, 0)

def fit():
    wrapper_func = partial(perform_resonance_fit, 
                           particle_config = Ks_PARTICLE_CONFIG,
                           fit_config = Ks_FIT_WINDOW,
                           which_bkg = 0,
                           which_order = 1,
                           which_sig = 2,
                           conv_gauss = False)

    quick_fit = QUICK_FIT(wrapper_func, [("Ks_z", 0, 1, 20), ("Ks_helicity_angle", -1, 1, 10)] )
    quick_fit.parse_arguments()

if __name__ == "__main__":
    fit()