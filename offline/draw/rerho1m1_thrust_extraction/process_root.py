#!/usr/bin/env python3

import sys
import os
import argparse
import datetime  # must be before ROOT: ROOT's import hook breaks datetime capsule
import numpy     # must be before ROOT: ROOT's import hook breaks numpy C extensions
import ROOT as R
from OFFLINE_PROCESS import RDF_process
from DRAW import style_draw, HistStyle
import glob
from math import pi
from typing import Optional, Tuple
sys.path.append("/home/belle2/wangz/Work/SpinAlignment/offline/draw/belle_kid_eff/")
from get_pid_eff_weight import get_weight
    

### reco and truth func are copied from phiSA_processor.py
def reco(df: R.RDataFrame, hasPhi: bool = False, truth_match: bool = False):
    """
    Inclusive phi / rho00 measurement: reconstruction-level processing.
    """
    phi_vars = []
    phi_vars_toSave = []
    kaon_vars = ["px_cms", "py_cms", "pz_cms", "E_cms"]
    kaon_vars_toSave = []

    is_MC = "kp_isSignal" in df.GetColumnNames()
    tools = RDF_process()

    # --- MC-only pre-selection ---
    if is_MC:
        kaon_vars += ["isSignal"]
        kaon_vars_toSave += ["isSignal"]
        if truth_match:
            for kaon in ["kp", "km"]:
                for var in kaon_vars:
                    df = df.Redefine(f"{kaon}_{var}", f"{kaon}_{var}[{kaon}_isSignal]")
        elif hasPhi:  # require at least one truth-level phi
            df = df.Filter("n_phi_truth > 0", "has truth phi")
        df = get_weight(df, period="svd2", is_pion=False, cut_value=6, particle_names=("kp", "km"))
        kaon_vars += ["pid_weight"]
        kaon_vars_toSave += ["pid_weight"]


    # --- Build all K+K- pairs and compute invariant mass ---
    df = tools.save_all_pairs(
        df,
        p1_branches=("kp_px_cms", "kp_py_cms", "kp_pz_cms"),
        p2_branches=("km_px_cms", "km_py_cms", "km_pz_cms"),
        mass=(0.493677, 0.493677),
        particle_name=("phi", "kp", "km"),
    )
    phi_vars += ["phi_E", "phi_px", "phi_py", "phi_pz", "kp_index", "km_index", "phi_helicity_angle", "phi_helicity_phi"]
    phi_vars_toSave += ["kp_index", "km_index", "phi_helicity_angle", "phi_helicity_phi"]

    ## Fold azimuthal angle by quadrant: result in [0, pi/2]
    df = df.Redefine("phi_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(phi_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")


    # --- Apply KK mass window (KK threshold ~0.987 GeV) --- cut 1
    df = df.Define("phi_M", "sqrt(phi_E*phi_E - phi_px*phi_px - phi_py*phi_py - phi_pz*phi_pz)")
    phi_vars.append("phi_M")
    phi_vars_toSave.append("phi_M")

    phi_mass_window = (0.99, 1.15)
    df = df.Define("mass_window_pass",
                   f"(phi_M > {phi_mass_window[0]}) && (phi_M < {phi_mass_window[1]})")
    for var in phi_vars:
        df = df.Redefine(var, f"{var}[mass_window_pass]")

    df = df.Filter("phi_M.size() > 0", "has phi after mass window cut")

    # --- Calculate angles related to thrust axis --- theta < 90^\circ cut 2 
    df = tools.calculate_pt_toAxis(
        df,
        particle=("phi_px", "phi_py", "phi_pz"),
        axis=("TMath::ACos(thrust[1])", "thrust[2]"),
        particle_name="phi", axis_name="thrust",
    )
    phi_vars += ["phi_thrust_pt", "phi_thrust_costheta"]
    phi_vars_toSave += ["phi_thrust_pt", "phi_thrust_costheta"]

    df = df.Define("phi_thrust_costheta_pass", "phi_thrust_costheta > 0")
    for var in phi_vars:
        df = df.Redefine(var, f"{var}[phi_thrust_costheta_pass]")

    df = df.Filter("phi_M.size() > 0", "has phi after thrust costheta cut")
    

    # --- calculate other variables

    ## --- thrust , polarization angles --- 
    df = tools.cal_pola_angles(df, ("kp_px_cms", "kp_py_cms", "kp_pz_cms"),
                            ("km_px_cms", "km_py_cms", "km_pz_cms"),
                            (0.493677, 0.493677), ("thrust[1]", "thrust[2]"), ("kp_index", "km_index"),
                            ("phi_thrust_helicity_costheta", "phi_thrust_helicity_phi"))
    phi_vars += ["phi_thrust_helicity_costheta", "phi_thrust_helicity_phi"]
    phi_vars_toSave += ["phi_thrust_helicity_costheta", "phi_thrust_helicity_phi"]

    df = df.Redefine("phi_thrust_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(phi_thrust_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")
    
    ## --- MC-only: PID efficiency weight per pair ---
    if is_MC:
        df = df.Define(
            "phi_weight",
            "ROOT::VecOps::Take(kp_pid_weight, kp_index) * "
            "ROOT::VecOps::Take(km_pid_weight, km_index)",
        )
        phi_vars.append("phi_weight")
        phi_vars_toSave.append("phi_weight")

    ## --- other Kinematic variables ---
    df = tools.convert_cartesian_to_spherical(df, particles=["phi"])
    df = tools.convert_cartesian_to_spherical(df, particles=["kp", "km"],
                                               px_branch="px_cms", py_branch="py_cms", pz_branch="pz_cms")
    df = df.Define("phi_z",  "2*phi_E/sqrts")
    phi_vars += ["phi_costheta", "phi_phi", "phi_p", "phi_z"]
    phi_vars_toSave += ["phi_costheta", "phi_phi", "phi_p", "phi_z"]
    kaon_vars += ["p", "costheta", "phi"]
    kaon_vars_toSave += ["p", "costheta", "phi"]

    # --- Branches to save ---
    Br2Save = phi_vars_toSave + \
        ["kp_" + var for var in kaon_vars_toSave] + ["km_" + var for var in kaon_vars_toSave] + \
        ["thrust"]

    df.Report().Print()

    return (df, Br2Save)


def truth(df:R.RDataFrame, kaon_inBarrel:bool=False):
    """
    process truth tree, retrieve coordinate generate leavel phi's z, cos*theta ...
    """
    phi_vars = []
    phi_vars_toSave = []
    kaon_vars = ["px", "py", "pz", "E"]
    kaon_vars = ["kp_" + var + "_cms_truth" for var in kaon_vars] + ["km_" + var + "_cms_truth" for var in kaon_vars]

    tools = RDF_process()

    
    # --- pre-selection ---
    if kaon_inBarrel:
        kaon_inBarrel = f"cos(kp_theta*{pi}/180) < 0.842 && cos(kp_theta*{pi}/180) > -0.511 && cos(km_theta*{pi}/180) < 0.842 && cos(km_theta*{pi}/180) > -0.511 && kp_pt_truth > 0.05 && km_pt_truth > 0.05" # type: ignore
        for var in ["px", "py", "pz"]:
            for kaon in ["kp", "km"]:
                df = df.Redefine(f"{kaon}_{var}_cms_truth", f"{kaon}_{var}_cms_truth[{kaon_inBarrel}]")
        

    # --- Build all K+K- pairs and compute invariant mass ---
    df = tools.save_all_pairs(df, 
                         p1_branches=("kp_px_cms_truth","kp_py_cms_truth","kp_pz_cms_truth"),
                         p2_branches=("km_px_cms_truth","km_py_cms_truth","km_pz_cms_truth"),
                         mass=(0.493677, 0.493677),
                         particle_name=("phi","kp","km"),
                         cross_mode= False)
    phi_vars += ["phi_E", "phi_px", "phi_py", "phi_pz", "kp_index", "km_index", "phi_helicity_angle", "phi_helicity_phi"]
    phi_vars_toSave += ["kp_index", "km_index", "phi_helicity_angle", "phi_helicity_phi"]

    ## fold by quadrant : pi/2 - |fmod(phi, pi) - pi/2|, result in [0, pi/2]
    df = df.Redefine("phi_helicity_phi",
                   "TMath::Pi()/2 - abs(fmod(phi_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")
    
    # --- calculate thrust pola angles ---
    df = tools.cal_pola_angles(df, ("kp_px_cms_truth", "kp_py_cms_truth", "kp_pz_cms_truth"),
                            ("km_px_cms_truth", "km_py_cms_truth", "km_pz_cms_truth"),
                            (0.493677, 0.493677), ("thrust_truth[1]", "thrust_truth[2]"), ("kp_index", "km_index"),
                            ("phi_thrust_helicity_costheta", "phi_thrust_helicity_phi"))
    phi_vars += ["phi_thrust_helicity_costheta", "phi_thrust_helicity_phi"]
    phi_vars_toSave += ["phi_thrust_helicity_costheta", "phi_thrust_helicity_phi"]


    # --- calculate other variables
    df = df.Define("phi_M", "sqrt(phi_E*phi_E - phi_px*phi_px - phi_py*phi_py - phi_pz*phi_pz)")
    df = df.Define("phi_z", f"2*phi_E/10.516469955444336")
    phi_vars += ["phi_M", "phi_z"]
    phi_vars_toSave += ["phi_M", "phi_z"]
    

    # --- Branches to save ---
    Br2Save = phi_vars_toSave + kaon_vars
    df.Report().Print()

    return (df, Br2Save)

if __name__ == "__main__":
    qqbarMC_path = "/gpfs/group/belle/users/dues/data_gMC_belle1/PhiSpinAlignment_v2.1.5_qqbar_svd2_duxs/steered.root"
    data_path = "/gpfs/group/belle/users/wangz/data_gMC_belle1/PhiSpinAlignment_v2.1.3_data_svd2/steered.root"

    """"    
    df_qqbar = R.RDataFrame("event", qqbarMC_path)
    df_data = R.RDataFrame("event", data_path)

    df_qqbar, Br2Save_qqbar = reco(df_qqbar)
    df_qqbar.Snapshot("event", "reco_output_qqbar.root", Br2Save_qqbar)

    df_data, Br2Save_data = reco(df_data)
    df_data.Snapshot("event", "reco_output_data.root", Br2Save_data)
    """

    df_qqbar_truth = R.RDataFrame("truth", qqbarMC_path)
    df_qqbar_truth, Br2Save_qqbar_truth = truth(df_qqbar_truth)
    df_qqbar_truth.Snapshot("truth", "truth_output_qqbar.root", Br2Save_qqbar_truth)