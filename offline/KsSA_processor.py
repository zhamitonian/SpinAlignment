#!/usr/bin/env python3

import sys
import os
import argparse
import ROOT as R
from OFFLINE_PROCESS import RDF_process
from DRAW import style_draw, HistStyle
import glob
from math import pi
from typing import Optional, Tuple
sys.path.append("/home/belle2/wangz/Work/SpinAlignment/offline/draw/belle_kid_eff/")
from get_pid_eff_weight import get_weight

    
def Ks_reco(df:R.RDataFrame, hasKs:bool=False, truth_match:bool=False):
    """
    inclusive Ks, rho00 measurement, reconstruction level processing
    """
    Ks_vars = []
    Ks_vars_toSave = []
    pion_vars = ["px_cms", "py_cms", "pz_cms", "E_cms"]
    pion_vars_toSave = []

    isMC = "pip_isSignal" in df.GetColumnNames()
    tools = RDF_process()

    # --- MC-only pre-selection ---
    if isMC:
        pion_vars += ["isSignal"]
        pion_vars_toSave += ["isSignal"]
        if truth_match:
            for pion in ["pip", "pim"]:
                for var in pion_vars:
                    df = df.Redefine(f"{pion}_{var}", f"{pion}_{var}[{pion}_isSignal]")
        elif hasKs:
            df = df.Filter("n_Ks_truth > 0", "has truth Ks")
        df = get_weight(df, period = "svd2", is_pion=True, cut_value=6, particle_names = ("pip", "pim"))
        df = df.Define("Ks_weight", "pip_pid_weight * pim_pid_weight")
        pion_vars += ["pid_weight"]
        pion_vars_toSave += ["pid_weight"]        
        Ks_vars += ["Ks_weight"]


    # --- test here
    """
    test = False
    if test:
        condition = "pip_p > 0.5 && pim_p > 0.5"
        for var in ["px_cms", "py_cms", "pz_cms", "isSignal"]:
            df = df.Redefine(f"pip_{var}", f"pip_{var}[{condition}]")
            df = df.Redefine(f"pim_{var}", f"pim_{var}[{condition}]")
    """

    # --- Build all pi+pi- pairs ---
    df = tools.save_all_pairs(df, 
                         p1_branches=("pip_px_cms","pip_py_cms","pip_pz_cms"),
                         p2_branches=("pim_px_cms","pim_py_cms","pim_pz_cms"),
                         mass=(0.13957061, 0.13957061),
                         particle_name=("Ks","pip","pim"),
                         cross_mode= False)
    Ks_vars += ["Ks_E", "Ks_px", "Ks_py", "Ks_pz", "pip_index", "pim_index", "Ks_helicity_angle", "Ks_helicity_phi"]
    Ks_vars_toSave += ["pip_index", "pim_index", "Ks_helicity_angle", "Ks_helicity_phi"]

    ## Fold azimuthal angle by quadrant: result in [0, pi/2]
    df = df.Redefine("Ks_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(Ks_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")


    # --- calculate other vars ---
    df = df.Define("Ks_M", "sqrt(Ks_E*Ks_E - Ks_px*Ks_px - Ks_py*Ks_py - Ks_pz*Ks_pz)")
    df = df.Define("Ks_z", "2*Ks_E/sqrts")
    df = tools.calculate_pt_toAxis(df, 
                                  particle=("Ks_px","Ks_py","Ks_pz"),
                                  axis = ("TMath::ACos(thrust[1])", "thrust[2]"),
                                  particle_name="Ks",axis_name="thrust")
    Ks_vars += ["Ks_M", "Ks_z", "Ks_thrust_pt", "Ks_thrust_costheta"]

    df = tools.convert_cartesian_to_spherical(df, particles=["Ks"])
    df = tools.convert_cartesian_to_spherical(df, particles=["pip", "pim"],
                                               px_branch = "px_cms", py_branch = "py_cms", pz_branch = "pz_cms")
    Ks_vars += ["Ks_costheta", "Ks_phi", "Ks_p"]
    Ks_vars_toSave += ["Ks_costheta", "Ks_phi", "Ks_p"]
    pion_vars += ["costheta", "phi", "p"]
    pion_vars_toSave += ["costheta", "phi", "p"]


    # --- Branches to save ---
    Br2Save = Ks_vars_toSave + ["pip_" + var for var in pion_vars_toSave] + ["pim_" + var for var in pion_vars_toSave]


    # --- Apply cuts ---
    df = df.Filter("Ks_M.size() > 0", "has Ks")
    df.Report().Print()

    return (df, Br2Save)

def Ks_reco_truth_matched(df:R.RDataFrame):
    """
    the situation is not same as phi_reco_truth_matched
    """
    tools = RDF_process()
    df = df.Filter("n_ks_truth > 0", "generate level has Ks")

    # ---------- reconstruct level -----------
    df = tools.save_all_pairs(df,
                            p1_branches=("pip_px_cms","pip_py_cms","pip_pz_cms"),
                            p2_branches=("pim_px_cms","pim_py_cms","pim_pz_cms"),
                            mass=(0.13957061, 0.13957061),
                            particle_name=("Ks","pip","pim"),
                            cross_mode= False)
    df = get_weight(df, period = "svd2", is_pion=True, cut_value=6, particle_names = ("pip", "pim"), eff_source="data")
    df = df.Define("Ks_eff", "pip_pid_eff * pim_pid_eff")
    df = tools.convert_cartesian_to_spherical(df, particles=["Ks"])
    df = tools.convert_cartesian_to_spherical(df, particles=["pip", "pim"],
                                               px_branch = "px_cms", py_branch = "py_cms", pz_branch = "pz_cms")
    
    df = df.Define("Ks_M", "sqrt(Ks_E*Ks_E - Ks_px*Ks_px - Ks_py*Ks_py - Ks_pz*Ks_pz)")
    df = df.Define("Ks_z", "2*Ks_E/sqrts")

    isSignal_condition = "pip_isSignal && pim_isSignal"
    for var in ["M", "z", "p", "costheta", "phi", "helicity_angle", "weight"]:
        df = df.Redefine(f"Ks_{var}", f"Ks_{var}[{isSignal_condition}]")
    for var in ["p", "costheta", "phi", "pid_eff"]:
        df = df.Redefine(f"pip_{var}", f"pip_{var}[{isSignal_condition}]")
        df = df.Redefine(f"pim_{var}", f"pim_{var}[{isSignal_condition}]")


    df = df.Filter("Ks_M.size() > 0", "has Ks after truth match")

    # ----------- generate level ---
    df = tools.save_all_pairs(df,
                              p1_branches=("pip_px_cms_gen","pip_py_cms_gen","pip_pz_cms_gen"),
                              p2_branches=("pim_px_cms_gen","pim_py_cms_gen","pim_pz_cms_gen"),
                              mass=(0.13957061, 0.13957061),
                              particle_name=("Ks_gen","pip_gen","pim_gen"),
                              cross_mode= False)
    df = df.Define("Ks_gen_M", "sqrt(Ks_gen_E*Ks_gen_E - Ks_gen_px*Ks_gen_px - Ks_gen_py*Ks_gen_py - Ks_gen_pz*Ks_gen_pz)")
    df = df.Define("Ks_gen_z", "2*Ks_gen_E/sqrts")


    Br2Save = ["Ks_M", "Ks_gen_M", "Ks_z", "Ks_gen_z", "Ks_helicity_angle", "Ks_gen_helicity_angle"]
    Br2Save += ["Ks_costheta", "Ks_phi", "Ks_p", "pip_p", "pip_costheta", "pip_phi", "pim_p", "pim_costheta", "pim_phi"]
    Br2Save += ["Ks_eff", "pip_pid_eff", "pim_pid_eff"]

    df.Report().Print()

    return (df, Br2Save)

def Ks_truth(df:R.RDataFrame):
    """
    generate level Ks, rho00 measurement, truth level processing
    """
    tools = RDF_process()
    sqrts = 10.516469955444336 
    df = tools.save_all_pairs(df,
                              p1_branches=("pip_px_cms_truth","pip_py_cms_truth","pip_pz_cms_truth"),
                              p2_branches=("pim_px_cms_truth","pim_py_cms_truth","pim_pz_cms_truth"),
                              mass=(0.13957061, 0.13957061),
                              particle_name=("Ks","pip","pim"),
                              cross_mode= False)
    df = df.Define("Ks_M", "sqrt(Ks_E*Ks_E - Ks_px*Ks_px - Ks_py*Ks_py - Ks_pz*Ks_pz)")
    df = df.Define("Ks_z", f"2*Ks_E/{sqrts}")
    df = df.Define("Ks_p", "sqrt(Ks_px*Ks_px + Ks_py*Ks_py + Ks_pz*Ks_pz)")
    
    Br2Save = ["Ks_M", "Ks_z", "Ks_helicity_angle", "Ks_p"]

    df.Report().Print()

    return (df, Br2Save)



def main():
    parser = argparse.ArgumentParser(description="set IO")
    parser.add_argument("which_channel",
                        choices=["Ks_reco", "Ks_reco_truth_matched", "Ks_truth"],
                        help="Channel to process")
    parser.add_argument("--input","-i", help="Input file path")
    parser.add_argument("--out","-o",help="Output file path, default: input_dir/processed_basename")
    parser.add_argument("--tree", "-t", default="event", help="Name of the TTree to process (default: event)")

    args = parser.parse_args()
    channel = args.which_channel 
    input_path = args.input
    if args.out:
        output_path = args.out
    else:
        input_dir = os.path.dirname(input_path)
        output_filename = f"processed_{os.path.basename(input_path)}"
        output_path = os.path.join(input_dir, output_filename)

    R.EnableImplicitMT(8)

    df  = R.RDataFrame(args.tree, input_path)
    result =globals()[channel](df) 
    if isinstance(result, tuple):
        df, branches_to_save = result
        print(f"Saving {len(branches_to_save)} branches: {branches_to_save}")
        df.Snapshot(args.tree , output_path, branches_to_save)
    else:
        df = result
        print("Saving all branches")
        df.Snapshot(args.tree, output_path)

if __name__ == "__main__":
    main()