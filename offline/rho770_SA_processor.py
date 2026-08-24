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
    
def reco(df:R.RDataFrame, hasrho:bool=False, truth_match:bool=False, cut_meson_thrust_costheta:bool=False):
    """
    inclusive rho(770), rho00 measurement, reconstruction level processing
    """
    rho_vars = []
    rho_vars_toSave = []
    pion_vars = ["px_cms", "py_cms", "pz_cms", "E_cms", "p", "costheta"]
    pion_vars_toSave = []

    is_MC = "pip_isSignal" in df.GetColumnNames()
    tools = RDF_process()

    # --- MC-only pre-selection ---
    if is_MC:
        pion_vars += ["isSignal"]
        pion_vars_toSave += ["isSignal"]
        if truth_match:
            for pion in ["pip", "pim"]:
                for var in pion_vars:
                    df = df.Redefine(f"{pion}_{var}", f"{pion}_{var}[{pion}_isSignal]")
        elif hasrho:  # require at least one truth-level rho
            df = df.Filter("n_rho_truth > 0", "has truth rho")
        df = get_weight(df, period="svd2", is_pion=True, cut_value=6, particle_names=("pip", "pim"))
        pion_vars += ["pid_weight"]
        pion_vars_toSave += ["pid_weight"]
  

    # --- pion momentum pre-selection ---
    for pion in ["pip", "pim"]:
        df = df.Define(f"{pion}_p_cut", f"{pion}_p > 0.5") # 0.5 GeV/c cut
        for var in pion_vars:
            df = df.Redefine(f"{pion}_{var}", f"{pion}_{var}[{pion}_p_cut]")

    df = df.Filter("pip_p.size() > 0 && pim_p.size() > 0", "has at least one pion pair after momentum cut")


    # --- Build all pi+pi- pairs ---
    df = tools.save_all_pairs(
        df,
        p1_branches=("pip_px_cms", "pip_py_cms", "pip_pz_cms"),
        p2_branches=("pim_px_cms", "pim_py_cms", "pim_pz_cms"),
        mass=(0.139570, 0.139570),
        particle_name=("rho", "pip", "pim"),
    )
    rho_vars += ["rho_E", "rho_px", "rho_py", "rho_pz", "pip_index", "pim_index", "rho_helicity_angle", "rho_helicity_phi"]
    rho_vars_toSave += ["pip_index", "pim_index", "rho_helicity_angle", "rho_helicity_phi"]

    ## Fold azimuthal angle by quadrant: result in [0, pi/2]
    df = df.Redefine("rho_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(rho_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")


    # --- Apply rho mass window (rho threshold ~0.770 GeV) --- cut 1
    df = df.Define("rho_M", "sqrt(rho_E*rho_E - rho_px*rho_px - rho_py*rho_py - rho_pz*rho_pz)")
    rho_vars.append("rho_M")
    rho_vars_toSave.append("rho_M")

    from DRAW import style_draw
    hist = df.Histo1D(("", ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});[MeV]", 100, 0.5, 1), "rho_M")
    style_draw([hist], "rho_mass.png")

    rho_mass_window = (0.770 - 0.150, 0.770 + 0.150)
    df = df.Define("mass_window_pass",
                   f"(rho_M > {rho_mass_window[0]}) && (rho_M < {rho_mass_window[1]})")
    for var in rho_vars:
        df = df.Redefine(var, f"{var}[mass_window_pass]")

    df = df.Filter("rho_M.size() > 0", "has rho after mass window cut")

    # --- Calculate angles related to thrust axis --- theta < 90^\circ cut 2 
    if cut_meson_thrust_costheta:
        df = tools.calculate_pt_toAxis(
            df,
            particle=("rho_px", "rho_py", "rho_pz"),
            axis=("TMath::ACos(thrust[1])", "thrust[2]"),
            particle_name="rho", axis_name="thrust",
        )
        rho_vars += ["rho_thrust_pt", "rho_thrust_costheta"]
        rho_vars_toSave += ["rho_thrust_pt", "rho_thrust_costheta"]

        df = df.Define("rho_thrust_costheta_pass", "rho_thrust_costheta > 0")
        for var in rho_vars:
            df = df.Redefine(var, f"{var}[rho_thrust_costheta_pass]")

        df = df.Filter("rho_M.size() > 0", "has rho after thrust costheta cut")


    # --- calculate other variables

    ## --- thrust , polarization angles --- 
    df = tools.cal_pola_angles(df, ("pip_px_cms", "pip_py_cms", "pip_pz_cms"),
                            ("pim_px_cms", "pim_py_cms", "pim_pz_cms"),
                            (0.139570, 0.139570), ("thrust[1]", "thrust[2]"), ("pip_index", "pim_index"),
                            ("rho_thrust_helicity_costheta", "rho_thrust_helicity_phi"))
    rho_vars += ["rho_thrust_helicity_costheta", "rho_thrust_helicity_phi"]
    rho_vars_toSave += ["rho_thrust_helicity_costheta", "rho_thrust_helicity_phi"]

    df = df.Redefine("rho_thrust_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(rho_thrust_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")
    
    ## --- MC-only: PID efficiency weight per pair ---
    if is_MC:
        df = df.Define(
            "rho_weight",
            "ROOT::VecOps::Take(pip_pid_weight, pip_index) * "
            "ROOT::VecOps::Take(pim_pid_weight, pim_index)",
        )
        rho_vars.append("rho_weight")
        rho_vars_toSave.append("rho_weight")

    ## --- other Kinematic variables ---
    df = tools.convert_cartesian_to_spherical(df, particles=["rho"])
    df = tools.convert_cartesian_to_spherical(df, particles=["pip", "pim"],
                                               px_branch="px_cms", py_branch="py_cms", pz_branch="pz_cms")
    df = df.Define("rho_z",  "2*rho_E/sqrts")
    rho_vars += ["rho_costheta", "rho_phi", "rho_p", "rho_z"]
    rho_vars_toSave += ["rho_costheta", "rho_phi", "rho_p", "rho_z"]
    pion_vars += ["p", "costheta", "phi"]
    pion_vars_toSave += ["p", "costheta", "phi"]

    # --- Branches to save ---
    Br2Save = rho_vars_toSave + \
        ["pip_" + var for var in pion_vars_toSave] + ["pim_" + var for var in pion_vars_toSave] + \
        ["thrust"]

    df.Report().Print()

    return (df, Br2Save)



def truth(df:R.RDataFrame):
    """
    process truth tree, retrieve coordinate generate leavel rho's z, cos*theta ...
    """
    rho_vars = []
    rho_vars_toSave = []
    pion_vars = ["px", "py", "pz", "E"]
    pion_vars = ["pip_" + var + "_cms_truth" for var in pion_vars] + ["pim_" + var + "_cms_truth" for var in pion_vars]

    tools = RDF_process()

    
    # --- Build all pi+pi- pairs and compute invariant mass ---
    df = tools.save_all_pairs(df, 
                         p1_branches=("pip_px_cms_truth","pip_py_cms_truth","pip_pz_cms_truth"),
                         p2_branches=("pim_px_cms_truth","pim_py_cms_truth","pim_pz_cms_truth"),
                         mass=(0.139570, 0.139570),
                         particle_name=("rho","pip","pim"),
                         cross_mode= False)
    rho_vars += ["rho_E", "rho_px", "rho_py", "rho_pz", "pip_index", "pim_index", "rho_helicity_angle", "rho_helicity_phi"]
    rho_vars_toSave += ["pip_index", "pim_index", "rho_helicity_angle", "rho_helicity_phi"]

    ## fold by quadrant : pi/2 - |fmod(rho, pi) - pi/2|, result in [0, pi/2]
    df = df.Redefine("rho_helicity_phi",
                   "TMath::Pi()/2 - abs(fmod(rho_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")
    
    # --- calculate thrust pola angles ---
    df = tools.cal_pola_angles(df, ("pip_px_cms_truth", "pip_py_cms_truth", "pip_pz_cms_truth"),
                            ("pim_px_cms_truth", "pim_py_cms_truth", "pim_pz_cms_truth"),
                            (0.139570, 0.139570), ("thrust_truth[1]", "thrust_truth[2]"), ("pip_index", "pim_index"),
                            ("rho_thrust_helicity_costheta", "rho_thrust_helicity_phi"))
    rho_vars += ["rho_thrust_helicity_costheta", "rho_thrust_helicity_phi"]
    rho_vars_toSave += ["rho_thrust_helicity_costheta", "rho_thrust_helicity_phi"]


    # --- calculate other variables
    df = df.Define("rho_M", "sqrt(rho_E*rho_E - rho_px*rho_px - rho_py*rho_py - rho_pz*rho_pz)")
    df = df.Define("rho_z", f"2*rho_E/10.516469955444336")
    rho_vars += ["rho_M", "rho_z"]
    rho_vars_toSave += ["rho_M", "rho_z"]
    

    # --- Branches to save ---
    Br2Save = rho_vars_toSave + pion_vars
    df.Report().Print()

    return (df, Br2Save)


# unclear truth match ,take care when using
def reco_truth_match(df:R.RDataFrame, phi_mass_window:Tuple[float,float]=(0.99, 1.15), cut_phi_thrust_costheta:bool=True):
    """
    using truth matched generic MC, for z, cos*theta, thrust, pt  resolution checking 
    save truth level info ,for resolution check
    """
    phi_vars = []
    phi_vars_toSave = []
    kaon_vars = ["px_cms", "py_cms", "pz_cms", "E_cms", "isSignal"]

    tools = RDF_process()

    # --- pre-selection: ---
    ## require at least one truth-level phi
    df = df.Filter("n_phi_truth > 0", "truth level has phi")

    ## truth match
    for var in kaon_vars:
        for kaon in ["kp", "km"]:
            df = df.Redefine(f"{kaon}_{var}", f"{kaon}_{var}[{kaon}_isSignal]")

    # --- Build all K+K- pairs ---
    df = tools.save_all_pairs(df, 
                         p1_branches=("kp_px_cms","kp_py_cms","kp_pz_cms"),
                         p2_branches=("km_px_cms","km_py_cms","km_pz_cms"),
                         mass=(0.493677, 0.493677),
                         particle_name=("phi","kp","km"))
    phi_vars += ["phi_E", "phi_px", "phi_py", "phi_pz", "kp_index", "km_index", "phi_helicity_angle", "phi_helicity_phi"]
    phi_vars_toSave += ["kp_index", "km_index", "phi_helicity_angle", "phi_helicity_phi"]

    ## Fold azimuthal angle by quadrant: result in [0, pi/2]
    df = df.Redefine("phi_helicity_phi",
                "TMath::Pi()/2 - abs(fmod(phi_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")


    # --- Apply phi_mass_window cut ---
    df = df.Define("phi_M", "sqrt(phi_E*phi_E - phi_px*phi_px - phi_py*phi_py - phi_pz*phi_pz)")
    phi_vars.append("phi_M")
    
    df = df.Define("mass_window_pass",
                   f"(phi_M > {phi_mass_window[0]}) && (phi_M < {phi_mass_window[1]})")
    for var in phi_vars:
        df = df.Redefine(var, f"{var}[mass_window_pass]")

    df = df.Filter("phi_M.size() > 0", "has phi after mass window cut")


    # --- Calculate angles related to thrust axis --- theta < 90^\circ cut 2 
    if cut_phi_thrust_costheta:
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
    

    # --- thrust , polarization angles --- 
    df = tools.cal_pola_angles(df, ("kp_px_cms", "kp_py_cms", "kp_pz_cms"),
                            ("km_px_cms", "km_py_cms", "km_pz_cms"),
                            (0.493677, 0.493677), ("thrust[1]", "thrust[2]"), ("kp_index", "km_index"),
                            ("phi_thrust_helicity_costheta", "phi_thrust_helicity_phi"))
    phi_vars += ["phi_thrust_helicity_costheta", "phi_thrust_helicity_phi"]
    phi_vars_toSave += ["phi_thrust_helicity_costheta", "phi_thrust_helicity_phi"]

    df = df.Redefine("phi_thrust_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(phi_thrust_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")
                     

    # --- calculate other variables
    df = df.Define("phi_z", "2*phi_E/sqrts")
    df = df.Redefine("thrust", "thrust[0]")
    phi_vars.append("phi_z")


    # ----------- generation info --------------------
    # --- Build all K+K- pairs ---
    df = tools.save_all_pairs(df, 
                        p1_branches=("kp_px_cms_gen","kp_py_cms_gen","kp_pz_cms_gen"),
                        p2_branches=("km_px_cms_gen","km_py_cms_gen","km_pz_cms_gen"),
                        mass=(0.493677, 0.493677),
                        particle_name=("phi_gen","kp_gen","km_gen"),
                        cross_mode= False) 
    phi_vars += ["phi_gen_E", "phi_gen_px", "phi_gen_py", "phi_gen_pz", "kp_gen_index", "km_gen_index", "phi_gen_helicity_angle", "phi_gen_helicity_phi"]
    phi_vars_toSave += ["kp_gen_index", "km_gen_index", "phi_gen_helicity_angle", "phi_gen_helicity_phi"]
    
    ## --- calculate thrust related variables ---
    df = tools.calculate_pt_toAxis(df, 
                        particle=("phi_gen_px","phi_gen_py","phi_gen_pz"),
                        axis = ("TMath::ACos(thrust_truth[1])", "thrust_truth[2]"),
                        particle_name="phi_gen",axis_name="thrust")
    phi_vars += ["phi_gen_thrust_pt", "phi_gen_thrust_costheta"]
    phi_vars_toSave += ["phi_gen_thrust_pt", "phi_gen_thrust_costheta"]
    
    df = df.Redefine("phi_gen_helicity_phi",
                "TMath::Pi()/2 - abs(fmod(phi_gen_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")
    
    ## --- thrust , polarization angles --- 
    df = tools.cal_pola_angles(df, ("kp_px_cms_gen", "kp_py_cms_gen", "kp_pz_cms_gen"),
                            ("km_px_cms_gen", "km_py_cms_gen", "km_pz_cms_gen"),
                            (0.493677, 0.493677), ("thrust_truth[1]", "thrust_truth[2]"), ("kp_index", "km_index"),
                            ("phi_gen_thrust_helicity_costheta", "phi_gen_thrust_helicity_phi"))
    phi_vars += ["phi_gen_thrust_helicity_costheta", "phi_gen_thrust_helicity_phi"]
    phi_vars_toSave += ["phi_gen_thrust_helicity_costheta", "phi_gen_thrust_helicity_phi"]

    df = df.Redefine("phi_gen_thrust_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(phi_gen_thrust_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")
    
    ## --- calculate other variables ---
    df = df.Redefine("thrust_truth", "thrust_truth[0]")
    df = df.Define("phi_gen_z", f"2*phi_gen_E/sqrts")
    df = df.Define("phi_gen_M", "sqrt(phi_gen_E*phi_gen_E - phi_gen_px*phi_gen_px - phi_gen_py*phi_gen_py - phi_gen_pz*phi_gen_pz)")
    phi_vars += ["phi_gen_M", "phi_gen_z"]
    phi_vars_toSave += ["phi_gen_M", "phi_gen_z"]


    # ----
    # kaon isSignal but there are some misidentified phi , suppressed by mass window and n phi cut , temporarily not care the still exist misidentified phi
    df = df.Filter("phi_M.size() == phi_gen_M.size()", "match reconstructed and generated phi count") 
    df = df.Filter("phi_M.size() > 0", "has phi after match")

    df.Report().Print()


    Br2Save = phi_vars_toSave + ["thrust", "thrust_truth"] 

    return (df, Br2Save)


def main():
    parser = argparse.ArgumentParser(description="set IO")
    parser.add_argument("which_channel",
                        choices=["reco", "truth", "reco_truth_match"],
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

    df  = R.RDataFrame(args.tree, input_path)
    result =globals()[channel](df) 
    if isinstance(result, tuple):
        df, branches_to_save = result
        existing_cols = set(df.GetColumnNames())
        branches_to_save = [b for b in branches_to_save if b in existing_cols]
        print(f"Saving {len(branches_to_save)} branches: {branches_to_save}")
        df.Snapshot(args.tree , output_path, branches_to_save)
    else:
        df = result
        print("Saving all branches")
        df.Snapshot(args.tree, output_path)

if __name__ == "__main__":
    main()
    #batch_process()