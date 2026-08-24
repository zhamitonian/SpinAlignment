#!/usr/bin/env python3
# K*(892)0 → K+π-   (kstar)
# K̄*(892)0 → K-π+   (akstar)
# PDG: M = 895.55 MeV, Γ = 47.3 MeV
# kaon mass = 493.677 MeV,  pion mass = 139.570 MeV

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

KAON_MASS   = 0.493677  # GeV
PION_MASS   = 0.139570  # GeV
KSTAR_MASS_WINDOW = (0.70, 1.1)  # ±3Γ


def reco(df: R.RDataFrame, hasKstar: bool = False, truth_match: bool = False,
         cut_kstar_thrust_costheta: bool = False):
    """
    Inclusive K*(892)0 / K̄*(892)0 reconstruction level processing.
    K*(892)0  → K⁺π⁻   (kstar)
    K̄*(892)0 → K⁻π⁺   (akstar)
    """
    # independent variable lists for each charge state
    kstar_vars       = []
    kstar_vars_toSave = []
    akstar_vars       = []
    akstar_vars_toSave = []

    kaon_vars        = ["px_cms", "py_cms", "pz_cms", "E_cms", "p", "costheta"]
    pion_vars        = ["px_cms", "py_cms", "pz_cms", "E_cms", "p", "costheta"]
    kaon_vars_toSave = []
    pion_vars_toSave = []

    is_MC = "kp_isSignal" in df.GetColumnNames()
    tools = RDF_process()

    # --- MC-only pre-selection ---
    if is_MC:
        kaon_vars += ["isSignal"]
        kaon_vars_toSave += ["isSignal"]
        pion_vars += ["isSignal"]
        pion_vars_toSave += ["isSignal"]
        if truth_match:
            for kaon in ["kp", "km"]:
                for var in kaon_vars:
                    df = df.Redefine(f"{kaon}_{var}", f"{kaon}_{var}[{kaon}_isSignal]")
            for pion in ["pip", "pim"]:
                for var in pion_vars:
                    df = df.Redefine(f"{pion}_{var}", f"{pion}_{var}[{pion}_isSignal]")
        elif hasKstar:
            df = df.Filter("n_kstar_truth > 0", "has truth Kstar")
        df = get_weight(df, period="svd2", is_pion=False, cut_value=6, particle_names=("kp", "km"))
        df = get_weight(df, period="svd2", is_pion=True,  cut_value=6, particle_names=("pip", "pim"))
        kaon_vars += ["pid_weight"]
        kaon_vars_toSave += ["pid_weight"]
        pion_vars += ["pid_weight"]
        pion_vars_toSave += ["pid_weight"]


    # -----------------------------------------------------------------------
    # K*(892)0: K⁺ π⁻
    # -----------------------------------------------------------------------
    df = tools.save_all_pairs(
        df,
        p1_branches=("kp_px_cms",  "kp_py_cms",  "kp_pz_cms"),
        p2_branches=("pim_px_cms", "pim_py_cms", "pim_pz_cms"),
        mass=(KAON_MASS, PION_MASS),
        particle_name=("kstar", "kp", "pim"),
    )
    kstar_vars += ["kstar_E", "kstar_px", "kstar_py", "kstar_pz",
                   "kp_index", "pim_index", "kstar_helicity_angle", "kstar_helicity_phi"]
    kstar_vars_toSave += ["kp_index", "pim_index", "kstar_helicity_angle", "kstar_helicity_phi"]

    ## Fold azimuthal angle by quadrant: result in [0, π/2]
    df = df.Redefine("kstar_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(kstar_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    ## Apply K*0 mass window --- cut 1
    df = df.Define("kstar_M",
                   "sqrt(kstar_E*kstar_E - kstar_px*kstar_px - kstar_py*kstar_py - kstar_pz*kstar_pz)")
    kstar_vars.append("kstar_M")
    kstar_vars_toSave.append("kstar_M")

    df = df.Define("kstar_mass_window_pass",
                   f"(kstar_M > {KSTAR_MASS_WINDOW[0]}) && (kstar_M < {KSTAR_MASS_WINDOW[1]})")
    for var in kstar_vars:
        df = df.Redefine(var, f"{var}[kstar_mass_window_pass]")

    df = df.Filter("kstar_M.size() > 0", "has K*0 after mass window cut")

    ## thrust costheta cut --- cut 2
    if cut_kstar_thrust_costheta:
        df = tools.calculate_pt_toAxis(
            df,
            particle=("kstar_px", "kstar_py", "kstar_pz"),
            axis=("TMath::ACos(thrust[1])", "thrust[2]"),
            particle_name="kstar", axis_name="thrust",
        )
        kstar_vars += ["kstar_thrust_pt", "kstar_thrust_costheta"]
        kstar_vars_toSave += ["kstar_thrust_pt", "kstar_thrust_costheta"]

        df = df.Define("kstar_thrust_costheta_pass", "kstar_thrust_costheta > 0")
        for var in kstar_vars:
            df = df.Redefine(var, f"{var}[kstar_thrust_costheta_pass]")
        df = df.Filter("kstar_M.size() > 0", "has K*0 after thrust costheta cut")

    ## thrust-axis polarization angles
    df = tools.cal_pola_angles(df,
                               ("kp_px_cms",  "kp_py_cms",  "kp_pz_cms"),
                               ("pim_px_cms", "pim_py_cms", "pim_pz_cms"),
                               (KAON_MASS, PION_MASS),
                               ("thrust[1]", "thrust[2]"),
                               ("kp_index", "pim_index"),
                               ("kstar_thrust_helicity_costheta", "kstar_thrust_helicity_phi"))
    kstar_vars += ["kstar_thrust_helicity_costheta", "kstar_thrust_helicity_phi"]
    kstar_vars_toSave += ["kstar_thrust_helicity_costheta", "kstar_thrust_helicity_phi"]

    df = df.Redefine("kstar_thrust_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(kstar_thrust_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    ## MC-only: per-pair PID efficiency weight
    if is_MC:
        df = df.Define(
            "kstar_weight",
            "ROOT::VecOps::Take(kp_pid_weight,  kp_index)  * "
            "ROOT::VecOps::Take(pim_pid_weight, pim_index)",
        )
        kstar_vars.append("kstar_weight")
        kstar_vars_toSave.append("kstar_weight")

    ## other kinematic variables
    df = tools.convert_cartesian_to_spherical(df, particles=["kstar"])
    df = df.Define("kstar_z", "2*kstar_E/sqrts")
    kstar_vars += ["kstar_costheta", "kstar_phi", "kstar_p", "kstar_z"]
    kstar_vars_toSave += ["kstar_costheta", "kstar_phi", "kstar_p", "kstar_z"]


    # -----------------------------------------------------------------------
    # K̄*(892)0: K⁻ π⁺
    # -----------------------------------------------------------------------
    df = tools.save_all_pairs(
        df,
        p1_branches=("km_px_cms",  "km_py_cms",  "km_pz_cms"),
        p2_branches=("pip_px_cms", "pip_py_cms", "pip_pz_cms"),
        mass=(KAON_MASS, PION_MASS),
        particle_name=("akstar", "km", "pip"),
    )
    akstar_vars += ["akstar_E", "akstar_px", "akstar_py", "akstar_pz",
                    "km_index", "pip_index", "akstar_helicity_angle", "akstar_helicity_phi"]
    akstar_vars_toSave += ["km_index", "pip_index", "akstar_helicity_angle", "akstar_helicity_phi"]

    df = df.Redefine("akstar_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(akstar_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    ## Apply K̄*0 mass window --- cut 1
    df = df.Define("akstar_M",
                   "sqrt(akstar_E*akstar_E - akstar_px*akstar_px - akstar_py*akstar_py - akstar_pz*akstar_pz)")
    akstar_vars.append("akstar_M")
    akstar_vars_toSave.append("akstar_M")

    df = df.Define("akstar_mass_window_pass",
                   f"(akstar_M > {KSTAR_MASS_WINDOW[0]}) && (akstar_M < {KSTAR_MASS_WINDOW[1]})")
    for var in akstar_vars:
        df = df.Redefine(var, f"{var}[akstar_mass_window_pass]")

    df = df.Filter("akstar_M.size() > 0", "has K̄*0 after mass window cut")

    ## thrust costheta cut --- cut 2
    if cut_kstar_thrust_costheta:
        df = tools.calculate_pt_toAxis(
            df,
            particle=("akstar_px", "akstar_py", "akstar_pz"),
            axis=("TMath::ACos(thrust[1])", "thrust[2]"),
            particle_name="akstar", axis_name="thrust",
        )
        akstar_vars += ["akstar_thrust_pt", "akstar_thrust_costheta"]
        akstar_vars_toSave += ["akstar_thrust_pt", "akstar_thrust_costheta"]

        df = df.Define("akstar_thrust_costheta_pass", "akstar_thrust_costheta > 0")
        for var in akstar_vars:
            df = df.Redefine(var, f"{var}[akstar_thrust_costheta_pass]")
        df = df.Filter("akstar_M.size() > 0", "has K̄*0 after thrust costheta cut")

    ## thrust-axis polarization angles
    df = tools.cal_pola_angles(df,
                               ("km_px_cms",  "km_py_cms",  "km_pz_cms"),
                               ("pip_px_cms", "pip_py_cms", "pip_pz_cms"),
                               (KAON_MASS, PION_MASS),
                               ("thrust[1]", "thrust[2]"),
                               ("km_index", "pip_index"),
                               ("akstar_thrust_helicity_costheta", "akstar_thrust_helicity_phi"))
    akstar_vars += ["akstar_thrust_helicity_costheta", "akstar_thrust_helicity_phi"]
    akstar_vars_toSave += ["akstar_thrust_helicity_costheta", "akstar_thrust_helicity_phi"]

    df = df.Redefine("akstar_thrust_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(akstar_thrust_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    ## MC-only: per-pair PID efficiency weight
    if is_MC:
        df = df.Define(
            "akstar_weight",
            "ROOT::VecOps::Take(km_pid_weight,  km_index)  * "
            "ROOT::VecOps::Take(pip_pid_weight, pip_index)",
        )
        akstar_vars.append("akstar_weight")
        akstar_vars_toSave.append("akstar_weight")

    ## other kinematic variables
    df = tools.convert_cartesian_to_spherical(df, particles=["akstar"])
    df = df.Define("akstar_z", "2*akstar_E/sqrts")
    akstar_vars += ["akstar_costheta", "akstar_phi", "akstar_p", "akstar_z"]
    akstar_vars_toSave += ["akstar_costheta", "akstar_phi", "akstar_p", "akstar_z"]

    ## daughter kinematics
    df = tools.convert_cartesian_to_spherical(df, particles=["kp", "km"],
                                               px_branch="px_cms", py_branch="py_cms", pz_branch="pz_cms")
    df = tools.convert_cartesian_to_spherical(df, particles=["pip", "pim"],
                                               px_branch="px_cms", py_branch="py_cms", pz_branch="pz_cms")
    kaon_vars += ["p", "costheta", "phi"]
    kaon_vars_toSave += ["p", "costheta", "phi"]
    pion_vars += ["p", "costheta", "phi"]
    pion_vars_toSave += ["p", "costheta", "phi"]

    # --- Combine K*0 and K̄*0 branches for fitting ---
    # Variables that exist for both charge states and should be merged
    combined_vars = ["M", "helicity_angle", "helicity_phi","z"]

    if cut_kstar_thrust_costheta:
        combined_vars += ["thrust_pt", "thrust_costheta"]
    if is_MC:
        combined_vars += ["weight"]

    for var in combined_vars:
        df = df.Define(f"Kstar_{var}",
                       f"ROOT::VecOps::Concatenate(kstar_{var}, akstar_{var})")

    combined_toSave = [f"Kstar_{var}" for var in combined_vars]

    # --- Branches to save ---
    Br2Save = kstar_vars_toSave + akstar_vars_toSave + combined_toSave + \
        ["kp_"  + var for var in kaon_vars_toSave] + \
        ["km_"  + var for var in kaon_vars_toSave] + \
        ["pip_" + var for var in pion_vars_toSave] + \
        ["pim_" + var for var in pion_vars_toSave] + \
        ["thrust"]

    df.Report().Print()

    return (df, Br2Save)


def truth(df: R.RDataFrame, pion_inBarrel: bool = False):
    """
    Process truth tree.
    K*(892)0  truth: kp + pim at generator level  (kstar)
    K̄*(892)0 truth: km + pip at generator level  (akstar)
    """
    kstar_vars        = []
    kstar_vars_toSave = []
    akstar_vars        = []
    akstar_vars_toSave = []

    truth_branches_kp  = ["kp_"  + v + "_cms_truth" for v in ["px", "py", "pz", "E"]]
    truth_branches_pim = ["pim_" + v + "_cms_truth" for v in ["px", "py", "pz", "E"]]
    truth_branches_km  = ["km_"  + v + "_cms_truth" for v in ["px", "py", "pz", "E"]]
    truth_branches_pip = ["pip_" + v + "_cms_truth" for v in ["px", "py", "pz", "E"]]

    tools = RDF_process()

    # -----------------------------------------------------------------------
    # K*(892)0: K⁺ π⁻  (truth)
    # -----------------------------------------------------------------------
    df = tools.save_all_pairs(df,
                              p1_branches=("kp_px_cms_truth",  "kp_py_cms_truth",  "kp_pz_cms_truth"),
                              p2_branches=("pim_px_cms_truth", "pim_py_cms_truth", "pim_pz_cms_truth"),
                              mass=(KAON_MASS, PION_MASS),
                              particle_name=("kstar", "kp", "pim"),
                              cross_mode=False)
    kstar_vars += ["kstar_E", "kstar_px", "kstar_py", "kstar_pz",
                   "kp_index", "pim_index", "kstar_helicity_angle", "kstar_helicity_phi"]
    kstar_vars_toSave += ["kp_index", "pim_index", "kstar_helicity_angle", "kstar_helicity_phi"]

    df = df.Redefine("kstar_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(kstar_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    ## thrust-axis polarization angles
    df = tools.cal_pola_angles(df,
                               ("kp_px_cms_truth",  "kp_py_cms_truth",  "kp_pz_cms_truth"),
                               ("pim_px_cms_truth", "pim_py_cms_truth", "pim_pz_cms_truth"),
                               (KAON_MASS, PION_MASS),
                               ("thrust_truth[1]", "thrust_truth[2]"),
                               ("kp_index", "pim_index"),
                               ("kstar_thrust_helicity_costheta", "kstar_thrust_helicity_phi"))
    kstar_vars += ["kstar_thrust_helicity_costheta", "kstar_thrust_helicity_phi"]
    kstar_vars_toSave += ["kstar_thrust_helicity_costheta", "kstar_thrust_helicity_phi"]

    df = df.Define("kstar_M",
                   "sqrt(kstar_E*kstar_E - kstar_px*kstar_px - kstar_py*kstar_py - kstar_pz*kstar_pz)")
    df = df.Define("kstar_z", f"2*kstar_E/10.516469955444336")
    kstar_vars += ["kstar_M", "kstar_z"]
    kstar_vars_toSave += ["kstar_M", "kstar_z"]

    # -----------------------------------------------------------------------
    # K̄*(892)0: K⁻ π⁺  (truth)
    # -----------------------------------------------------------------------
    df = tools.save_all_pairs(df,
                              p1_branches=("km_px_cms_truth",  "km_py_cms_truth",  "km_pz_cms_truth"),
                              p2_branches=("pip_px_cms_truth", "pip_py_cms_truth", "pip_pz_cms_truth"),
                              mass=(KAON_MASS, PION_MASS),
                              particle_name=("akstar", "km", "pip"),
                              cross_mode=False)
    akstar_vars += ["akstar_E", "akstar_px", "akstar_py", "akstar_pz",
                    "km_index", "pip_index", "akstar_helicity_angle", "akstar_helicity_phi"]
    akstar_vars_toSave += ["km_index", "pip_index", "akstar_helicity_angle", "akstar_helicity_phi"]

    df = df.Redefine("akstar_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(akstar_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    ## thrust-axis polarization angles
    df = tools.cal_pola_angles(df,
                               ("km_px_cms_truth",  "km_py_cms_truth",  "km_pz_cms_truth"),
                               ("pip_px_cms_truth", "pip_py_cms_truth", "pip_pz_cms_truth"),
                               (KAON_MASS, PION_MASS),
                               ("thrust_truth[1]", "thrust_truth[2]"),
                               ("km_index", "pip_index"),
                               ("akstar_thrust_helicity_costheta", "akstar_thrust_helicity_phi"))
    akstar_vars += ["akstar_thrust_helicity_costheta", "akstar_thrust_helicity_phi"]
    akstar_vars_toSave += ["akstar_thrust_helicity_costheta", "akstar_thrust_helicity_phi"]

    df = df.Define("akstar_M",
                   "sqrt(akstar_E*akstar_E - akstar_px*akstar_px - akstar_py*akstar_py - akstar_pz*akstar_pz)")
    df = df.Define("akstar_z", f"2*akstar_E/10.516469955444336")
    akstar_vars += ["akstar_M", "akstar_z"]
    akstar_vars_toSave += ["akstar_M", "akstar_z"]


    # Variables that exist for both charge states and should be merged
    combined_vars = ["M", "helicity_angle", "helicity_phi","z"]

    for var in combined_vars:
        df = df.Define(f"Kstar_{var}",
                       f"ROOT::VecOps::Concatenate(kstar_{var}, akstar_{var})")

    combined_toSave = [f"Kstar_{var}" for var in combined_vars]

    # --- Branches to save ---
    Br2Save = kstar_vars_toSave + akstar_vars_toSave + \
        truth_branches_kp + truth_branches_pim + \
        truth_branches_km + truth_branches_pip + \
        combined_toSave

    df.Report().Print()

    return (df, Br2Save)


# unclear truth match, take care when using
def reco_truth_match(df: R.RDataFrame,
                     kstar_mass_window: Tuple[float, float] = KSTAR_MASS_WINDOW,
                     cut_kstar_thrust_costheta: bool = True):
    """
    Truth-matched generic MC: resolution studies for z, cosθ, thrust, pt.
    Saves both reco-level and gen-level K*(892)0 / K̄*(892)0 info.
    """
    kstar_vars        = []
    kstar_vars_toSave = []
    akstar_vars        = []
    akstar_vars_toSave = []

    kaon_vars = ["px_cms", "py_cms", "pz_cms", "E_cms", "isSignal"]
    pion_vars = ["px_cms", "py_cms", "pz_cms", "E_cms", "isSignal"]

    tools = RDF_process()

    # --- pre-selection: require at least one truth-level K* ---
    df = df.Filter("n_kstar_truth > 0", "truth level has K*")

    ## truth match
    for kaon in ["kp", "km"]:
        for var in kaon_vars:
            df = df.Redefine(f"{kaon}_{var}", f"{kaon}_{var}[{kaon}_isSignal]")
    for pion in ["pip", "pim"]:
        for var in pion_vars:
            df = df.Redefine(f"{pion}_{var}", f"{pion}_{var}[{pion}_isSignal]")


    # -----------------------------------------------------------------------
    # RECO: K*(892)0 → K⁺ π⁻
    # -----------------------------------------------------------------------
    df = tools.save_all_pairs(df,
                              p1_branches=("kp_px_cms",  "kp_py_cms",  "kp_pz_cms"),
                              p2_branches=("pim_px_cms", "pim_py_cms", "pim_pz_cms"),
                              mass=(KAON_MASS, PION_MASS),
                              particle_name=("kstar", "kp", "pim"))
    kstar_vars += ["kstar_E", "kstar_px", "kstar_py", "kstar_pz",
                   "kp_index", "pim_index", "kstar_helicity_angle", "kstar_helicity_phi"]
    kstar_vars_toSave += ["kp_index", "pim_index", "kstar_helicity_angle", "kstar_helicity_phi"]

    df = df.Redefine("kstar_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(kstar_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    ## mass window
    df = df.Define("kstar_M",
                   "sqrt(kstar_E*kstar_E - kstar_px*kstar_px - kstar_py*kstar_py - kstar_pz*kstar_pz)")
    kstar_vars.append("kstar_M")

    df = df.Define("kstar_mass_window_pass",
                   f"(kstar_M > {kstar_mass_window[0]}) && (kstar_M < {kstar_mass_window[1]})")
    for var in kstar_vars:
        df = df.Redefine(var, f"{var}[kstar_mass_window_pass]")
    df = df.Filter("kstar_M.size() > 0", "has K*0 after mass window cut")

    ## thrust costheta cut
    if cut_kstar_thrust_costheta:
        df = tools.calculate_pt_toAxis(df,
                                       particle=("kstar_px", "kstar_py", "kstar_pz"),
                                       axis=("TMath::ACos(thrust[1])", "thrust[2]"),
                                       particle_name="kstar", axis_name="thrust")
        kstar_vars += ["kstar_thrust_pt", "kstar_thrust_costheta"]
        kstar_vars_toSave += ["kstar_thrust_pt", "kstar_thrust_costheta"]

        df = df.Define("kstar_thrust_costheta_pass", "kstar_thrust_costheta > 0")
        for var in kstar_vars:
            df = df.Redefine(var, f"{var}[kstar_thrust_costheta_pass]")
        df = df.Filter("kstar_M.size() > 0", "has K*0 after thrust costheta cut")

    ## thrust-axis polarization angles
    df = tools.cal_pola_angles(df,
                               ("kp_px_cms",  "kp_py_cms",  "kp_pz_cms"),
                               ("pim_px_cms", "pim_py_cms", "pim_pz_cms"),
                               (KAON_MASS, PION_MASS),
                               ("thrust[1]", "thrust[2]"),
                               ("kp_index", "pim_index"),
                               ("kstar_thrust_helicity_costheta", "kstar_thrust_helicity_phi"))
    kstar_vars += ["kstar_thrust_helicity_costheta", "kstar_thrust_helicity_phi"]
    kstar_vars_toSave += ["kstar_thrust_helicity_costheta", "kstar_thrust_helicity_phi"]

    df = df.Redefine("kstar_thrust_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(kstar_thrust_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    df = df.Define("kstar_z", "2*kstar_E/sqrts")
    kstar_vars.append("kstar_z")
    kstar_vars_toSave.append("kstar_M")
    kstar_vars_toSave.append("kstar_z")


    # -----------------------------------------------------------------------
    # RECO: K̄*(892)0 → K⁻ π⁺
    # -----------------------------------------------------------------------
    df = tools.save_all_pairs(df,
                              p1_branches=("km_px_cms",  "km_py_cms",  "km_pz_cms"),
                              p2_branches=("pip_px_cms", "pip_py_cms", "pip_pz_cms"),
                              mass=(KAON_MASS, PION_MASS),
                              particle_name=("akstar", "km", "pip"))
    akstar_vars += ["akstar_E", "akstar_px", "akstar_py", "akstar_pz",
                    "km_index", "pip_index", "akstar_helicity_angle", "akstar_helicity_phi"]
    akstar_vars_toSave += ["km_index", "pip_index", "akstar_helicity_angle", "akstar_helicity_phi"]

    df = df.Redefine("akstar_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(akstar_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    ## mass window
    df = df.Define("akstar_M",
                   "sqrt(akstar_E*akstar_E - akstar_px*akstar_px - akstar_py*akstar_py - akstar_pz*akstar_pz)")
    akstar_vars.append("akstar_M")

    df = df.Define("akstar_mass_window_pass",
                   f"(akstar_M > {kstar_mass_window[0]}) && (akstar_M < {kstar_mass_window[1]})")
    for var in akstar_vars:
        df = df.Redefine(var, f"{var}[akstar_mass_window_pass]")
    df = df.Filter("akstar_M.size() > 0", "has K̄*0 after mass window cut")

    ## thrust costheta cut
    if cut_kstar_thrust_costheta:
        df = tools.calculate_pt_toAxis(df,
                                       particle=("akstar_px", "akstar_py", "akstar_pz"),
                                       axis=("TMath::ACos(thrust[1])", "thrust[2]"),
                                       particle_name="akstar", axis_name="thrust")
        akstar_vars += ["akstar_thrust_pt", "akstar_thrust_costheta"]
        akstar_vars_toSave += ["akstar_thrust_pt", "akstar_thrust_costheta"]

        df = df.Define("akstar_thrust_costheta_pass", "akstar_thrust_costheta > 0")
        for var in akstar_vars:
            df = df.Redefine(var, f"{var}[akstar_thrust_costheta_pass]")
        df = df.Filter("akstar_M.size() > 0", "has K̄*0 after thrust costheta cut")

    ## thrust-axis polarization angles
    df = tools.cal_pola_angles(df,
                               ("km_px_cms",  "km_py_cms",  "km_pz_cms"),
                               ("pip_px_cms", "pip_py_cms", "pip_pz_cms"),
                               (KAON_MASS, PION_MASS),
                               ("thrust[1]", "thrust[2]"),
                               ("km_index", "pip_index"),
                               ("akstar_thrust_helicity_costheta", "akstar_thrust_helicity_phi"))
    akstar_vars += ["akstar_thrust_helicity_costheta", "akstar_thrust_helicity_phi"]
    akstar_vars_toSave += ["akstar_thrust_helicity_costheta", "akstar_thrust_helicity_phi"]

    df = df.Redefine("akstar_thrust_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(akstar_thrust_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    df = df.Define("akstar_z", "2*akstar_E/sqrts")
    akstar_vars.append("akstar_z")
    akstar_vars_toSave.append("akstar_M")
    akstar_vars_toSave.append("akstar_z")

    df = df.Redefine("thrust", "thrust[0]")


    # -----------------------------------------------------------------------
    # GEN: K*(892)0 → K⁺ π⁻  (generator level)
    # -----------------------------------------------------------------------
    df = tools.save_all_pairs(df,
                              p1_branches=("kp_px_cms_gen",  "kp_py_cms_gen",  "kp_pz_cms_gen"),
                              p2_branches=("pim_px_cms_gen", "pim_py_cms_gen", "pim_pz_cms_gen"),
                              mass=(KAON_MASS, PION_MASS),
                              particle_name=("kstar_gen", "kp_gen", "pim_gen"),
                              cross_mode=False)
    kstar_vars += ["kstar_gen_E", "kstar_gen_px", "kstar_gen_py", "kstar_gen_pz",
                   "kp_gen_index", "pim_gen_index", "kstar_gen_helicity_angle", "kstar_gen_helicity_phi"]
    kstar_vars_toSave += ["kp_gen_index", "pim_gen_index",
                          "kstar_gen_helicity_angle", "kstar_gen_helicity_phi"]

    df = tools.calculate_pt_toAxis(df,
                                   particle=("kstar_gen_px", "kstar_gen_py", "kstar_gen_pz"),
                                   axis=("TMath::ACos(thrust_truth[1])", "thrust_truth[2]"),
                                   particle_name="kstar_gen", axis_name="thrust")
    kstar_vars += ["kstar_gen_thrust_pt", "kstar_gen_thrust_costheta"]
    kstar_vars_toSave += ["kstar_gen_thrust_pt", "kstar_gen_thrust_costheta"]

    df = df.Redefine("kstar_gen_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(kstar_gen_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    df = tools.cal_pola_angles(df,
                               ("kp_px_cms_gen",  "kp_py_cms_gen",  "kp_pz_cms_gen"),
                               ("pim_px_cms_gen", "pim_py_cms_gen", "pim_pz_cms_gen"),
                               (KAON_MASS, PION_MASS),
                               ("thrust_truth[1]", "thrust_truth[2]"),
                               ("kp_gen_index", "pim_gen_index"),
                               ("kstar_gen_thrust_helicity_costheta", "kstar_gen_thrust_helicity_phi"))
    kstar_vars += ["kstar_gen_thrust_helicity_costheta", "kstar_gen_thrust_helicity_phi"]
    kstar_vars_toSave += ["kstar_gen_thrust_helicity_costheta", "kstar_gen_thrust_helicity_phi"]

    df = df.Redefine("kstar_gen_thrust_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(kstar_gen_thrust_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    df = df.Define("kstar_gen_M",
                   "sqrt(kstar_gen_E*kstar_gen_E - kstar_gen_px*kstar_gen_px - kstar_gen_py*kstar_gen_py - kstar_gen_pz*kstar_gen_pz)")
    df = df.Define("kstar_gen_z", "2*kstar_gen_E/sqrts")
    kstar_vars += ["kstar_gen_M", "kstar_gen_z"]
    kstar_vars_toSave += ["kstar_gen_M", "kstar_gen_z"]


    # -----------------------------------------------------------------------
    # GEN: K̄*(892)0 → K⁻ π⁺  (generator level)
    # -----------------------------------------------------------------------
    df = tools.save_all_pairs(df,
                              p1_branches=("km_px_cms_gen",  "km_py_cms_gen",  "km_pz_cms_gen"),
                              p2_branches=("pip_px_cms_gen", "pip_py_cms_gen", "pip_pz_cms_gen"),
                              mass=(KAON_MASS, PION_MASS),
                              particle_name=("akstar_gen", "km_gen", "pip_gen"),
                              cross_mode=False)
    akstar_vars += ["akstar_gen_E", "akstar_gen_px", "akstar_gen_py", "akstar_gen_pz",
                    "km_gen_index", "pip_gen_index", "akstar_gen_helicity_angle", "akstar_gen_helicity_phi"]
    akstar_vars_toSave += ["km_gen_index", "pip_gen_index",
                           "akstar_gen_helicity_angle", "akstar_gen_helicity_phi"]

    df = tools.calculate_pt_toAxis(df,
                                   particle=("akstar_gen_px", "akstar_gen_py", "akstar_gen_pz"),
                                   axis=("TMath::ACos(thrust_truth[1])", "thrust_truth[2]"),
                                   particle_name="akstar_gen", axis_name="thrust")
    akstar_vars += ["akstar_gen_thrust_pt", "akstar_gen_thrust_costheta"]
    akstar_vars_toSave += ["akstar_gen_thrust_pt", "akstar_gen_thrust_costheta"]

    df = df.Redefine("akstar_gen_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(akstar_gen_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    df = tools.cal_pola_angles(df,
                               ("km_px_cms_gen",  "km_py_cms_gen",  "km_pz_cms_gen"),
                               ("pip_px_cms_gen", "pip_py_cms_gen", "pip_pz_cms_gen"),
                               (KAON_MASS, PION_MASS),
                               ("thrust_truth[1]", "thrust_truth[2]"),
                               ("km_gen_index", "pip_gen_index"),
                               ("akstar_gen_thrust_helicity_costheta", "akstar_gen_thrust_helicity_phi"))
    akstar_vars += ["akstar_gen_thrust_helicity_costheta", "akstar_gen_thrust_helicity_phi"]
    akstar_vars_toSave += ["akstar_gen_thrust_helicity_costheta", "akstar_gen_thrust_helicity_phi"]

    df = df.Redefine("akstar_gen_thrust_helicity_phi",
                     "TMath::Pi()/2 - abs(fmod(akstar_gen_thrust_helicity_phi, TMath::Pi()) - TMath::Pi()/2)")

    df = df.Define("akstar_gen_M",
                   "sqrt(akstar_gen_E*akstar_gen_E - akstar_gen_px*akstar_gen_px - akstar_gen_py*akstar_gen_py - akstar_gen_pz*akstar_gen_pz)")
    df = df.Define("akstar_gen_z", "2*akstar_gen_E/sqrts")
    akstar_vars += ["akstar_gen_M", "akstar_gen_z"]
    akstar_vars_toSave += ["akstar_gen_M", "akstar_gen_z"]

    df = df.Redefine("thrust_truth", "thrust_truth[0]")

    # require matching reco/gen multiplicity (per charge state independently)
    df = df.Filter("kstar_M.size()  == kstar_gen_M.size()",  "K*0  reco/gen count match")
    df = df.Filter("akstar_M.size() == akstar_gen_M.size()", "K̄*0 reco/gen count match")
    df = df.Filter("kstar_M.size() > 0 || akstar_M.size() > 0", "has K*0 or K̄*0")

    df.Report().Print()

    Br2Save = kstar_vars_toSave + akstar_vars_toSave + ["thrust", "thrust_truth"]

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

    df = R.RDataFrame(args.tree, input_path)
    result = globals()[channel](df)
    if isinstance(result, tuple):
        df, branches_to_save = result
        existing_cols = set(df.GetColumnNames())
        branches_to_save = [b for b in branches_to_save if b in existing_cols]
        print(f"Saving {len(branches_to_save)} branches: {branches_to_save}")
        df.Snapshot(args.tree, output_path, branches_to_save)
    else:
        df = result
        print("Saving all branches")
        df.Snapshot(args.tree, output_path)

if __name__ == "__main__":
    main()
    #batch_process()