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

def reco(df:R.RDataFrame):
    """
    inclusive phi/ Kstar , rho00 measurement, reconstruction level processing
    """
    # phi
    tools = RDF_process()
    df = tools.save_all_pairs(df, 
                         p1_branches=["kp_px_cms","kp_py_cms","kp_pz_cms"],
                         p2_branches=["km_px_cms","km_py_cms","km_pz_cms"],
                         mass=(0.493677, 0.493677),
                         particle_name=("phi","phi_kp","phi_km"))

    df = df.Define("phi_M", "sqrt(phi_E*phi_E - phi_px*phi_px - phi_py*phi_py - phi_pz*phi_pz)")
    phi_mass_window = (0.99, 1.15) ## KK threshold 0.98736
    mass_condition = f"(phi_M > {phi_mass_window[0]}) && (phi_M < {phi_mass_window[1]})"
    df = df.Redefine("phi_M", f"phi_M[{mass_condition}]")

    # Kstar
    df = tools.save_all_pairs(df,
                              p1_branches=["kp_px_cms","kp_py_cms","kp_pz_cms"],
                              p2_branches=["pim_px_cms","pim_py_cms","pim_pz_cms"],
                              mass=(0.493677, 0.13957039),
                              particle_name=("Kstar","Kstar_kp","Kstar_pim"))

    df = tools.save_all_pairs(df,
                                p1_branches=["km_px_cms","km_py_cms","km_pz_cms"],
                                p2_branches=["pip_px_cms","pip_py_cms","pip_pz_cms"],
                                mass=(0.493677, 0.13957039),
                                particle_name=("Kstar_bar","Ksatr_bar_km","Kstar_bar_pip"))

    df = df.Define("Kstar_M", "sqrt(Kstar_E*Kstar_E - Kstar_px*Kstar_px - Kstar_py*Kstar_py - Kstar_pz*Kstar_pz)")
    df = df.Define("Kstar_bar_M", "sqrt(Kstar_bar_E*Kstar_bar_E - Kstar_bar_px*Kstar_bar_px - Kstar_bar_py*Kstar_bar_py - Kstar_bar_pz*Kstar_bar_pz)")
    df = df.Redefine("Kstar_M", "Kstar_M[0.69 < Kstar_M && Kstar_M < 1.09]")
    df = df.Redefine("Kstar_bar_M", "Kstar_bar_M[0.69 < Kstar_bar_M && Kstar_bar_M < 1.09]")

    # D0
    df = tools.save_all_pairs(df, 
                              p1_branches=["km_px_cms","km_py_cms","km_pz_cms"],
                              p2_branches=["pip_px_cms","pip_py_cms","pip_pz_cms"],
                              mass=(0.493677, 0.13957039),
                              particle_name=("D0","D0_km","D0_pip"))
    
    df = tools.save_all_pairs(df,
                              p1_branches=["kp_px_cms","kp_py_cms","kp_pz_cms"],
                              p2_branches=["pim_px_cms","pim_py_cms","pim_pz_cms"],
                              mass=(0.493677, 0.13957039),
                              particle_name=("D0_bar","D0_bar_kp","D0_bar_pim"))    
    # D0 1.8648
    df = df.Define("D0_M", "sqrt(D0_E*D0_E - D0_px*D0_px - D0_py*D0_py - D0_pz*D0_pz)")
    df = df.Define("D0_bar_M", "sqrt(D0_bar_E*D0_bar_E - D0_bar_px*D0_bar_px - D0_bar_py*D0_bar_py - D0_bar_pz*D0_bar_pz)")
    df = df.Redefine("D0_M", "D0_M[1.6 < D0_M && D0_M < 2.15]")
    df = df.Redefine("D0_bar_M", "D0_bar_M[1.6 < D0_bar_M && D0_bar_M < 2.15]")

    Br2Save = ["phi_M", "Kstar_M", "Kstar_bar_M", "D0_M", "D0_bar_M"]

    df.Report().Print()

    return (df, Br2Save)

def main():
    parser = argparse.ArgumentParser(description="set IO")
    parser.add_argument("which_channel",
                        choices=["reco"],
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
        print(f"Saving {len(branches_to_save)} branches: {branches_to_save}")
        df.Snapshot(args.tree , output_path, branches_to_save)
    else:
        df = result
        print("Saving all branches")
        df.Snapshot(args.tree, output_path)

if __name__ == "__main__":
    main()
    #batch_process()