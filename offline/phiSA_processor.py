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

def SA(df:R.RDataFrame)->R.RDataFrame: # SpinAlignment
    #df = df.Filter("nGood > 2", "good trk")
    sqrts = 10.516469955444336 
    #df = df.Filter(f"Evis_cms > {sqrts}*0.2", "Evis cut")
    df = df.Filter("Evis_cms > 7", "Evis cut")
    df = df.Filter(f"Esum_intheta_cms > {sqrts}*0.1 && Esum_intheta_cms < {sqrts}*0.8", "Esum cut")
    df = df.Filter(f"abs(PzSum_cms) < 0.5*{sqrts}","Pz Balance")
    df = df.Filter(f"PrimeVr <1.5", "primary Vr cut")
    df = df.Filter("nECL_inCostheta > 1", "nECL cut") 
    df = df.Filter(f"Esum_cms > 0.18*{sqrts} || HeavyJetMass > 1.8 || HeavyJetMass > 0.25*Evis_cms", "HJM cut")
    df = df.Filter("Esum_cms < nECL", "Esum_cms/nECL < 1")

    df = df.Filter("abs(PrimeVz) < 3.5", "primary Vz cut")

    df.Report().Print()
    
    return df

def test_cut(df:R.RDataFrame, cut:Optional[str]=None)->R.RDataFrame:
     df = df.Filter("Evis_cms > 7", "Evis cut")
     df = df.Filter(f"HeavyJetMass > 1.8 || HeavyJetMass > 0.25*Evis_cms", "HJM cut")
     df = df.Filter("abs(PrimeVz) < 3.5", "primary Vz cut")
     df = df.Filter("PrimeVr <1.5", "primary Vr cut")

     df.Report().Print()

     return df
    
#def phi_test(df:R.RDataFrame)->Tuple[R.RDataFrame: list[str]]:
def phi_reco(df:R.RDataFrame, hasPhi:bool=False):
    """
    inclusive phi , rho00 measurement, reconstruction level processing
    """
    tools = RDF_process()
    df = tools.save_all_pairs(df, 
                         p1_branches=["kp_px_cms","kp_py_cms","kp_pz_cms"],
                         p2_branches=["km_px_cms","km_py_cms","km_pz_cms"],
                         mass=(0.493677, 0.493677),
                         particle_name=("phi","kp","km"))

    df = df.Define("phi_M", "sqrt(phi_E*phi_E - phi_px*phi_px - phi_py*phi_py - phi_pz*phi_pz)")

    phi_mass_window = (0.99, 1.15) ## KK threshold 0.98736
    mass_condition = f"(phi_M > {phi_mass_window[0]}) && (phi_M < {phi_mass_window[1]})"
    
    df = df.Redefine("phi_E", f"phi_E[{mass_condition}]")
    df = df.Redefine("phi_px", f"phi_px[{mass_condition}]")
    df = df.Redefine("phi_py", f"phi_py[{mass_condition}]")
    df = df.Redefine("phi_pz", f"phi_pz[{mass_condition}]")
    df = df.Redefine("kp_index", f"kp_index[{mass_condition}]")
    df = df.Redefine("km_index", f"km_index[{mass_condition}]")
    df = df.Redefine("phi_helicity_angle", f"phi_helicity_angle[{mass_condition}]")

    # the variable used to perform selction showld be redefined at last
    df = df.Redefine("phi_M", f"phi_M[{mass_condition}]")
    
    #df = df.Define("phi_xp", "sqrt(phi_px*phi_px + phi_py*phi_py)/sqrt(sqrts*sqrts*0.25 - 1.0195*1.0195)") # previous calculation wrong !!!
    df = df.Define("phi_xp", "sqrt(phi_E*phi_E - phi_M*phi_M)/sqrt(sqrts*sqrts*0.25 - 1.0195*1.0195)")
    df = df.Define("phi_z", "2*phi_E/sqrts")
    df = tools.calculate_pt_toAxis(df, 
                                  particle=["phi_px","phi_py","phi_pz"],
                                  axis = ["TMath::ACos(thrust[1])", "thrust[2]"],
                                  particle_name="phi",axis_name="thrust")

    df = tools.convert_cartesian_to_spherical(df, particles=["phi"])
    df = tools.convert_cartesian_to_spherical(df, particles=["kp", "km"], px_branch = "px_cms", py_branch = "py_cms", pz_branch = "pz_cms")

    Br2Save = ["phi_M", "phi_z", "phi_xp","phi_helicity_angle", "phi_thrust_pt"]
    Br2Save += ["phi_costheta", "phi_phi", "phi_p", "kp_p", "kp_costheta", "kp_phi", "km_p", "km_costheta", "km_phi"]
    Br2Save += ["kp_index", "km_index"]

    if hasPhi:
        df = df.Filter("n_phi_truth > 0", "has truth phi") 

    df.Report().Print()

    return (df, Br2Save)
    #return df

def phi_truth(df:R.RDataFrame, kaon_inBarrel:bool=False):
    """
    process truth tree, retrieve coordinate generate leavel phi's z, cos*theta ...
    """
    tools = RDF_process()

    if kaon_inBarrel:
        kaon_inBarrel = f"cos(kp_theta*{pi}/180) < 0.842 && cos(kp_theta*{pi}/180) > -0.511 && cos(km_theta*{pi}/180) < 0.842 && cos(km_theta*{pi}/180) > -0.511 && kp_pt_truth > 0.05 && km_pt_truth > 0.05"
        for var in ["px", "py", "pz"]:
            for kaon in ["kp", "km"]:
                df = df.Redefine(f"{kaon}_{var}_cms_truth", f"{kaon}_{var}_cms_truth[{kaon_inBarrel}]")
        
    df = tools.save_all_pairs(df, 
                         p1_branches=["kp_px_cms_truth","kp_py_cms_truth","kp_pz_cms_truth"],
                         p2_branches=["km_px_cms_truth","km_py_cms_truth","km_pz_cms_truth"],
                         mass=(0.493677, 0.493677),
                         particle_name=("phi","kp","km"),
                         cross_mode= False)

    df = df.Define("phi_M", "sqrt(phi_E*phi_E - phi_px*phi_px - phi_py*phi_py - phi_pz*phi_pz)")
    sqrts = 10.516469955444336 
    df = df.Define("phi_z", f"2*phi_E/{sqrts}")
    df = df.Define("phi_xp", f"sqrt(phi_E*phi_E - phi_M*phi_M)/sqrt({sqrts}*{sqrts}*0.25 - 1.0195*1.0195)")
    
    Br2Save = ["phi_M", "phi_z", "phi_xp", "phi_helicity_angle"]
    Br2Save += ["kp_px_cms_truth", "kp_py_cms_truth", "kp_pz_cms_truth", "kp_E_cms_truth"]
    Br2Save += ["km_px_cms_truth", "km_py_cms_truth", "km_pz_cms_truth", "km_E_cms_truth"]
    df.Report().Print()

    return (df, Br2Save)


def phi_reco_truth_matched(df:R.RDataFrame, phi_mass_window:Tuple[float,float]=(0.99, 1.15)):
    """
    using truth matched generic MC, for z, cos*theta, thrust, pt  resolution checking 
    """
    tools = RDF_process()
    df = df.Filter("n_phi_truth > 0", "has phi")

    for var in ["px_cms", "py_cms", "pz_cms", "E_cms"]:
        for kaon in ["kp", "km"]:
            df = df.Redefine(f"{kaon}_{var}", f"{kaon}_{var}[{kaon}_isSignal]")

    df = tools.save_all_pairs(df, 
                         p1_branches=["kp_px_cms","kp_py_cms","kp_pz_cms"],
                         p2_branches=["km_px_cms","km_py_cms","km_pz_cms"],
                         mass=(0.493677, 0.493677),
                         particle_name=("phi","kp","km"))

    df = df.Define("phi_M", "sqrt(phi_E*phi_E - phi_px*phi_px - phi_py*phi_py - phi_pz*phi_pz)")
    
    #phi_mass_window = (0.99, 1.15) ## KK threshold 0.98736
    mass_condition = f"(phi_M > {phi_mass_window[0]}) && (phi_M < {phi_mass_window[1]})"
    
    # Apply mass window filter to all phi-related variables
    df = df.Redefine("phi_E", f"phi_E[{mass_condition}]")
    df = df.Redefine("phi_px", f"phi_px[{mass_condition}]")
    df = df.Redefine("phi_py", f"phi_py[{mass_condition}]")
    df = df.Redefine("phi_pz", f"phi_pz[{mass_condition}]")
    df = df.Redefine("phi_helicity_angle", f"phi_helicity_angle[{mass_condition}]")
    df = df.Redefine("phi_M", f"phi_M[{mass_condition}]")
    
    df = df.Define("phi_z", "2*phi_E/sqrts")
    df = df.Define("phi_xp", "sqrt(phi_E*phi_E - phi_M*phi_M)/sqrt(sqrts*sqrts*0.25 - 1.0195*1.0195)")
    df = tools.calculate_pt_toAxis(df, 
                                  particle=["phi_px","phi_py","phi_pz"],
                                  axis = ["TMath::ACos(thrust[1])", "thrust[2]"],
                                  particle_name="phi",axis_name="thrust")
    df = df.Redefine("thrust", "thrust[0]")

    # ----------- generation info --------------------
    df = tools.save_all_pairs(df, 
                        p1_branches=["kp_px_cms_gen","kp_py_cms_gen","kp_pz_cms_gen"],
                        p2_branches=["km_px_cms_gen","km_py_cms_gen","km_pz_cms_gen"],
                        mass=(0.493677, 0.493677),
                        particle_name=("phi_gen","kp_gen","km_gen"),
                        cross_mode= False) 
    df = tools.calculate_pt_toAxis(df, 
                        particle=["phi_gen_px","phi_gen_py","phi_gen_pz"],
                        axis = ["TMath::ACos(thrust_truth[1])", "thrust_truth[2]"],
                        particle_name="phi_gen",axis_name="thrust")

    df = df.Redefine("thrust_truth", "thrust_truth[0]")
    df = df.Define("phi_gen_z", f"2*phi_gen_E/sqrts")
    df = df.Define("phi_gen_M", "sqrt(phi_gen_E*phi_gen_E - phi_gen_px*phi_gen_px - phi_gen_py*phi_gen_py - phi_gen_pz*phi_gen_pz)")
    df = df.Define("phi_gen_xp", "sqrt(phi_gen_E*phi_gen_E - phi_gen_M*phi_gen_M)/sqrt(sqrts*sqrts*0.25 - 1.0195*1.0195)")

    # kaon isSignal but there are some misidentified phi , suppressed by mass window and n phi cut , temporarily not care the still exist misidentified phi
    df = df.Filter("phi_M.size() == phi_gen_M.size()", "match reconstructed and generated phi count") 
    df = df.Filter("phi_M.size() > 0", "has phi after match")

    df = tools.convert_cartesian_to_spherical(df, particles=["phi"])
    df = tools.convert_cartesian_to_spherical(df, particles=["kp", "km"], px_branch = "px_cms", py_branch = "py_cms", pz_branch = "pz_cms")

    Br2Save = ["phi_M", "phi_gen_M", "phi_z", "phi_gen_z", "phi_xp", "phi_gen_xp", "phi_helicity_angle", "phi_gen_helicity_angle", "phi_thrust_pt"]
    #Br2Save += ["phi_thrust_pt", "phi_gen_thrust_pt","thrust", "thrust_truth"]
    #Br2Save += ["kp_E_cms", "kp_E_cms_gen", "km_E_cms", "km_E_cms_gen"]
    Br2Save += ["phi_costheta", "phi_phi", "phi_p", "kp_p", "kp_costheta", "kp_phi", "km_p", "km_costheta", "km_phi"]

    df.Report().Print()

    return (df, Br2Save)
    #return df


def phi_truth_check(df:R.RDataFrame, kaon_inBarrel:bool=False):
    """
    process truth tree, retrieve coordinate generate leavel phi's z, cos*theta ...
    """
    tools = RDF_process()

    if kaon_inBarrel:
        kaon_inBarrel = f"cos(kp_theta*{pi}/180) < 0.842 && cos(kp_theta*{pi}/180) > -0.511 && cos(km_theta*{pi}/180) < 0.842 && cos(km_theta*{pi}/180) > -0.511 && kp_pt_truth > 0.05 && km_pt_truth > 0.05"
        for var in ["px", "py", "pz"]:
            for kaon in ["kp", "km"]:
                df = df.Redefine(f"{kaon}_{var}_cms_truth", f"{kaon}_{var}_cms_truth[{kaon_inBarrel}]")
        
    df = tools.save_all_pairs(df, 
                         p1_branches=["kp_px_cms_truth","kp_py_cms_truth","kp_pz_cms_truth"],
                         p2_branches=["km_px_cms_truth","km_py_cms_truth","km_pz_cms_truth"],
                         mass=(0.493677, 0.493677),
                         particle_name=("phi","kp","km"),
                         cross_mode= True)

    df = df.Define("phi_M", "sqrt(phi_E*phi_E - phi_px*phi_px - phi_py*phi_py - phi_pz*phi_pz)")

    phi_mass_window = (0.99, 1.15) ## KK threshold 0.98736
    mass_condition = f"(phi_M > {phi_mass_window[0]}) && (phi_M < {phi_mass_window[1]})"
    
    df = df.Redefine("phi_E", f"phi_E[{mass_condition}]")
    df = df.Redefine("phi_px", f"phi_px[{mass_condition}]")
    df = df.Redefine("phi_py", f"phi_py[{mass_condition}]")
    df = df.Redefine("phi_pz", f"phi_pz[{mass_condition}]")
    df = df.Redefine("kp_index", f"kp_index[{mass_condition}]")
    df = df.Redefine("km_index", f"km_index[{mass_condition}]")
    df = df.Redefine("phi_helicity_angle", f"phi_helicity_angle[{mass_condition}]")

    # the variable used to perform selction showld be redefined at last
    df = df.Redefine("phi_M", f"phi_M[{mass_condition}]")

    sqrts = 10.516469955444336 
    df = df.Define("phi_z", f"2*phi_E/{sqrts}")
    df = df.Define("phi_xp", f"sqrt(phi_E*phi_E - phi_M*phi_M)/sqrt({sqrts}*{sqrts}*0.25 - 1.0195*1.0195)")
    
    Br2Save = ["phi_M", "phi_z", "phi_xp", "phi_helicity_angle"]
    df.Report().Print()

    return (df, Br2Save)


def temp_test(df:R.RDataFrame):
    """
    inclusive phi , rho00 measurement, reconstruction level processing
    """
    tools = RDF_process()
    df = tools.save_all_pairs(df, 
                         p1_branches=["kp_px_cms","kp_py_cms","kp_pz_cms"],
                         p2_branches=["km_px_cms","km_py_cms","km_pz_cms"],
                         mass=(0.493677, 0.493677),
                         particle_name=("phi","kp","km"))

    df = df.Define("phi_M", "sqrt(phi_E*phi_E - phi_px*phi_px - phi_py*phi_py - phi_pz*phi_pz)")

    phi_mass_window = (0.99, 1.15) ## KK threshold 0.98736
    mass_condition = f"(phi_M > {phi_mass_window[0]}) && (phi_M < {phi_mass_window[1]})"
    
    df = df.Redefine("phi_E", f"phi_E[{mass_condition}]")
    df = df.Redefine("phi_px", f"phi_px[{mass_condition}]")
    df = df.Redefine("phi_py", f"phi_py[{mass_condition}]")
    df = df.Redefine("phi_pz", f"phi_pz[{mass_condition}]")
    df = df.Redefine("kp_index", f"kp_index[{mass_condition}]")
    df = df.Redefine("km_index", f"km_index[{mass_condition}]")
    df = df.Redefine("phi_helicity_angle", f"phi_helicity_angle[{mass_condition}]")

    # the variable used to perform selction showld be redefined at last
    df = df.Redefine("phi_M", f"phi_M[{mass_condition}]")
    
    df = df.Define("phi_xp", "sqrt(phi_E*phi_E - phi_M*phi_M)/sqrt(sqrts*sqrts*0.25 - 1.0195*1.0195)")
    df = df.Define("phi_z", "2*phi_E/sqrts")

    Br2Save = ["phi_M", "phi_z", "phi_xp", "phi_helicity_angle"]

    df.Report().Print()

    return (df, Br2Save)

def main():
    parser = argparse.ArgumentParser(description="set IO")
    parser.add_argument("which_channel",
                        choices=["SA", "test_cut", "phi_reco", "phi_truth", "phi_reco_truth_matched", "temp_test", "phi_truth_check"],
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