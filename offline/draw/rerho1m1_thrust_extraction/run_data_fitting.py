#!/usr/bin/env python3
import subprocess

def main():
    first_run = False

    rootFile = "./reco_output_data.root"
    output_dir = "./fit_results/data_fit/"
    fitting_script = "./data_fitting.py"

    vec_branches = ["phi_M", "phi_p", "phi_costheta", "phi_phi", "phi_z", "phi_helicity_phi", "phi_helicity_angle",  "phi_thrust_helicity_costheta", "phi_thrust_helicity_phi"]

    def submit_cmd(bin_list, if_lsf=True, signle_thread=False):
        if if_lsf:
            cmd = ["bsub", "-q", "s"]
        else:
            cmd = []

        bins_args = [",".join(str(i) for i in bin_list)] if signle_thread else [str(i) for i in bin_list]
        base_cmd = cmd + [
            "python3", fitting_script,
            "-i", rootFile,
            "-od", output_dir,
            "--batch", 
            "-vecBr", ",".join(vec_branches),
            "--binned"
        ]
        for bins_arg in bins_args:
            subprocess.run(base_cmd + ["--bins", bins_arg], check=True)

    if first_run:
        submit_cmd(list(range(0,300)), if_lsf=True)
    else:
        refit_bins = [130, 235, 281, 301]
        submit_cmd(refit_bins, if_lsf=False, signle_thread=True)

    return

if __name__ == "__main__":
    main()