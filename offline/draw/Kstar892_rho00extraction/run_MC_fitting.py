#!/usr/bin/env python3
import subprocess

def main():
    first_run = False

    rootFile = "/gpfs/group/belle/users/dues/data_gMC_belle1/Kstar892SpinAlignment_v1.0.0_qqbar_svd2_duxs/reco.root"
    output_dir = "./fit_results/MC_fit/"
    fitting_script = "./MC_fitting.py"

    vec_branches = ["Kstar_z", "Kstar_M", "Kstar_helicity_angle", "Kstar_weight"]

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
        submit_cmd(list(range(50, 170)), if_lsf=True)
    else:
        refit_bins = [114]
        submit_cmd(refit_bins, if_lsf=False, signle_thread=True)

    return

if __name__ == "__main__":
    main()