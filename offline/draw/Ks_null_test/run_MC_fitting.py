#!/usr/bin/env python3
import subprocess

def main():
    first_run = True

    # MC fitting(with bkg)
    rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.3_qqbar_svd2/svd2_only/svd2_st0_reco_processed.root"
    output_dir = "./fit_results/MC_fit/"
    fitting_script = "./MC_fitting.py"

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
            "--batch", "-BrN", "Ks_M",
            "-WgN", "Ks_weight",
            "--binned"
        ]
        for bins_arg in bins_args:
            subprocess.run(base_cmd + ["--bins", bins_arg], check=True)

    if first_run:
        submit_cmd(list(range(150,160)), if_lsf=True)
    else:
        refit_bins = [152]
        submit_cmd(refit_bins, if_lsf=False, signle_thread=True)

    return

if __name__ == "__main__":
    main()