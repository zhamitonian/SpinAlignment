#!/usr/bin/env python3
import subprocess

def main():
    first_run = False

    rootFile = "rootFiles/sig_isSignal_v2.3.3.root"

    def submit_cmd(bin_list, if_lsf=True, signle_thread=False):
        if if_lsf:
            cmd = ["bsub", "-q", "s"]
        else:
            cmd = []

        bins_args = [",".join(str(i) for i in bin_list)] if signle_thread else [str(i) for i in bin_list]
        base_cmd = cmd + [
            #"./fitting_func.py",
            "./data_fitting.py",
            "-i", rootFile,
            "-od", "fit_results/test/",
            "--batch", "-BrN", "Ks_M",
            #"-WgN", "Ks_weight",
            "--binned"
        ]
        for bins_arg in bins_args:
            subprocess.run(base_cmd + ["--bins", bins_arg], check=True)

    if first_run:
        submit_cmd(list(range(10,150), if_lsf=True))
    else:
        #refit_bins = list(range(30,40))
        refit_bins = [31]
        submit_cmd(refit_bins, if_lsf=False, signle_thread=True)

    return

if __name__ == "__main__":
    main()