import subprocess


def main():
    success = []
    success += []

    refit_bins = list(range(60, 70))
    refit_bins = [81]

    first_run = True
    refit = True

    rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar_svd2/svd2_reco_truth_matched.root"

    def submit_cmd(i, if_lsf=True):
        if if_lsf:
            cmd = ["bsub", "-q", "s"]
        else:
            cmd = []

        cmd += [
            "./start_fitting.py",
            "-i", rootFile,
            "-od", "output/",
            "--batch", "--bins", str(i),
            "-BrN", "Ks_M",
            "-WgN", "Ks_weight",
            "--binned"
        ]
        subprocess.run(cmd, check=True)

    if first_run :
        for i in range(0, 200):
            submit_cmd(i)
    elif refit == False:
        for i in range(60,130):
            if i in success:
                continue
            submit_cmd(i, if_lsf=False)
    else:
        for i in refit_bins:
            submit_cmd(i, if_lsf=False)

    return

if __name__ == "__main__":
    main()