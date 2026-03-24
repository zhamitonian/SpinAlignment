import subprocess


def main():
    success = [95,106,]

    refit_bins = list(range(90, 100))
    refit_bins = [95]

    first_run = False
    refit = True

    #rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar_svd2/svd2_reco_processed.root"
    rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar_svd2/svd2_st0_reco_processed.root"

    def submit_cmd(i, if_lsf=True):
        if if_lsf:
            cmd = ["bsub", "-q", "s"]
        else:
            cmd = []

        # info will be used in batch fitting
        cmd += [
            "./start_fitting.py",
            "-i", rootFile,
            "-od", "output_st0/",
            "--batch", "--bins", str(i),
            "-BrN", "Ks_M",
            "-WgN", "Ks_weight",
            "--binned"
        ]
        subprocess.run(cmd, check=True)

    if first_run :
        for i in range(10, 160):
            submit_cmd(i)
    elif refit == False:
        for i in range(10,160):
            if i in success:
                continue
            submit_cmd(i, if_lsf=True)
    else:
        for i in refit_bins:
            submit_cmd(i, if_lsf=False)

    return

if __name__ == "__main__":
    main()