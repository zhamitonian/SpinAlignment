import subprocess


def main():
    success = [10,110,117,120,126,130,131,132,134,135,136,138,139,140,142,143,144,147,148,149,151,152,153,154,155,156,157,158,159,]
    success += [12,13,18,19,121,122,123,124,126,128,133,137,141,145,146,150,]
    success += [17,39,93,100,101,103,104,105,107,111,112,113,114,115,116,118,125,127,129,]
    success += [95,106,]
    success = []

    refit_bins = list(range(70,80)) 
    refit_bins = [72,73,74,75,76]
    refit_bins = [86]

    first_run = False
    refit = True

    rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_data_svd2/svd2_reco_processed.root"

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
            "--binned"
        ]
        subprocess.run(cmd, check=True)

    if first_run :
        for i in range(0, 200):
            submit_cmd(i)
    elif refit == False:
        for i in range(10,160):
            if i in success:
                continue
            submit_cmd(i, if_lsf=False)
    else:
        for i in refit_bins:
            submit_cmd(i, if_lsf=False)

    return

if __name__ == "__main__":
    main()