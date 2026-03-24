import subprocess


def main():
    success = [70, 71, 77, 79, 80, 81, 83, 84, 85, 87, 88, 89, 90, 91, 92, 94, 96, 97, 98, 100, 101, 102, 105, 106, 107, 109, 111, 112, 113, 114, 115, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129]
    success += [78, 82, 95, 99, 103, ]
    success += [61, 62, 66, 67, 68, 69, 72, 75, 82, 86, 93, 104, 108, 110, 116]
    success+= [60, 63, 65, 73, 74, 76]

    first_run = False
    refit = True

    refit_bins = [75]

    if first_run :
        for i in range(200):
            cmd = [
                "bsub", "-q", "s",
                "./start_fitting.py",
                "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar/exp55_reco_processed_sig.root",
                "-od", "output/",
                "--batch", "--bins", str(i),
                "-BrN", "Ks_M",
                "-WgN", "Ks_weight",
                #"--binned"
            ]
            subprocess.run(cmd, check=True)
    elif refit == False:
        for i in range(60,130):
            if i in success:
                continue
            cmd = [
                    "bsub", "-q", "s",
                    "./start_fitting.py",
                    "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar/exp55_reco_processed_sig.root",
                    "-od", "output/",
                    #"--batch", "--bins", ",".join([str(x) for x in range(60,130) if x not in success]),
                    "--batch", "--bins", str(i),
                    "-BrN", "Ks_M",
                    "-WgN", "Ks_weight",
                    #"--binned"
                ]
            subprocess.run(cmd, check=True)
    else:
        for i in refit_bins:
            cmd = [
                "bsub", "-q", "s",
                "./start_fitting.py",
                "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar/exp55_reco_processed_sig.root",
                "-od", "output/",
                #"--batch", "--bins", ",".join([str(x) for x in range(60,130) if x not in success]),
                "--batch", "--bins", str(i),
                "-BrN", "Ks_M",
                "-WgN", "Ks_weight",
                #"--binned"
            ]
            subprocess.run(cmd, check=True)

    return

if __name__ == "__main__":
    main()