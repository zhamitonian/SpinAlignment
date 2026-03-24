import subprocess


def main():
    success = [70, 71, 77, 78, 79, 80, 81, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 104, 105, 106, 107, 110, 112, 113, 114, 115 ,116, 117, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129]
    success += [103, 108, 109, 111, 118,]
    success += [61, 62, 63, 66, 69, 72, 73, 82, 103, 108, 109, 111, 118]
    success+= [60, 64, 65, 68, 75, 76]
    success += [74]

    first_run = False
    refit = True
    refit_bins = [63]

    if first_run :
        for i in range(200):
            cmd = [
                "bsub", "-q", "s",
                "./start_fitting.py",
                "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar/exp55_reco_processed_sig.root",
                "-od", "output/",
                "--batch", "--bins", str(i),
                "-BrN", "Ks_M",
                #"--binned"
            ]
            subprocess.run(cmd, check=True)
    elif refit == False:
        for i in range(60,130):
            if i in success:
                continue
            cmd = [
                    #"bsub", "-q", "s",
                    "./start_fitting.py",
                    "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar/exp55_reco_processed_sig.root",
                    "-od", "output/",
                    "--batch", "--bins", str(i),
                    "-BrN", "Ks_M",
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
                    "--batch", "--bins", str(i),
                    "-BrN", "Ks_M",
                    #"--binned"
                ]
            subprocess.run(cmd, check=True)

    return

if __name__ == "__main__":
    main()