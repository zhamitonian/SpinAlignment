import subprocess


def main():
    success = [60, 62, 63, 64, 65, 66, 67, 68, 70, 74, 75, 77, 78, 82, 83, 85, 87, 88, 89, 90, 91, 92, 93,94, 95, 96, 97, 99, 100, 101, 102, 103, 104, 105, 107, 109, 110, 112, 113, 114, 116, 118, 119, 121, 122, 126, 127, 128, 129]
    success += [69, 71, 76, 86, 98, 124, 125]
    success+= [72, 73, 79, 80, 81, 106, 108, 111, 117]
    success+= [61, 84, 115, ]
    success += [120]
    success += [123]

    refit_bins = [72, 73, 79, 80, 81, 106, 108, 111, 117]

    first_run = False
    refit = True

    if first_run :
        for i in range(200):
            cmd = [
                "bsub", "-q", "s",
                "./start_fitting.py",
                "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_data/exp55_reco_processed.root",
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
                    "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_data/exp55_reco_processed.root",
                    "-od", "output/",
                    "--batch", "--bins", str(i),
                    "-BrN", "Ks_M",
                    #"--binned"
                ]
            subprocess.run(cmd, check=True)
    else:
        for i in refit_bins:
            cmd = [
                #"bsub", "-q", "s",
                "./start_fitting.py",
                "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_data/exp55_reco_processed.root",
                "-od", "output/",
                "--batch", "--bins", str(i),
                "-BrN", "Ks_M",
                #"--binned"
                ]
            subprocess.run(cmd, check=True) 

    return

if __name__ == "__main__":
    main()