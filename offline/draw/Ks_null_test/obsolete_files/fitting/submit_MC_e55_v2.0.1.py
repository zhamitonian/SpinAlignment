import subprocess

failed_bins = [29]

for i in range(150):
    i += 10
#if 1:
    cmd = [
        "bsub", "-q", "s",
        "./start_fitting.py",
        "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.1.0_qqbar/exp55_reco_processed.root",
        "-od", "../images/fit_MC_e55_v2.1.0/",
        "--batch", "--bins", str(i),
        #"--batch", "--bins", ",".join([str(x) for x in failed_bins]),
        "-BrN", "Ks_M",
        #"--binned"
    ]
    subprocess.run(cmd, check=True)
    