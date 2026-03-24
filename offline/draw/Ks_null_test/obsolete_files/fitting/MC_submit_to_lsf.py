import subprocess

failed_bins = [i+10 for i in range(150)]
failed_bins = [11,12,20,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,61,69,76,81,82,87,92,94,98,101,104,105,106,107,108,141,143,144,145,148,149,151,153,156,159]
failed_bins = [11,20,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,61,82,92,143,144,149,151]
failed_bins = [11,20,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,61,82,144,151]
#for i in range(150):
    #i += 10
if 1:
    cmd = [
        #"bsub", "-q", "s",
        "./start_fitting_MC.py",
        "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.0.0_qqbar/exp55_reco_processed.root",
        "-od", "../images/MC_fit/",
        #"--batch", "--bins", str(i),
        "--batch", "--bins", ",".join([str(x) for x in failed_bins]),
        "-BrN", "Ks_M",
        "--binned"
    ]
    subprocess.run(cmd, check=True)
    