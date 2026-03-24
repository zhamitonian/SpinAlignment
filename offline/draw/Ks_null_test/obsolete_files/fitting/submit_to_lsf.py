import subprocess

#failed_bins = [20,21,22,23,28,29,30,32,33,40,41,42,45,46,48,53,57,62,68,70,75,77,79,83,84,96,98,108,121,134,139,140,152,154,155,157]
#failed_bins = [20,21,22,23,28,29,30,32,33,41,45,53,68,70,77,79,84,96,108,121,134,139,152,154,155]
#failed_bins = [20,21,22,28,29,30,32,33,45,68,108]
#failed_bins = [20,21,22,28,29,30,32,33,45]
#failed_bins = [20,21,29,30,32]
#failed_bins = [20]
#failed_bins = [29,32]
failed_bins = [29]

#for i in range(150):
    #i += 10
if 1:
    cmd = [
        #"bsub", "-q", "s",
        "./start_fitting.py",
        "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.0.0_data/exp55_reco_processed.root",
        "-od", "../images/fit/",
        #"--batch", "--bins", str(i),
        "--batch", "--bins", ",".join([str(x) for x in failed_bins]),
        "-BrN", "Ks_M",
        "--binned"
    ]
    subprocess.run(cmd, check=True)
    