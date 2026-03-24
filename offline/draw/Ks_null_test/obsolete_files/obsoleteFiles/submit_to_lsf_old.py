import subprocess

#for i in range(200):
for i in [24]:
    cmd = [
        #"bsub", "-q", "s",
        "./start_fitting_old.py",
        "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.0.0_data/exp55_reco_processed.root",
        "-od", "./images/fit/",
        "--batch", "--bins", str(i),
        "--binned"
    ]
    #print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    