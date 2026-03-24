import subprocess

#for i in range(1000):
#for i in [2,5]:
#for i in [0, 9, 10, 19, 20, 29, 30, 39, 49, 50, 59, 60]:
for i in [0, 39, 50]:
    cmd = [
        "bsub", "-q", "s",
        "./start_fitting.py",
        "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/PhiSpinAlignment_v2.1.0_data/continuum_reco_processed.root",
        "-od", "./fit_result/",
        "-BrN", "phi_M",
        "--batch", "--bins", str(i)
    ]
    #print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    