import subprocess

for i in range(10):
    cmd = [
        "bsub", "-q", "s",
        "./start_fitting.py",
        "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/PhiSpinAlignment_v2.1.0_data/continuum_reco_processed.root",
        "-od", "./test_images/",
        "-BrN", "phi_M",
        "--batch", "--bins", str(i),
        "--binned"
    ]
    #print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    