import subprocess

for i in range(200):
    cmd = [
        "bsub", "-q", "s",
        "./start_fitting.py",
        "-i", "/gpfs/group/belle2/users2022/luruihua/for_wangz/data_gMC_belle1/2025-11-12_SpinAlignment_gMC/continuum_reco_processed.root",
        "-od", "./new_images/fitting/",
        "-BrN", "phi_M",
        "--batch", "--bins", str(i),
    ]
    #print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    