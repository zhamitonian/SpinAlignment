import subprocess

for i in range(40):
    cmd = [
        "bsub", "-q", "s",
        "./start_fitting.py",
        "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/2025-12-03_SpinAlignment_qqbarMC/continuum_reco_processed.root",
        "-od", "./eff_Dec03/fitting/",
        "-BrN", "phi_M",
        "--batch", "--bins", str(i),
        "--binned"
    ]
    #print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    #"-i", "/gpfs/group/belle2/users2022/luruihua/for_wangz/data_gMC_belle1/2025-11-25_SpinAlignment_gMC/continuum_reco_processed.root",
    #"-od", "./eff_Nov26/fitting/",

    #"-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/2025-12-02_SpinAlignment_qqbarMC/continuum_reco_processed.root",
    #"-od", "./eff_Dec02/fitting/",
    
    #"-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/2025-12-02_SpinAlignment_qqbarMC2/continuum_reco_processed.root",
    #"-od", "./eff_Dec02_2/fitting/",