import subprocess

success = [60,70,71,73,74,75,76,77,78,79,80,82,83,84,85,86,148,149] + list(range(89, 146))
success += [12,13,14,15,16,17,]
success += [50,51,52,64,67,69,81,82,87,88]
success += [11,30,39,40,49,57,58,61,63,65,66,72,]
success += [10,]
success += [20,25,28,29,31,35,38,53,54,59,]
success += [18,19,]
success += [47,48,62,68,]
success += [41,]
success += [32,42,45,]
success += [21,22,23,26,37,44,55,]
success += [24,27,56]
success += [33,34,36,43,46]


cmd = ["./start_fitting.py",
        "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.1.0_qqbar/exp55_reco_processed.root",
        "-od", "../images/fit_mc_e55_v2.1.0/",
        #"--batch", "--bins", ",".join([str(x) for x in range(10,140) if x not in success]),
        "--batch", "--bins", "147",
        "-BrN", "Ks_M",
        "--binned"
        ]
subprocess.run(cmd, check=True)   


'''
for i in [x+10 for x in range(140)]:
    cmd = [
        "bsub", "-q", "s",
        "./start_fitting.py",
        "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.1.0_qqbar/exp55_reco_processed.root",
        "-od", "../images/fit_mc_e55_v2.1.0/",
        "--batch", "--bins", str(i),
        "-BrN", "Ks_M",
        #"--binned"
    ]
    subprocess.run(cmd, check=True)
'''

