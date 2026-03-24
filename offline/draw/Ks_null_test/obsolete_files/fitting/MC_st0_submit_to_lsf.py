import subprocess

#failed_bins = [i+10 for i in range(150)]
#for i in range(150):
    #i += 10
failed_bins =[10,11,12,13,14,15,16,18,19,20,21,22,24,25,26,27,29,32,33,34,35,36,37,38,40,41,42,43,44,45,46,47,48,52,54,66,72,79,86,91,101,124,141,142,145,150,151,152,154,157,]
failed_bins =[11,13,18,19,20,21,22,24,25,26,27,29,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,52,54,66,72,79,86,101,141,142,145,150,151,152,154,157,]
failed_bins =[13,18,19,20,21,22,24,25,26,27,29,32,33,34,35,36,37,38,44,46,47,48,141,142,145,150,151,152,154,157,]
failed_bins = [20,24,25,26,27,29,32,33,34,35,36,37,44,46,47,48,141,142,145,150,151,152,154,157]
failed_bins = [20,24,25,26,29,32,34,35,36,37,141,142,145,150,151,157]
failed_bins = [24,25,26,29,32,34,35,36,141,142,145,150,151,157,]
failed_bins = [24,36,141,142,145,150,151,157]
failed_bins = [24,36]
    
if 1:
    cmd = [
        #"bsub", "-q", "s",
        "./start_fitting_MC_st0.py",
        "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.0.0_qqbar/4S_offres_exp55_st0_reco_processed.root",
        "-od", "../images/fit_exp55_4S_offres_st0/",
        #"--batch", "--bins", str(i),
        "--batch", "--bins", ",".join([str(x) for x in failed_bins]),
        "-BrN", "Ks_M",
        "--binned"
    ]
    subprocess.run(cmd, check=True)
    