import subprocess

failed_bins = [i + 10 for i in range(10)]  +\
      [20,21,22,28,29,72,78,80,84,89,103,104,114,120,123,125,127,130,133,135,137,138,139,142,143,] +\
      [150+i for i in range(10)] 
failed_bins =  [10,13,14,16,17] + [20,21,22,28,29,114,120,123,138,] + [150+ i for i in range(10)] 
failed_bins =  [10,13,14,16,17] + [20,28,29,120,123,138,] + [150+ i for i in range(10)] 
failed_bins =  [10,17] + [20,28,29,120,123,138,] + [150+ i for i in range(10)] 
failed_bins =  [20,28,29,120,123,138,] + [150+ i for i in range(10)] 
failed_bins =  [20,120,123,138,] + [150+ i for i in range(10)] 
failed_bins =  [120,123,138,] 


success = []

cmd = ["./start_fitting.py",
        "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.1.0_data/exp55_reco_processed.root",
        "-od", "../images/fit_dt_e55_v2.1.0/",
        "--batch", "--bins", ",".join([str(x) for x in range(10,140) if x not in success]),
        "-BrN", "Ks_M",
        ]
subprocess.run(cmd, check=True)   

'''
cmd = ["./start_fitting.py",
        "-i", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.1.0_data/exp55_reco_processed.root",
        "-od", "../images/fit_dt_e55_v2.1.0/",
        "--batch", "--bins", ",".join([str(x) for x in failed_bins]),
        "-BrN", "Ks_M",
        ]
subprocess.run(cmd, check=True)   
'''
