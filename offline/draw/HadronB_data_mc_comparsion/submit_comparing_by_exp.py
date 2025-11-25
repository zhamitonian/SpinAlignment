import subprocess
import os

for exp in [7, 11, 13, 15, 17, 19, 23, 25, 27, 31, 33, 35,
                      37, 39, 41, 43, 45, 47, 49, 51, 55, 61, 63, 65, 67, 69, 71]:
    output_dir = "./"
    os.makedirs(os.path.join(output_dir, "log"), exist_ok=True)
    log_file = os.path.join(output_dir, "log", f"exp{exp}.log")
    cmd = f"bsub -q s -oo {log_file} python3 by_exp.py {output_dir} {exp}"
    subprocess.run(cmd, shell=True, check=True)