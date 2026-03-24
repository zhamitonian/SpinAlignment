import subprocess

for i in range(1000):
    cmd = [
        "bsub", "-q", "s",
        "./phi_data_mc_comparsion.py",
        "--bin", str(i)
    ]
    subprocess.run(cmd, check=True)