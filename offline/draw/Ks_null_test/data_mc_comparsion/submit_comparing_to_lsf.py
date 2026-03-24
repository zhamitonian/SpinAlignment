import subprocess

for i in range(150, 160):
    cmd = [
        "bsub", "-q", "s",
        "./start_comparing.py",
        "--bin", str(i)
    ]
    subprocess.run(cmd, check=True)