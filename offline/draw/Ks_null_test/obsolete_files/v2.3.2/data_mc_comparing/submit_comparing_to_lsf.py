import subprocess

for i in range(60, 130):
    cmd = [
        "bsub", "-q", "s",
        "./start_comparing.py",
        "--bin", str(i)
    ]
    subprocess.run(cmd, check=True)