import subprocess

for i in [i+10 for i in range(140)]:
    cmd = [
        "bsub", "-q", "s",
        "./start_comparing.py",
        "--bin", str(i)
    ]
    subprocess.run(cmd, check=True)