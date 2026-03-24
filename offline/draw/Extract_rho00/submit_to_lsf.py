import subprocess

for i in range(200):
    cmd = [
        "bsub", "-q", "s",
        "./start_fitting.py",
        "-i", "data_processed_test.root",
        "-od", "./images_test/",
        "-BrN", "phi_M",
        "--batch", "--bins", str(i),
    ]
    #print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    