#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_nsig(txt_file, output_file="nsig.png"):
    df = pd.read_csv(txt_file, sep='\s+', comment='#', 
                    names=['center', 'width', 'nsig', 'nsig_err', 'nsig_err2'])

    plt.figure(figsize=(10, 6))
    plt.errorbar(df['center'], df['nsig'], xerr=df['width'], yerr=df['nsig_err'], fmt='o', markersize=3, capsize=2)#, color='black')
    plt.yscale('log')
    #plt.xlabel(r'$p_t$', fontsize=20)
    plt.xlabel(r'$z$', fontsize=20)
    plt.ylabel(r'$N_{\mathrm{signal}}$', fontsize=20)
    plt.grid(True, alpha=0.3)
    plt.savefig(output_file)

    print(df)
    df.describe()
    print(f"Plot saved to: {output_file}")

def main():
    import sys
    if len(sys.argv) < 2:
        print("Usage: python plot_nsig.py <nsig_results.txt> [output_file.png]")
        return
    txt_file = sys.argv[1]
    output_file = "nsig.png"
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    plot_nsig(txt_file, output_file)


if __name__ == "__main__":
    main()