import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def get_rebinning_boundaries(start_bin_width, start_bin_count, n_exp_bins, total_range=2.5, min_value=0.0):
    """
    Calculate bin boundaries with exponential growth
    
    Args:
        start_bin_width: Width of uniform bins at the beginning
        start_bin_count: Number of uniform bins at the beginning
        n_exp_bins: Number of exponentially growing bins
        total_range: Total range to cover
    
    Returns:
        list of bin boundaries
    """
    
    # Calculate range covered by uniform bins
    uniform_range = start_bin_width * start_bin_count
    exp_range = total_range - uniform_range
    
    print(f"Configuration:")
    print(f"  Uniform bins: {start_bin_count} bins × {start_bin_width} = {uniform_range}")
    print(f"  Exponential range: {exp_range}")
    print(f"  Total range: {total_range}\n")
    
    # Solve for exponential growth factor a
    # Sum of exponential widths = exp_range
    # width_i = start_bin_width * 10^(a*i) for i=1,2,...,n_exp_bins
    def eq(a):
        return np.sum([start_bin_width * 10**(a*i) for i in range(1, n_exp_bins + 1)]) - exp_range

    a_solution = fsolve(eq, x0=0.1)[0]
    print(f"Found exponential factor a: {a_solution:.6f}")

    # Verify
    exp_widths_sum = np.sum([start_bin_width * 10**(a_solution*i) for i in range(1, n_exp_bins + 1)])
    print(f"Check: exp widths sum = {exp_widths_sum:.6f} (target: {exp_range:.6f})")
    print(f"Total range check: {uniform_range + exp_widths_sum:.6f} (target: {total_range:.6f})\n")

    # Calculate all bin widths
    all_widths = []
    
    # Uniform widths
    for i in range(start_bin_count):
        all_widths.append(start_bin_width)
        print(f"Bin {i+1} (uniform): width = {start_bin_width:.6f}")
    
    # Exponential widths
    for i in range(1, n_exp_bins + 1):
        width = start_bin_width * 10**(a_solution * i)
        all_widths.append(width)
        print(f"Bin {start_bin_count + i} (exp): width = {width:.6f}")

    # Convert widths to boundaries (cumulative sum)
    boundaries = [min_value]
    for width in all_widths:
        boundaries.append(boundaries[-1] + width)
    
    print(f"\n{'='*60}")
    print(f"Total bins: {len(all_widths)}")
    print(f"  Uniform: {start_bin_count}")
    print(f"  Exponential: {n_exp_bins}")
    print(f"\nBin boundaries:")
    print(f"{[round(b, 4) for b in boundaries]}")
    print(f"\nFor fit_tools.py:")
    print(f"bin_var_config = [(\"variable_name\", {[round(b, 4) for b in boundaries]})]")
    
    return boundaries


if __name__ == "__main__":
    #get_rebinning_boundaries(0.125, 4, 6)
    get_rebinning_boundaries(0.075, 4, 6, 0.8, 0.2)


'''
def plot_after_rebinning():
    """Rebin data and plot comparison"""
    boundaries = get_rebinning_boundaries()

    # Read original data
    txt_file = "./test_images/nsig_results_binning_v0.txt"
    df = pd.read_csv(txt_file, sep='\s+', comment='#',
                     names=['center', 'width', 'nsig', 'nsig_err', 'nsig_err2'])
    
    # Calculate original bin edges
    df['bin_min'] = df['center'] - df['width']
    df['bin_max'] = df['center'] + df['width']

    # Split data at 0.5
    df_before = df[df['center'] < 0.5].copy()
    df_after = df[df['center'] >= 0.5].copy()
    
    print("\n" + "="*60)
    print(f"Data before 0.5: {len(df_before)} bins, {df_before['nsig'].sum():.0f} events")
    print(f"Data after 0.5: {len(df_after)} bins, {df_after['nsig'].sum():.0f} events")
    
    # Rebin data after 0.5 with exponential bins
    rebinned_data = []
    print("\n" + "="*60)
    print("Rebinning results (x >= 0.5):")
    
    # Shift boundaries to start from 0.5
    boundaries_shifted = [b + 0.5 for b in boundaries]
    
    for i in range(len(boundaries_shifted)-1):
        bin_min = boundaries_shifted[i]
        bin_max = boundaries_shifted[i+1]
        
        # Find all original bins that overlap with this new bin
        mask = (df_after['bin_max'] > bin_min) & (df_after['bin_min'] < bin_max)
        df_in_bin = df_after[mask]
        
        if len(df_in_bin) > 0:
            nsig_sum = df_in_bin['nsig'].sum()
            nsig_err_sum = np.sqrt(np.sum(df_in_bin['nsig_err']**2))
            
            rebinned_data.append({
                'bin_min': bin_min,
                'bin_max': bin_max,
                'center': (bin_min + bin_max) / 2,
                'width': (bin_max - bin_min) / 2,
                'nsig': nsig_sum,
                'nsig_err': nsig_err_sum,
                'region': 'exponential'
            })
            
            print(f"Bin {i+1}: [{bin_min:.4f}, {bin_max:.4f}], "
                  f"width={bin_max-bin_min:.4f}, N_sig = {nsig_sum:.0f} ± {nsig_err_sum:.0f}") 
    
    df_rebinned_after = pd.DataFrame(rebinned_data)
    
    # Add region label to before data
    df_before['region'] = 'uniform'
    
    # Combine both regions
    df_combined = pd.concat([
        df_before[['bin_min', 'bin_max', 'center', 'width', 'nsig', 'nsig_err', 'region']],
        df_rebinned_after
    ], ignore_index=True)
    
    print("\n" + "="*60)
    print(f"Combined bins: {len(df_combined)}")
    print(f"  Before 0.5: {len(df_before)} bins")
    print(f"  After 0.5: {len(df_rebinned_after)} bins")
    
    # Save combined data
    output_file = "./test_images/nsig_results_exp_rebinned.txt"
    with open(output_file, 'w') as f:
        f.write("# center  width  nsig  nsig_err  nsig_err\n")
        for _, row in df_combined.iterrows():
            f.write(f"{row['center']:.4f}  {row['width']:.4f}  "
                   f"{row['nsig']:.4f}  {row['nsig_err']:.4f}  {row['nsig_err']:.4f}\n")
    print(f"\nSaved rebinned data to: {output_file}")
    
    # Plot comparison
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Original data (linear)
    ax1 = axes[0, 0]
    ax1.errorbar(df['center'], df['nsig'], yerr=df['nsig_err'],
                fmt='o-', markersize=3, capsize=2, alpha=0.7, label='Original')
    ax1.set_xlabel('p_T [GeV/c]', fontsize=12)
    ax1.set_ylabel('$N_{sig}$', fontsize=12)
    ax1.set_title('Original Binning (Linear)', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot 2: Original data (log-y)
    ax2 = axes[0, 1]
    ax2.errorbar(df['center'], df['nsig'], yerr=df['nsig_err'],
                fmt='o-', markersize=3, capsize=2, alpha=0.7, label='Original')
    ax2.set_xlabel('p_T [GeV/c]', fontsize=12)
    ax2.set_ylabel('$N_{sig}$', fontsize=12)
    ax2.set_title('Original Binning (Log)', fontsize=14)
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.legend()
    
    # Plot 3: Combined rebinned data (linear)
    ax3 = axes[1, 0]
    ax3.errorbar(df_combined['center'], df_combined['nsig'], 
                yerr=df_combined['nsig_err'],
                xerr=df_combined['width'],
                fmt='s-', markersize=5, capsize=3, 
                color='red', alpha=0.7, label='Combined')
    ax3.axvline(x=0.5, color='green', linestyle='--', linewidth=2, alpha=0.5, label='Split at 0.5')
    ax3.set_xlabel('p_T [GeV/c]', fontsize=12)
    ax3.set_ylabel('$N_{sig}$', fontsize=12)
    ax3.set_title('Combined: Original (<0.5) + Exp Rebin (≥0.5)', fontsize=14)
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # Plot 4: Combined rebinned data (log-y)
    ax4 = axes[1, 1]
    ax4.errorbar(df_combined['center'], df_combined['nsig'], 
                yerr=df_combined['nsig_err'],
                xerr=df_combined['width'],
                fmt='s-', markersize=5, capsize=3, 
                color='red', alpha=0.7, label='Combined')
    ax4.axvline(x=0.5, color='green', linestyle='--', linewidth=2, alpha=0.5, label='Split at 0.5')
    ax4.set_xlabel('p_T [GeV/c]', fontsize=12)
    ax4.set_ylabel('$N_{sig}$', fontsize=12)
    ax4.set_title('Combined: Original (<0.5) + Exp Rebin (≥0.5) [Log]', fontsize=14)
    ax4.set_yscale('log')
    ax4.grid(True, alpha=0.3, which='both')
    ax4.legend()
    
    plt.tight_layout()
    plot_file = "./test_images/exp_rebinning_comparison.png"
    plt.savefig(plot_file, dpi=300)
    print(f"Plot saved to: {plot_file}")
    plt.show()
    
    # Print statistics
    print("\n" + "="*60)
    print("Final Statistics:")
    print(f"Original total bins: {len(df)}")
    print(f"Combined bins: {len(df_combined)}")
    print(f"  - Uniform bins (x < 0.5): {len(df_before)}")
    print(f"  - Exponential bins (x ≥ 0.5): {len(df_rebinned_after)}")
    print(f"Total events (original): {df['nsig'].sum():.0f}")
    print(f"Total events (combined): {df_combined['nsig'].sum():.0f}")
    print(f"Average events per bin: {df_combined['nsig'].mean():.0f}")
    print(f"Min events: {df_combined['nsig'].min():.0f}")
    print(f"Max events: {df_combined['nsig'].max():.0f}")
    
    # Print bin boundaries for fit_tools.py
    print("\n" + "="*60)
    print("Bin boundaries for fit_tools.py:")
    boundaries_list = list(df_combined['bin_min'].values) + [df_combined.iloc[-1]['bin_max']]
    print(f"boundaries = {[round(b, 4) for b in boundaries_list]}")
    print(f"\nUsage:")
    print(f"bin_var_config = [(\"pt\", {[round(b, 4) for b in boundaries_list]})]")
'''    
