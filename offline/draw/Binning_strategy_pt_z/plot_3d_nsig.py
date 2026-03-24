#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def plot_3d_nsig(txt_file, output_dir="./"):
    """Plot 3D nsig data with multiple visualizations"""
    
    # Read 3D binning data
    df = pd.read_csv(txt_file, sep='\s+', comment='#',
                    names=['z_center', 'z_width', 'pt_center', 'pt_width', 
                          'cos_center', 'cos_width', 'nsig', 'nsig_err', 'nsig_err2'])
    
    print(f"Loaded {len(df)} bins (before filtering)")
    
    # Remove bins with nsig = 0
    df = df[df['nsig'] > 0].reset_index(drop=True)
    
    print(f"After removing zero bins: {len(df)} bins")
    print(f"Data shape: {df.shape}")
    print(f"\nData ranges:")
    print(f"  z: [{df['z_center'].min():.3f}, {df['z_center'].max():.3f}]")
    print(f"  pt: [{df['pt_center'].min():.3f}, {df['pt_center'].max():.3f}]")
    print(f"  cos(θ): [{df['cos_center'].min():.3f}, {df['cos_center'].max():.3f}]")
    print(f"  N_sig: [{df['nsig'].min():.0f}, {df['nsig'].max():.0f}]")
    
    # ============ Plot 1: 3D Scatter Plot ============
    fig1 = plt.figure(figsize=(10, 8))
    ax1 = fig1.add_subplot(111, projection='3d')
    scatter = ax1.scatter(df['z_center'], df['pt_center'], df['cos_center'],
                         c=df['nsig'], cmap='viridis', s=20, alpha=0.6)
    ax1.set_xlabel('z', fontsize=14)
    ax1.set_ylabel('$p_T$ [GeV/c]', fontsize=14)
    ax1.set_zlabel(r'cos($\theta$)', fontsize=14)
    ax1.set_title('3D Scatter: N_sig', fontsize=16)
    plt.colorbar(scatter, ax=ax1, label='$N_{sig}$', shrink=0.7)
    plt.tight_layout()
    fig1.savefig(f"{output_dir}/plot1_3d_scatter.png", dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/plot1_3d_scatter.png")
    plt.close(fig1)
    
    # ============ Plot 2: 3D Scatter (log scale) ============
    fig2 = plt.figure(figsize=(10, 8))
    ax2 = fig2.add_subplot(111, projection='3d')
    scatter2 = ax2.scatter(df['z_center'], df['pt_center'], df['cos_center'],
                          c=np.log10(df['nsig']+1), cmap='plasma', s=20, alpha=0.6)
    ax2.set_xlabel('z', fontsize=14)
    ax2.set_ylabel('$p_T$ [GeV/c]', fontsize=14)
    ax2.set_zlabel(r'cos($\theta$)', fontsize=14)
    ax2.set_title('3D Scatter: log(N_sig)', fontsize=16)
    plt.colorbar(scatter2, ax=ax2, label='log$_{10}(N_{sig})$', shrink=0.7)
    plt.tight_layout()
    fig2.savefig(f"{output_dir}/plot2_3d_scatter_log.png", dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/plot2_3d_scatter_log.png")
    plt.close(fig2)
    
    # ============ Plot 3: Projection onto z-pt plane ============
    fig3 = plt.figure(figsize=(10, 8))
    ax3 = fig3.add_subplot(111)
    # Sum over cos(theta) dimension
    df_proj_zpt = df.groupby(['z_center', 'pt_center'])['nsig'].sum().reset_index()
    
    # Filter z range [0.2, 1.0]
    df_proj_zpt = df_proj_zpt[(df_proj_zpt['z_center'] >= 0.2) & (df_proj_zpt['z_center'] <= 1.0)]
    
    z_unique = sorted(df_proj_zpt['z_center'].unique())
    pt_unique = sorted(df_proj_zpt['pt_center'].unique())
    Z, PT = np.meshgrid(z_unique, pt_unique)
    Nsig_zpt = np.zeros_like(Z)
    
    for i, pt in enumerate(pt_unique):
        for j, z in enumerate(z_unique):
            val = df_proj_zpt[(df_proj_zpt['z_center']==z) & (df_proj_zpt['pt_center']==pt)]['nsig']
            Nsig_zpt[i, j] = val.values[0] if len(val) > 0 else 0
    
    im3 = ax3.pcolormesh(Z, PT, Nsig_zpt, cmap='viridis', shading='auto')
    
    # Mark bins with N_sig < 1000
    low_stat_bins_zpt = df_proj_zpt[df_proj_zpt['nsig'] < 1000]
    if len(low_stat_bins_zpt) > 0:
        from matplotlib.patches import Rectangle
        
        for _, row in low_stat_bins_zpt.iterrows():
            z_center = row['z_center']
            pt_center = row['pt_center']
            
            # Get actual bin widths for this specific bin from original data
            bin_info = df[(df['z_center'] == z_center) & (df['pt_center'] == pt_center)]
            if len(bin_info) > 0:
                z_width = bin_info['z_width'].iloc[0] * 2  # Full bin width
                pt_width = bin_info['pt_width'].iloc[0] * 2  # Full bin width (variable)
                
                # Draw rectangle around low stat bin
                rect = Rectangle((z_center - z_width/2, pt_center - pt_width/2),
                               z_width, pt_width,
                               linewidth=2, edgecolor='red', facecolor='none',
                               linestyle='--', alpha=0.8)
                ax3.add_patch(rect)
    
    ax3.set_xlabel('z', fontsize=14)
    ax3.set_ylabel('$p_T$ [GeV/c]', fontsize=14)
    ax3.set_title('Projection: z vs $p_T$ (red box: $\Sigma N_{sig} < 1000$)', fontsize=16)
    plt.colorbar(im3, ax=ax3, label='$\Sigma N_{sig}$')
    plt.tight_layout()
    fig3.savefig(f"{output_dir}/plot3_proj_z_pt.png", dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/plot3_proj_z_pt.png")
    plt.close(fig3)
    
    # ============ Plot 4: Projection onto pt-cos plane ============
    fig4 = plt.figure(figsize=(10, 8))
    ax4 = fig4.add_subplot(111)
    df_proj_ptcos = df.groupby(['pt_center', 'cos_center'])['nsig'].sum().reset_index()
    
    cos_unique = sorted(df_proj_ptcos['cos_center'].unique())
    PT2, COS = np.meshgrid(pt_unique, cos_unique)
    Nsig_ptcos = np.zeros_like(PT2)
    
    for i, cos in enumerate(cos_unique):
        for j, pt in enumerate(pt_unique):
            val = df_proj_ptcos[(df_proj_ptcos['pt_center']==pt) & (df_proj_ptcos['cos_center']==cos)]['nsig']
            Nsig_ptcos[i, j] = val.values[0] if len(val) > 0 else 0
    
    im4 = ax4.pcolormesh(PT2, COS, Nsig_ptcos, cmap='plasma', shading='auto')
    ax4.set_xlabel('$p_T$ [GeV/c]', fontsize=14)
    ax4.set_ylabel(r'cos($\theta$)', fontsize=14)
    ax4.set_title('Projection: $p_T$ vs cos($\\theta$)', fontsize=16)
    plt.colorbar(im4, ax=ax4, label='$\Sigma N_{sig}$')
    plt.tight_layout()
    fig4.savefig(f"{output_dir}/plot4_proj_pt_cos.png", dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/plot4_proj_pt_cos.png")
    plt.close(fig4)
    
    # ============ Plot 5: Projection onto z-cos plane ============
    fig5 = plt.figure(figsize=(10, 8))
    ax5 = fig5.add_subplot(111)
    df_proj_zcos = df.groupby(['z_center', 'cos_center'])['nsig'].sum().reset_index()
    
    # Filter z range [0.2, 1.0]
    df_proj_zcos = df_proj_zcos[(df_proj_zcos['z_center'] >= 0.2) & (df_proj_zcos['z_center'] <= 1.0)]
    z_unique_filtered = sorted(df_proj_zcos['z_center'].unique())
    
    Z2, COS2 = np.meshgrid(z_unique_filtered, cos_unique)
    Nsig_zcos = np.zeros_like(Z2)
    
    for i, cos in enumerate(cos_unique):
        for j, z in enumerate(z_unique_filtered):
            val = df_proj_zcos[(df_proj_zcos['z_center']==z) & (df_proj_zcos['cos_center']==cos)]['nsig']
            Nsig_zcos[i, j] = val.values[0] if len(val) > 0 else 0
    
    im5 = ax5.pcolormesh(Z2, COS2, Nsig_zcos, cmap='coolwarm', shading='auto')
    ax5.set_xlabel('z', fontsize=14)
    ax5.set_ylabel(r'cos($\theta$)', fontsize=14)
    ax5.set_title('Projection: z vs cos($\\theta$)', fontsize=16)
    plt.colorbar(im5, ax=ax5, label='$\Sigma N_{sig}$')
    plt.tight_layout()
    fig5.savefig(f"{output_dir}/plot5_proj_z_cos.png", dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/plot5_proj_z_cos.png")
    plt.close(fig5)
    
    # ============ Plot 6: 1D projections ============
    fig6 = plt.figure(figsize=(10, 8))
    ax6 = fig6.add_subplot(111)
    
    # Filter z range for 1D projections
    df_filtered_z = df[(df['z_center'] >= 0.2) & (df['z_center'] <= 1.0)]
    
    # Sum over other dimensions
    nsig_vs_z = df_filtered_z.groupby('z_center')['nsig'].sum()
    nsig_vs_pt = df.groupby('pt_center')['nsig'].sum()
    nsig_vs_cos = df.groupby('cos_center')['nsig'].sum()
    
    ax6_z = ax6.twinx()
    ax6_pt = ax6.twinx()
    
    # Offset the right spine
    ax6_pt.spines['right'].set_position(('outward', 60))
    
    p1, = ax6.plot(nsig_vs_z.index, nsig_vs_z.values, 'b-o', markersize=4, label='z')
    p2, = ax6_z.plot(nsig_vs_pt.index, nsig_vs_pt.values, 'r-s', markersize=4, label='$p_T$')
    p3, = ax6_pt.plot(nsig_vs_cos.index, nsig_vs_cos.values, 'g-^', markersize=4, label=r'cos($\theta$)')
    
    ax6.set_xlabel('Variable Value', fontsize=12)
    ax6.set_ylabel('$\Sigma N_{sig}$ (z)', color='b', fontsize=12)
    ax6_z.set_ylabel('$\Sigma N_{sig}$ ($p_T$)', color='r', fontsize=12)
    ax6_pt.set_ylabel('$\Sigma N_{sig}$ (cos)', color='g', fontsize=12)
    
    ax6.tick_params(axis='y', labelcolor='b')
    ax6_z.tick_params(axis='y', labelcolor='r')
    ax6_pt.tick_params(axis='y', labelcolor='g')
    
    ax6.set_title('1D Projections', fontsize=14)
    ax6.legend(handles=[p1, p2, p3], loc='upper right')
    
    plt.tight_layout()
    output_file = f"{output_dir}/nsig_3d_visualization.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")
    plt.show()
    
    # Print statistics
    print("\n" + "="*60)
    print("Statistics:")
    print(f"Total events: {df['nsig'].sum():.0f}")
    print(f"Mean per bin: {df['nsig'].mean():.0f}")
    print(f"Std per bin: {df['nsig'].std():.0f}")
    
    # Find bins with low statistics
    print("\n" + "="*60)
    print("Bins with N_sig < 100:")
    low_stat_bins = df[df['nsig'] < 100]
    print(f"Number of bins: {len(low_stat_bins)}")
    
    if len(low_stat_bins) > 0:
        print("\nDetails:")
        for idx, row in low_stat_bins.iterrows():
            print(f"  z={row['z_center']:.3f}, pt={row['pt_center']:.3f}, "
                  f"cos={row['cos_center']:.3f}: N_sig={row['nsig']:.1f} ± {row['nsig_err']:.1f}")
        
        # Save to file
        low_stat_file = f"{output_dir}/low_stat_bins.txt"
        with open(low_stat_file, 'w') as f:
            f.write("# Bins with N_sig < 100\n")
            f.write("# z_center  pt_center  cos_center  nsig  nsig_err\n")
            for idx, row in low_stat_bins.iterrows():
                f.write(f"{row['z_center']:.4f}  {row['pt_center']:.4f}  "
                       f"{row['cos_center']:.4f}  {row['nsig']:.2f}  {row['nsig_err']:.2f}\n")
        print(f"\nLow statistics bins saved to: {low_stat_file}")
    else:
        print("  None found!")
    
    # Additional thresholds
    for threshold in [50, 200, 500]:
        count = len(df[df['nsig'] < threshold])
        print(f"Bins with N_sig < {threshold}: {count}")
    
    print(f"\nBins distribution:")
    print(f"  N_sig < 100: {len(df[df['nsig'] < 100])}")
    print(f"  100 ≤ N_sig < 500: {len(df[(df['nsig'] >= 100) & (df['nsig'] < 500)])}")
    print(f"  500 ≤ N_sig < 1000: {len(df[(df['nsig'] >= 500) & (df['nsig'] < 1000)])}")
    print(f"  N_sig ≥ 1000: {len(df[df['nsig'] >= 1000])}")


def main():
    import sys
    if len(sys.argv) < 2:
        print("Usage: python plot_3d_nsig.py <nsig_results.txt> [output_dir]")
        return
    
    txt_file = sys.argv[1]
    output_dir = "./"
    if len(sys.argv) >= 3:
        output_dir = sys.argv[2]
    
    plot_3d_nsig(txt_file, output_dir)


if __name__ == "__main__":
    main()