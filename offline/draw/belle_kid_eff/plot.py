import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#plt.rcParams['font.family'] = 'Times New Roman'

def plot_pid_eff(period, in_barrel_only=False, is_pion=True, cut_value=6, output_file='pid_efficiency.png'):
    plab_boundaries = np.array([
        0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
        2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 4.0, 4.5, 100.0
    ])
    max_plot_p = 5.0
    plab_boundaries[-1] = max_plot_p #for plotting purposes only

    costheta_boundaries = np.array([
        -1, -0.612, -0.511, -0.300, -0.152, 0.017,
        0.209, 0.355, 0.435, 0.542, 0.692, 0.842, 1
    ])

    plab_widths = plab_boundaries[1:] - plab_boundaries[:-1]
    costheta_widths = costheta_boundaries[1:] - costheta_boundaries[:-1]

    def load_kid_eff(filename):
        df = pd.read_csv(filename, sep=r'\s+', comment='#',
                names=['kind', 'pid', 'map', 'eff_dt', 'statErr_dt', 'systErr_dt', 
                       'eff_mc', 'statErr_mc', 'ratio', 'statErr_ratio', 'systErr_ratio' ,'flag' ]) 
        if is_pion:
            df = df[(df['kind']==2) & (df['pid'] == cut_value) & (df['flag']==0) ]
        else:
            df = df[(df['kind']==0) & (df['pid'] == cut_value) & (df['flag']==0) ]
        df['map'] = df['map'].astype(str).str.zfill(4)
        df['pbin'] = df['map'].str[:2].astype(int) - 1
        df['ctbin'] = df['map'].str[2:].astype(int) - 1
        if in_barrel_only:
            df = df[(df['ctbin']>=2) & (df['ctbin'] <= len(costheta_boundaries)-2)]
        return df

    if period == "svd1":
        df = load_kid_eff('kideff-2006-svd1-all.dat')
    elif period == "svd2":
        df = load_kid_eff('kideff-2010.dat')
    else: # combined
        df_svd1 = load_kid_eff('kideff-2006-svd1-all.dat')
        df_svd2 = load_kid_eff('kideff-2010.dat')
        df = pd.concat([df_svd1, df_svd2], ignore_index=True)

    def combine_bins(effs, stat_errs, syst_err=None):
        if syst_err is None:
            syst_err = np.array([0.0]*len(effs))

        mask = stat_errs > 0
        if not np.any(mask):
            return np.nan, np.nan
        effs = effs[mask]
        stat_errs = stat_errs[mask]
        syst_err = syst_err[mask]

        #weights = 1.0 / (stat_errs**2 + syst_err**2)
        weights = 1.0 / stat_errs**2
        #eff_combined = np.sum(effs * weights) / np.sum(weights) # same as below
        eff_combined = np.average(effs, weights=weights)
        stat_err_combined = np.sqrt(1.0 / np.sum(weights))
        syst_err_combined = np.average(syst_err, weights=weights)
        total_err_combined = np.sqrt(stat_err_combined**2 + syst_err_combined**2)

        return eff_combined, total_err_combined

    grouped = df.groupby('pbin')
    eff_mc, err_mc = zip(*[combine_bins(g['eff_mc'].values, g['statErr_mc'].values ) for _, g in grouped])
    eff_mc = np.array(eff_mc)
    err_mc = np.array(err_mc)

    eff_data, err_data = zip(*[combine_bins(g['eff_dt'].values, g['statErr_dt'].values, g['systErr_dt'].values) for _, g in grouped])
    #eff_ratio, err_ratio = zip(*[combine_bins(g['ratio'].values, g['statErr_ratio'].values, g['systErr_ratio'].values) for _, g in grouped])
    eff_data = np.array(eff_data)
    err_data = np.array(err_data)
    
    # Calculate ratio from combined efficiencies (more accurate than averaging ratios)
    eff_ratio = eff_data / eff_mc
    # Error propagation: sigma_r = r * sqrt((sigma_d/d)^2 + (sigma_m/m)^2)
    err_ratio = eff_ratio * np.sqrt((err_data/eff_data)**2 + (err_mc/eff_mc)**2)

    p_centers = 0.5 * (plab_boundaries[:-1] + plab_boundaries[1:])

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, 
                gridspec_kw={'height_ratios': [3, 1]}, figsize=(8, 6) )

    ax1.errorbar(p_centers, eff_data, yerr=err_data, fmt='o', label='Data', color='blue')
    ax1.errorbar(p_centers, eff_mc, yerr=err_mc, fmt='o', label='MC', color='orange')

    label = 'PionID' if is_pion else 'KaonID'
    ax1.set_ylabel(f'{label} Efficieny')
    label = 'SVD1 ' if period=='svd1' else ('SVD2 ' if period=='svd2' else '')  + label
    ax1.set_title(f'{label} Efficiency vs. p')
    if in_barrel_only:
        ax1.legend(title=r'$-0.511 < \cos\theta < 0.842$')
    ax1.legend()
    ax1.set_ylim(0.6, 1.1)
    ax1.grid()

    from matplotlib.ticker import MultipleLocator
    ax2.errorbar(p_centers, eff_ratio, yerr=err_ratio, fmt='o', color='green')
    ax2.axhline(1, color='gray', linestyle='--')
    ax2.set_xlabel('Momentum (GeV/c)')
    ax2.set_ylabel(r'$\varepsilon_{Data}/\varepsilon_{MC}$')
    ax2.set_ylim(0.9, 1.1)
    ax2.set_xlim(0, max_plot_p)
    ax2.yaxis.set_minor_locator(MultipleLocator(0.02)) 
    for y in np.arange(0.9, 1.1+0.02, 0.02):
        ax2.axhline(y, color='gray', linestyle=':', linewidth=0.5, zorder=0)
    ax2.grid()

    ax1.xaxis.set_minor_locator(MultipleLocator(0.2))
    ax2.xaxis.set_minor_locator(MultipleLocator(0.2))
    ax1.xaxis.set_major_locator(MultipleLocator(1.0))
    ax2.xaxis.set_major_locator(MultipleLocator(1.0))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax1.tick_params(axis='x', which='minor', length=4, direction='in', top=True, bottom=True)
    ax2.tick_params(axis='x', which='minor', length=4, direction='in', top=True, bottom=True)
    ax1.tick_params(axis='x', which='major', length=8, direction='in', top=True, bottom=True)
    ax2.tick_params(axis='x', which='major', length=8, direction='in', top=True, bottom=True)
    ax1.tick_params(axis='y', which='minor', length=4, direction='in', top=True, bottom=True)
    ax1.tick_params(axis='y', which='major', length=8, direction='in', top=True, bottom=True)


    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

    return

def eff_proj_2d(period, is_pion=True, cut_value=6, output_file='ratio_2d.png'):
    import matplotlib.patheffects as patheffects
    # Load data
    plab_boundaries = np.array([
        0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
        2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 4.0, 4.5, 100.0
    ])
    max_plot_p = 5.0
    plab_boundaries[-1] = max_plot_p
    costheta_boundaries = np.array([
        -1, -0.612, -0.511, -0.300, -0.152, 0.017,
        0.209, 0.355, 0.435, 0.542, 0.692, 0.842, 1
    ])

    def load_kid_eff(filename):
        df = pd.read_csv(filename, sep=r'\s+', comment='#',
                names=['kind', 'pid', 'map', 'eff_dt', 'statErr_dt', 'systErr_dt',
                       'eff_mc', 'statErr_mc', 'ratio', 'statErr_ratio', 'systErr_ratio' ,'flag' ])
        if is_pion:
            df = df[(df['kind']==2) & (df['pid'] == cut_value) & (df['flag']==0) ]
        else:
            df = df[(df['kind']==0) & (df['pid'] == cut_value) & (df['flag']==0) ]
        df['map'] = df['map'].astype(str).str.zfill(4)
        df['pbin'] = df['map'].str[:2].astype(int) - 1
        df['ctbin'] = df['map'].str[2:].astype(int) - 1
        return df

    if period == "svd1":
        df = load_kid_eff('kideff-2006-svd1-all.dat')
    elif period == "svd2":
        df = load_kid_eff('kideff-2010.dat')
    else: # combined
        df_svd1 = load_kid_eff('kideff-2006-svd1-all.dat')
        df_svd2 = load_kid_eff('kideff-2010.dat')
        df = pd.concat([df_svd1, df_svd2], ignore_index=True)

    n_p = len(plab_boundaries) - 1
    n_ct = len(costheta_boundaries) - 1

    # Prepare 2D arrays for efficiency
    eff_data_2d = np.full((n_p, n_ct), np.nan)
    eff_mc_2d = np.full((n_p, n_ct), np.nan)
    eff_ratio_2d = np.full((n_p, n_ct), np.nan)

    for _, row in df.iterrows():
        pbin = int(row['pbin'])
        ctbin = int(row['ctbin'])
        if 0 <= pbin < n_p and 0 <= ctbin < n_ct:
            eff_data_2d[pbin, ctbin] = row['eff_dt']
            eff_mc_2d[pbin, ctbin] = row['eff_mc']
            eff_ratio_2d[pbin, ctbin] = row['ratio']

    # Plot 2D efficiency maps
    # --- Data/MC Ratio ---
    fig2, ax2 = plt.subplots(figsize=(20, 15), dpi=200)
    # pcolormesh expects 2D arrays for X, Y bin edges
    X, Y = np.meshgrid(costheta_boundaries, plab_boundaries)
    pcm = ax2.pcolormesh(X, Y, eff_ratio_2d, vmin=0.8, vmax=1.2, cmap='coolwarm', shading='auto')
    label =  'SVD1, ' if period=='svd1' else ('SVD2, ' if period=='svd2' else '')  
    label += f'PionID > 0.{cut_value}' if is_pion else f'KaonID > 0.{cut_value}'

    ax2.set_title(r'$\varepsilon(Data)/\varepsilon(MC)$ distribution ' + f' ({label})')

    fig2.colorbar(pcm, ax=ax2, fraction=0.046, pad=0.04)
    # Annotate at bin centers
    for i in range(eff_ratio_2d.shape[0]):
        for j in range(eff_ratio_2d.shape[1]):
            val = eff_ratio_2d[i, j]
            if not np.isnan(val):
                p_center = 0.5 * (plab_boundaries[i] + plab_boundaries[i+1])
                ct_center = 0.5 * (costheta_boundaries[j] + costheta_boundaries[j+1])
                ax2.text(ct_center, p_center, '{:.3f}'.format(val),
                        ha='center', va='center', color='white', fontsize=8, fontweight='bold',
                        path_effects=[patheffects.withStroke(linewidth=1, foreground='black')])
    # 设置x轴和y轴主刻度为bin边界
    ax2.set_xticks(costheta_boundaries)
    ax2.set_xticklabels([f"{x:.3f}" for x in costheta_boundaries], rotation=45, ha='right')
    ax2.set_yticks(plab_boundaries)
    ax2.set_yticklabels([f"{y:.1f}" for y in plab_boundaries[0:-1]] + [r'$\infty$'])
    ax2.set_xlabel(r'$\cos\theta$')
    ax2.set_xlim(costheta_boundaries[0], costheta_boundaries[-1])
    ax2.set_ylabel('Momentum (GeV/c)')
    ax2.set_ylim(plab_boundaries[0], plab_boundaries[-1])
    ax2.grid(True, which='both', linestyle=':', linewidth=0.7)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    #plt.show()

    # ...existing code...

if __name__ == "__main__":
    plot_pid_eff(period='svd2', in_barrel_only=True, 
                 is_pion=True, cut_value=6, output_file="images/PiID_eff_svd2_6_barrel.png")

    """
    for is_pion in [True, False]:
        for cut_value in [6, 8]:
            for period in ['svd1', 'svd2']:
                output_file = f"images/{'PiID' if is_pion else 'KID'}_ratio_2d_{period}_{cut_value}.png"
                eff_proj_2d(period=period, is_pion=is_pion, cut_value=cut_value, output_file=output_file)
    """        