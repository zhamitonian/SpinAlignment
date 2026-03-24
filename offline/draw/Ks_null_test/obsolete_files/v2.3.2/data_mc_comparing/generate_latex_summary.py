#!/usr/bin/env python3
"""
Generate LaTeX document to summarize data/MC comparison plots
Each bin gets one page with:
- 1 fit plot (from fit_result)
- N comparison plots (user-configurable)

Usage:
    # With default phi settings
    generate_latex_document()
    
    # With custom configuration
    binning_config = {"Ks_z": [...], "Ks_pt": [...]}
    var_labels = {"Ks_z": r"z", "Ks_pt": r"p_T"}
    plot_labels = {"Ks_p": r"$p(K_S^0)$", "Ks_costheta": r"$\cos\theta(K_S^0)$", ...}
    generate_latex_document(binning_config, plot_labels, var_labels=var_labels)
"""

import os
import glob
import numpy as np

def generate_latex_document(binning_config=None, plot_labels=None, var_labels=None,
                           title=None, output_file=None, fit_plot_dir=None, 
                           comparison_plot_dir=None, valid_bins=None):
    """
    Generate LaTeX document summarizing all bins
    
    Parameters:
    -----------
    binning_config : dict, optional
        Dictionary mapping variable names to bin edge lists.
        Default: phi particle 3D binning in z, thrust_pt, helicity_angle
    plot_labels : dict, optional
        Dictionary mapping plot names to LaTeX labels for captions.
        Plot names are derived from the keys of this dictionary.
        Default: labels for phi and kaon kinematic plots (9 plots)
    var_labels : dict, optional
        Dictionary mapping variable names to their LaTeX representations.
        Used for displaying bin definitions nicely.
        Default: uses variable names as-is (e.g., 'phi_z' for phi particles)
    title : str, optional
        Document title. Default: "Data/MC Comparison for $\\phi$ Spin Alignment Analysis"
    output_file : str, optional
        Output .tex filename. Default: "data_mc_comparison_summary.tex"
    fit_plot_dir : str or list, optional
        Directory or list of directories to search for fit plots
        Default: ["./fit_result/fit_plots", "./fit_result", "./fit_result/plots"]
    comparison_plot_dir : str, optional
        Directory containing comparison plot subdirectories (bin_XXX)
        Default: "./batch_comparison_images"
    valid_bins : list, optional
        List of bin indices to process. If None, processes all bins.
        Default: None (process all bins that have plots)
    """
    
    # Default binning configuration for phi
    if binning_config is None:
        binning_config = {
            "phi_z": np.linspace(0.2, 1.0, 11).tolist(),
            "phi_thrust_pt": [0.0, 0.125, 0.25, 0.375, 0.5, 0.6611, 0.8688, 1.1366, 1.4817, 1.9265, 2.5],
            "phi_helicity_angle": np.linspace(-1.0, 1.0, 11).tolist()
        }
    
    # Default plot labels (excluding phi_M)
    if plot_labels is None:
        plot_labels = {
            "phi_p": r"$p(\phi)$",
            "phi_costheta": r"$\cos\theta(\phi)$",
            "phi_phi": r"$\varphi(\phi)$",
            "kp_p": r"$p(K^+)$",
            "km_p": r"$p(K^-)$",
            "kp_costheta": r"$\cos\theta(K^+)$",
            "km_costheta": r"$\cos\theta(K^-)$",
            "kp_phi": r"$\varphi(K^+)$",
            "km_phi": r"$\varphi(K^-)$"
        }
    
    # Extract comparison plots from plot_labels keys
    comparison_plots = list(plot_labels.keys())
    
    # Default variable labels for LaTeX presentation
    if var_labels is None:
        var_labels = {
            "phi_z": r"z",
            "phi_thrust_pt": r"p_T^{\rm thrust}",
            "phi_helicity_angle": r"\cos\theta_{\rm hel}"
        }
    
    # Default title
    if title is None:
        title = r"Data/MC Comparison for $\phi$ Spin Alignment Analysis"
    
    # Default output filename
    if output_file is None:
        output_file = "data_mc_comparison_summary.tex"
    
    # Default fit plot directories (will search in order)
    if fit_plot_dir is None:
        fit_plot_dir = ["./fit_result/fit_plots", "./fit_result", "./fit_result/plots"]
    elif isinstance(fit_plot_dir, str):
        fit_plot_dir = [fit_plot_dir]
    
    # Default comparison plot directory
    if comparison_plot_dir is None:
        comparison_plot_dir = "./batch_comparison_images"
    
    bins = [len(binning_list) - 1 for binning_list in binning_config.values()]
    total_bins = np.prod(bins)
    pad_width = len(str(total_bins - 1))
    
    # Extract variable names for bin labels
    var_names = list(binning_config.keys())
    
    latex_content = []
    
    # LaTeX header
    latex_content.append(r"\documentclass[11pt,a4paper]{article}")
    latex_content.append(r"\usepackage{graphicx}")
    latex_content.append(r"\usepackage{geometry}")
    latex_content.append(r"\usepackage{subcaption}")
    latex_content.append(r"\usepackage{amsmath}")
    latex_content.append(r"\geometry{margin=1cm}")
    latex_content.append(r"")
    latex_content.append(r"\begin{document}")
    latex_content.append(r"")
    latex_content.append(r"\title{" + title + r"}")
    
    # Create subtitle from binning dimensions - use var_labels for nice display
    n_dims = len(var_names)
    # Use var_labels if available, wrap in $ for math mode
    display_vars = []
    for var_name in var_names:
        if var_name in var_labels:
            # Strip existing $ and re-wrap to ensure math mode
            clean_label = var_labels[var_name].strip('$')
            display_vars.append(f"${clean_label}$")
        else:
            # Escape underscores for non-math variable names
            display_vars.append(var_name.replace('_', r'\_'))
    subtitle = f"{n_dims}D Binning in " + ", ".join(display_vars)
    latex_content.append(r"\author{" + subtitle + r"}")
    latex_content.append(r"\date{\today}")
    latex_content.append(r"\maketitle")
    latex_content.append(r"")
    latex_content.append(r"\tableofcontents")
    latex_content.append(r"\clearpage")
    latex_content.append(r"")
    
    processed_bins = 0
    skipped_bins = []

    # Loop over all bins
    for flat_bin_idx in range(total_bins):
        
        # Skip bins not in valid_bins list if provided
        if valid_bins is not None and flat_bin_idx not in valid_bins:
            skipped_bins.append(flat_bin_idx)
            continue
        
        # Calculate N-dimensional bin indices
        remaining = flat_bin_idx
        bin_idx = []
        for i in range(len(bins)):
            divisor = int(np.prod(bins[i+1:])) if i < len(bins) - 1 else 1
            bin_idx.append(remaining // divisor)
            remaining = remaining % divisor
        
        # Get bin boundaries for each dimension
        bin_ranges = {}
        for i, var_name in enumerate(var_names):
            idx = bin_idx[i]
            bin_ranges[var_name] = (binning_config[var_name][idx], binning_config[var_name][idx + 1])
        
        # Check if bin directory exists
        bin_dir = os.path.join(comparison_plot_dir, f"bin_{flat_bin_idx:0{pad_width}d}")
        
        # Search for fit plot in multiple locations
        fit_plot = None
        for fit_dir in fit_plot_dir:
            # Try specific fit plot naming patterns
            potential_paths = [
                os.path.join(fit_dir, f"bin_{flat_bin_idx:0{pad_width}d}_fit.png"),
                os.path.join(fit_dir, f"bin_{flat_bin_idx:0{pad_width}d}.png")
            ]
            for path in potential_paths:
                if os.path.exists(path):
                    fit_plot = path
                    break
            
            # If not found, try wildcard pattern
            if not fit_plot:
                pattern = os.path.join(fit_dir, f"bin_{flat_bin_idx:0{pad_width}d}_*.png")
                matches = glob.glob(pattern)
                if matches:
                    fit_plot = matches[0]  # Use first match
                    break
            
            if fit_plot:
                break
        
        # Check if at least some comparison plots exist
        has_plots = os.path.exists(bin_dir) and any(
            os.path.exists(os.path.join(bin_dir, f"{plot}.png")) 
            for plot in comparison_plots
        )
        
        if not has_plots:
            skipped_bins.append(flat_bin_idx)
            continue  # Skip bins without plots

        processed_bins += 1
        
        # Start new page for this bin
        latex_content.append(r"\clearpage")
        latex_content.append(r"\section*{Bin " + f"{flat_bin_idx:0{pad_width}d}" + r"}")
        latex_content.append(r"")
        
        # Bin information and fit plot side by side
        latex_content.append(r"\subsection*{Bin Definition and Mass Fit}")
        
        if fit_plot and os.path.exists(fit_plot):
            # Two-column layout: bin definition on left, fit plot on right
            latex_content.append(r"\begin{minipage}[t]{0.48\textwidth}")
            latex_content.append(r"\textbf{Bin Definition:}")
            latex_content.append(r"\begin{itemize}")
        else:
            # No fit plot, use full width for bin definition
            latex_content.append(r"\textbf{Bin Definition:}")
            latex_content.append(r"\begin{itemize}")
        
        # Add bin ranges for each variable
        for var_name in var_names:
            low, high = bin_ranges[var_name]
            # Use var_labels if provided, otherwise use variable name as-is
            display_name = var_labels.get(var_name, var_name)
            # Strip any existing $ signs from display_name to avoid double-wrapping
            display_name_clean = display_name.strip('$')
            latex_content.append(f"  \\item ${low:.3f} < {display_name_clean} < {high:.4f}$")
        
        # Add bin indices
        bin_idx_parts = []
        for i in range(len(var_names)):
            var_display = var_labels.get(var_names[i], var_names[i]).strip('$')
            bin_idx_parts.append(f"{var_display}_{{\\rm bin}}={bin_idx[i]}")
        bin_idx_str = ", ".join(bin_idx_parts)
        latex_content.append(f"  \\item Bin indices: ${bin_idx_str}$")
        latex_content.append(r"\end{itemize}")
        
        # Fit plot on the right side (if exists)
        if fit_plot and os.path.exists(fit_plot):
            latex_content.append(r"\end{minipage}")
            latex_content.append(r"\hfill")
            latex_content.append(r"\begin{minipage}[t]{0.48\textwidth}")
            latex_content.append(r"\textbf{Mass Fit:}")
            latex_content.append(r"\vspace{0.3cm}")
            latex_content.append(r"\centering")
            latex_content.append(r"\includegraphics[width=\textwidth]{" + fit_plot + r"}")
            latex_content.append(r"\end{minipage}")
        
        latex_content.append(r"")
        
        # Comparison plots
        latex_content.append(r"\subsection*{Data/MC Comparisons}")
        latex_content.append(r"\begin{figure}[h!]")
        latex_content.append(r"  \centering")
        
        # Determine grid layout (3 plots per row)
        plots_per_row = 3
        width_per_plot = 0.32  # 32% width for 3 plots per row
        
        for i, plot_name in enumerate(comparison_plots):
            plot_file = os.path.join(bin_dir, f"{plot_name}.png")
            
            if os.path.exists(plot_file):
                latex_content.append(f"  \\begin{{subfigure}}[b]{{{width_per_plot}\\textwidth}}")
                latex_content.append(r"    \centering")
                latex_content.append(r"    \includegraphics[width=\textwidth]{" + plot_file + r"}")
                latex_content.append(r"    \caption{" + plot_labels.get(plot_name, plot_name) + r"}")
                latex_content.append(r"  \end{subfigure}")
                if (i + 1) % plots_per_row == 0 and i < len(comparison_plots) - 1:
                    latex_content.append(r"  \\")  # New row after every N plots
            else:
                # Empty placeholder if plot doesn't exist
                latex_content.append(f"  \\begin{{subfigure}}[b]{{{width_per_plot}\\textwidth}}")
                latex_content.append(r"    \centering")
                latex_content.append(r"    \vspace{3cm}")
                latex_content.append(r"    \caption{N/A}")
                latex_content.append(r"  \end{subfigure}")
                if (i + 1) % plots_per_row == 0 and i < len(comparison_plots) - 1:
                    latex_content.append(r"  \\")
        
        latex_content.append(r"  \caption{Data/MC comparison plots for bin " + f"{flat_bin_idx:0{pad_width}d}" + r"}")
        latex_content.append(r"\end{figure}")
        latex_content.append(r"")
    
    # End document
    latex_content.append(r"\clearpage")
    latex_content.append(r"\section*{Summary}")
    latex_content.append(f"Total bins processed: {processed_bins} out of {total_bins}")
    if skipped_bins:
        latex_content.append(f"\\\\Bins skipped: {len(skipped_bins)}")
    latex_content.append(r"")
    latex_content.append(r"\end{document}")
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(latex_content))
    
    print(f"LaTeX document generated: {output_file}")
    print(f"Bins processed: {processed_bins}/{total_bins}")
    if skipped_bins:
        print(f"Bins skipped: {len(skipped_bins)}")
    print(f"\nTo compile:")
    print(f"  pdflatex {output_file}")
    print(f"  pdflatex {output_file}  # Run twice for TOC")
    
    return output_file


if __name__ == "__main__":
    # Example usage with default phi configuration
    generate_latex_document()
    
    # Example for Ks:
    # binning_config = {
    #     "Ks_z": np.linspace(0.2, 1.0, 11).tolist(),
    #     "Ks_pt": [0.0, 0.5, 1.0, 1.5, 2.0, 2.5],
    #     "Ks_helicity_angle": np.linspace(-1.0, 1.0, 11).tolist()
    # }
    # var_labels = {
    #     "Ks_z": r"z",
    #     "Ks_pt": r"p_T",
    #     "Ks_helicity_angle": r"\cos\theta_{\rm hel}"
    # }
    # plot_labels = {
    #     "Ks_p": r"$p(K_S^0)$",
    #     "Ks_costheta": r"$\cos\theta(K_S^0)$",
    #     "Ks_phi": r"$\varphi(K_S^0)$",
    #     "pip_p": r"$p(\pi^+)$",
    #     "pim_p": r"$p(\pi^-)$"
    # }
    # generate_latex_document(
    #     binning_config=binning_config,
    #     plot_labels=plot_labels,
    #     var_labels=var_labels,
    #     title=r"Data/MC Comparison for $K_S^0$ Spin Alignment Analysis",
    #     output_file="Ks_data_mc_comparison.tex"
    # )

