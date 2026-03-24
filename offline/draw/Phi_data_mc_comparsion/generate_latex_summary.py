#!/usr/bin/env python3
"""
Generate LaTeX document to summarize data/MC comparison plots
Each bin gets one page with:
- 1 fit plot (from fit_result)
- 9 comparison plots (excluding phi_M)
"""

import os
import numpy as np

def generate_latex_document():
    """Generate LaTeX document summarizing all bins"""
    
    binning_config = {
        "phi_z": np.linspace(0.2, 1.0, 11).tolist(),
        "phi_thrust_pt": [0.0, 0.125, 0.25, 0.375, 0.5, 0.6611, 0.8688, 1.1366, 1.4817, 1.9265, 2.5],
        "phi_helicity_angle": np.linspace(-1.0, 1.0, 11).tolist()
    }
    
    bins = [len(binning_list) - 1 for binning_list in binning_config.values()]
    total_bins = bins[0] * bins[1] * bins[2]
    pad_width = len(str(total_bins - 1))
     
    # Comparison plots to include (excluding phi_M)
    comparison_plots = [
        "phi_p", "phi_costheta", "phi_phi",
        "kp_p", "km_p", 
        "kp_costheta", "km_costheta",
        "kp_phi", "km_phi"
    ]
    
    # Plot labels for caption
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
    latex_content.append(r"\title{Data/MC Comparison for $\phi$ Spin Alignment Analysis}")
    latex_content.append(r"\author{3D Binning in $z$, $p_T$, and Helicity Angle}")
    latex_content.append(r"\date{\today}")
    latex_content.append(r"\maketitle")
    latex_content.append(r"")
    latex_content.append(r"\tableofcontents")
    latex_content.append(r"\clearpage")
    latex_content.append(r"")
    
    valid_bins = 0
    unvalid_bins = []

    valid_nsig_bins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 413, 414, 415, 416, 417, 418, 419, 420, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 476, 477, 479, 480, 481, 482, 483, 485, 486, 487, 488, 489, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 714, 715, 716, 717, 718, 719, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 800, 801, 802, 803, 804, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 817, 819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 900, 901, 902, 903, 904, 905, 906, 907, 908, 909]


    # Loop over all bins
    for flat_bin_idx in range(total_bins):
        
        if flat_bin_idx not in  valid_nsig_bins:
            unvalid_bins.append(flat_bin_idx)
            continue

    #for flat_bin_idx in range(3):
        # Calculate 3D bin indices
        bin_idx = [
            flat_bin_idx // bins[2] // bins[1] % bins[0],
            flat_bin_idx // bins[2] % bins[1],
            flat_bin_idx % bins[2]
        ]
        
        # Bin boundaries
        z_low = binning_config['phi_z'][bin_idx[0]]
        z_high = binning_config['phi_z'][bin_idx[0] + 1]
        pt_low = binning_config['phi_thrust_pt'][bin_idx[1]]
        pt_high = binning_config['phi_thrust_pt'][bin_idx[1] + 1]
        helicity_low = binning_config['phi_helicity_angle'][bin_idx[2]]
        helicity_high = binning_config['phi_helicity_angle'][bin_idx[2] + 1]
        
        # Check if bin directory exists
        bin_dir = f"./batch_comparison_images/bin_{flat_bin_idx:0{pad_width}d}"
        fit_plot = f"./fit_result/fit_plots/bin_{flat_bin_idx:0{pad_width}d}_fit.png"
        
        # Alternative fit plot locations
        if not os.path.exists(fit_plot):
            fit_plot = f"./fit_result/bin_{flat_bin_idx:0{pad_width}d}_fit.png"
        if not os.path.exists(fit_plot):
            fit_plot = f"./fit_result/plots/bin_{flat_bin_idx:0{pad_width}d}.png"
        
        # Check if at least some comparison plots exist
        has_plots = os.path.exists(bin_dir) and any(
            os.path.exists(os.path.join(bin_dir, f"{plot}.png")) 
            for plot in comparison_plots
        )
        
        if not has_plots:
            #print(f"Skipping bin {flat_bin_idx:0{pad_width}d}: No plots found.")
            unvalid_bins.append(flat_bin_idx)
            continue  # Skip bins without plots
        
        valid_bins += 1
        
        # Start new page for this bin
        latex_content.append(r"\clearpage")
        latex_content.append(r"\section*{Bin " + f"{flat_bin_idx:0{pad_width}d}" + r"}")
        latex_content.append(r"")
        
        # Bin information
        latex_content.append(r"\subsection*{Bin Definition}")
        latex_content.append(r"\begin{itemize}")
        latex_content.append(f"  \\item $z$: ${z_low:.3f} < z < {z_high:.3f}$")
        latex_content.append(f"  \\item $p_T$ (thrust): ${pt_low:.4f} < p_T < {pt_high:.4f}$ GeV/c")
        latex_content.append(f"  \\item Helicity angle: ${helicity_low:.2f} < \\cos\\theta_{{\\rm hel}} < {helicity_high:.2f}$")
        latex_content.append(f"  \\item Bin indices: $z_{{\\rm bin}}={bin_idx[0]}$, $p_T^{{\\rm bin}}={bin_idx[1]}$, $\\theta_{{\\rm bin}}={bin_idx[2]}$")
        latex_content.append(r"\end{itemize}")
        latex_content.append(r"")
        
        # Fit plot (if exists)
        if os.path.exists(fit_plot):
            latex_content.append(r"\subsection*{Mass Fit}")
            latex_content.append(r"\begin{figure}[h!]")
            latex_content.append(r"  \centering")
            latex_content.append(r"  \includegraphics[width=0.7\textwidth]{" + fit_plot + r"}")
            latex_content.append(r"  \caption{Mass fit for bin " + f"{flat_bin_idx:0{pad_width}d}" + r"}")
            latex_content.append(r"\end{figure}")
            latex_content.append(r"")
        
        # Comparison plots
        latex_content.append(r"\subsection*{Data/MC Comparisons}")
        latex_content.append(r"\begin{figure}[h!]")
        latex_content.append(r"  \centering")
        
        # Create 3x3 grid for 9 comparison plots
        for i, plot_name in enumerate(comparison_plots):
            plot_file = os.path.join(bin_dir, f"{plot_name}.png")
            
            if os.path.exists(plot_file):
                latex_content.append(r"  \begin{subfigure}[b]{0.32\textwidth}")
                latex_content.append(r"    \centering")
                latex_content.append(r"    \includegraphics[width=\textwidth]{" + plot_file + r"}")
                latex_content.append(r"    \caption{" + plot_labels.get(plot_name, plot_name) + r"}")
                latex_content.append(r"  \end{subfigure}")
                if (i + 1) % 3 == 0 and i < len(comparison_plots) - 1:
                    latex_content.append(r"  \\")  # New row after every 3 plots
            else:
                # Empty placeholder if plot doesn't exist
                latex_content.append(r"  \begin{subfigure}[b]{0.32\textwidth}")
                latex_content.append(r"    \centering")
                latex_content.append(r"    \vspace{3cm}")
                latex_content.append(r"    \caption{N/A}")
                latex_content.append(r"  \end{subfigure}")
                if (i + 1) % 3 == 0 and i < len(comparison_plots) - 1:
                    latex_content.append(r"  \\")
        
        latex_content.append(r"  \caption{Data/MC comparison plots for bin " + f"{flat_bin_idx:0{pad_width}d}" + r"}")
        latex_content.append(r"\end{figure}")
        latex_content.append(r"")
    
    # End document
    latex_content.append(r"\clearpage")
    latex_content.append(r"\section*{Summary}")
    latex_content.append(f"Total valid bins processed: {valid_bins} out of {total_bins}")
    latex_content.append(r"")
    latex_content.append(r"\end{document}")
    
    # Write to file
    output_file = "data_mc_comparison_summary.tex"
    with open(output_file, 'w') as f:
        f.write('\n'.join(latex_content))
    
    print(f"LaTeX document generated: {output_file}")
    print(f"Valid bins found: {valid_bins}/{total_bins}")
    print(f"\nTo compile:")
    print(f"  pdflatex {output_file}")
    print(f"  pdflatex {output_file}  # Run twice for TOC")
    print("\n unvalid bins: " , len(unvalid_bins))
    print(f"\nBins without valid plots: {unvalid_bins}")
    
    return output_file

if __name__ == "__main__":
    generate_latex_document()
