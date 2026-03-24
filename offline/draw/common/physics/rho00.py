"""rho00 extraction utilities: efficiency map and rho00 vs z.

Currently only support 2d z and helicity angle.
"""

import os
import math
import pandas as pd
import ROOT as R

from common.config.rho00 import (
    DEFAULT_BRANCH_MAP,
    DEFAULT_WEIGHT_BRANCH,
    DEFAULT_CSV_COL_MAP,
    RHO00_BINNING_CONFIG,
)


def _normalize_bin_cfg(cfg, name):
    """Allow binning to be provided as (nbins, min, max) or an explicit edge list."""
    if isinstance(cfg, tuple) and len(cfg) == 3:
        nbins, xmin, xmax = cfg
        return int(nbins), float(xmin), float(xmax)

    if isinstance(cfg, (list, tuple)):
        boundaries = list(cfg)
        if len(boundaries) < 2:
            raise ValueError(f"Invalid binning for {name}: need at least two boundaries, got {len(boundaries)}")
        return len(boundaries) - 1, float(boundaries[0]), float(boundaries[-1])

    raise ValueError(f"Invalid binning format for {name}: {cfg}")

def get_efficiency_hist(
    reco_rootfile,
    truth_rootfile,
    output_dir,
    branch_map=None,
    weight_branch=DEFAULT_WEIGHT_BRANCH,
    binning_config=None,
    cache_filename="efficiency_2D.root",
    particle_tag="particle",
):
    """Build or load efficiency TH2 (cosθ vs z) and save diagnostic plots.
    
    Only support 2D efficiency map for now, with helicity angle on x and z on y. 
    
    Args:
        branch_map: {"z": branch_name, "helicity_angle": branch_name}; falls back to defaults
        weight_branch: optional weight branch for reco; set None for unweighted
        binning_config: dict with keys "z" and "helicity_angle" mapping to (nbins, min, max)
        particle_tag: used in output figure names
    """
    os.makedirs(output_dir, exist_ok=True)
    eff_rootfile = os.path.join(output_dir, cache_filename)
    if branch_map is None:
        branch_map = DEFAULT_BRANCH_MAP
    if binning_config is None:
        binning_config = RHO00_BINNING_CONFIG

    z_cfg_raw = binning_config.get("z", RHO00_BINNING_CONFIG["z"])
    hel_cfg_raw = binning_config.get("helicity_angle", RHO00_BINNING_CONFIG["helicity_angle"])
    z_cfg = _normalize_bin_cfg(z_cfg_raw, "z")
    hel_cfg = _normalize_bin_cfg(hel_cfg_raw, "helicity_angle")

    z_branch = branch_map.get("z", DEFAULT_BRANCH_MAP["z"])
    helicity_branch = branch_map.get("helicity_angle", DEFAULT_BRANCH_MAP["helicity_angle"])

    if os.path.exists(eff_rootfile):
        eff_file = R.TFile.Open(eff_rootfile, "READ")
        hist_eff = eff_file.Get("hist_eff")
        hist_eff.SetDirectory(0)
        hist_reco = eff_file.Get("hist_reco")
    else:
        rdf_reco = R.RDataFrame("event", reco_rootfile)
        rdf_truth = R.RDataFrame("truth", truth_rootfile)

        hist_model = R.RDF.TH2DModel(
            "hist_eff",
            ";cos#theta;z",
            hel_cfg[0], hel_cfg[1], hel_cfg[2],
            z_cfg[0], z_cfg[1], z_cfg[2],
        )

        if weight_branch:
            hist_reco = rdf_reco.Histo2D(hist_model, helicity_branch, z_branch, weight_branch).GetValue()
        else:
            hist_reco = rdf_reco.Histo2D(hist_model, helicity_branch, z_branch).GetValue()
        hist_truth = rdf_truth.Histo2D(hist_model, helicity_branch, z_branch).GetValue()

        eff = R.TEfficiency(hist_reco, hist_truth)
        eff.SetName("eff_costheta_z")
        eff.SetTitle("Efficiency vs cos#theta and z;cos#theta;z;Efficiency")

        hist_eff = hist_reco.Clone("hist_eff")
        for ix in range(1, hist_eff.GetNbinsX() + 1):
            for iy in range(1, hist_eff.GetNbinsY() + 1):
                global_bin = eff.GetGlobalBin(ix, iy)
                eff_val = eff.GetEfficiency(global_bin)
                err_lo = eff.GetEfficiencyErrorLow(global_bin)
                err_hi = eff.GetEfficiencyErrorUp(global_bin)
                eff_err = (err_lo + err_hi) / 2.0
                hist_eff.SetBinContent(ix, iy, eff_val)
                hist_eff.SetBinError(ix, iy, eff_err)

        out_file = R.TFile(eff_rootfile, "RECREATE")
        hist_eff.Write()
        hist_reco.Write("hist_reco")
        out_file.Close()

    R.gStyle.SetOptStat(0)
    canvas = R.TCanvas("c_eff_final", "c_eff_final", 800, 600)
    R.gStyle.SetPaintTextFormat(".3f")
    hist_eff.Draw("COLZ TEXT")
    canvas.SaveAs(os.path.join(output_dir, f"{particle_tag}_efficiency_costheta_vs_z.png"))

    canvas.cd()
    R.gStyle.SetPaintTextFormat(".0f")
    hist_reco.Draw("COLZ TEXT")
    canvas.SaveAs(os.path.join(output_dir, f"{particle_tag}_reco_count_costheta_vs_z.png"))

    return hist_eff


def extract_rho00(
    nsig_csv,
    output_dir,
    eff_hist,
    fit_func,
    csv_col_map=None,
    binning_config=None,
    bin_half_width=0.025,
    particle_tag="particle",
):
    """Fit rho00 vs z using nsig CSV and efficiency histogram.

    Args:
        csv_col_map: {"z": col, "helicity_angle": col, "nsig": col, "nsig_err": col}; defaults provided
        binning_config: dict with keys "z" and "helicity_angle" mapping to (nbins, min, max)
        particle_tag: used in output figure names and logs
    """
    from DRAW import style_draw, HistStyle  # local import to avoid hard dependency at import time

    os.makedirs(output_dir, exist_ok=True)
    if csv_col_map is None:
        csv_col_map = DEFAULT_CSV_COL_MAP
    if binning_config is None:
        binning_config = RHO00_BINNING_CONFIG

    z_cfg_raw = binning_config.get("z", RHO00_BINNING_CONFIG["z"])
    hel_cfg_raw = binning_config.get("helicity_angle", RHO00_BINNING_CONFIG["helicity_angle"])
    z_cfg = _normalize_bin_cfg(z_cfg_raw, "z")
    hel_cfg = _normalize_bin_cfg(hel_cfg_raw, "helicity_angle")

    z_center_col = csv_col_map.get("z", DEFAULT_CSV_COL_MAP["z"])
    helicity_center_col = csv_col_map.get("helicity_angle", DEFAULT_CSV_COL_MAP["helicity_angle"])
    nsig_col = csv_col_map.get("nsig", DEFAULT_CSV_COL_MAP["nsig"])
    nsig_err_col = csv_col_map.get("nsig_err", DEFAULT_CSV_COL_MAP["nsig_err"])

    df = pd.read_csv(nsig_csv)

    df_valid = df[(df[nsig_col] > 0) | (df[z_center_col] > 0)].reset_index(drop=True)
    df_valid = df_valid[(df_valid[z_center_col] >= z_cfg[1]) & (df_valid[z_center_col] <= z_cfg[2])].reset_index(drop=True)

    grouped = df_valid.groupby(z_center_col)
    unique_z = sorted(x for x in df_valid[z_center_col].unique() if x > 0)

    z_values, rho00_values, rho00_errors = [], [], []

    for idx, z_val in enumerate(unique_z):
        group = grouped.get_group(z_val)
        hist = R.TH1F(
            f"hist_z_{z_val:.3f}",
            f"z = {z_val:.3f};cos#theta^{{*}};N_{'{'}sig{'}'}/#varepsilon",
            hel_cfg[0], hel_cfg[1], hel_cfg[2],
        )

        for _, row in group.iterrows():
            cos_theta = row[helicity_center_col]
            z_center = row[z_center_col]
            nsig = row[nsig_col]
            nsig_err = row[nsig_err_col]

            bin_x = eff_hist.GetXaxis().FindBin(cos_theta)
            bin_y = eff_hist.GetYaxis().FindBin(z_center)
            efficiency = eff_hist.GetBinContent(bin_x, bin_y)
            eff_error = eff_hist.GetBinError(bin_x, bin_y)
            if efficiency > 1:
                print(f"Warning: Efficiency > 1 at cos_theta={cos_theta}, z={z_center} : {efficiency}")
            if math.isnan(eff_error):
                eff_error = 0.001

            if efficiency > 0:
                corrected_nsig = nsig / efficiency
                corrected_err = corrected_nsig * R.TMath.Sqrt((nsig_err / nsig) ** 2 + (eff_error / efficiency) ** 2) if nsig > 0 else 0
            else:
                corrected_nsig = 0
                corrected_err = 0

            bin_idx = hist.FindBin(cos_theta)
            hist.SetBinContent(bin_idx, corrected_nsig)
            hist.SetBinError(bin_idx, corrected_err)

        value, err = fit_func(
            hist,
            os.path.join(output_dir, f"fit_{idx:02d}.png"),
            True,
            extra_legend=f"{z_val - bin_half_width:.2f} < z < {z_val + bin_half_width:.2f}",
        )
        z_values.append(z_val)
        rho00_values.append(value)
        rho00_errors.append(err)

    hist_rho00_z = R.TH1F("hist_rho00_z", ";z;#rho_{00}", 20, 0, 1)
    for z_val, rho_val, rho_err in zip(z_values, rho00_values, rho00_errors):
        bin_x = hist_rho00_z.GetXaxis().FindBin(z_val)
        hist_rho00_z.SetBinContent(bin_x, rho_val)
        hist_rho00_z.SetBinError(bin_x, rho_err)

    canvas = style_draw([hist_rho00_z], "", styles=[HistStyle.error_bars(R.kBlack)], y_min=0, y_max=0.6, use_user_y_range=True, save=False)
    line = R.TLine(0, 1 / 3, 1, 1 / 3)
    line.SetLineColor(R.kRed)
    line.SetLineStyle(R.kDashed)
    line.SetLineWidth(2)
    line.Draw("same")

    leg = R.TLegend(0.6, 0.3, 0.9, 0.55)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(line, "#rho_{00} = 1/3 (unpolarized)", "l")
    leg.Draw("same")
    canvas.SaveAs(os.path.join(output_dir, f"rho00_vs_z_{particle_tag}.png"))

    with open(os.path.join(output_dir, f"rho00_results_{particle_tag}.txt"), "w") as fout:
        fout.write("# z_center  rho00  rho00_error\n")
        for z_val, rho_val, rho_err in zip(z_values, rho00_values, rho00_errors):
            fout.write(f"{z_val:.4f}  {rho_val:.6f}  {rho_err:.6f}\n")
    print(f"fitting result has been saved to rho00_results_{particle_tag}.txt")

    hist_rho00_z.SetDirectory(0)
    return hist_rho00_z
