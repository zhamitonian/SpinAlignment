import pandas as pd 
import os 
import sys
import ROOT as R
from DRAW import style_draw, HistStyle
from FIT import fit_rerho1m1
import math

Z_BIN_CONFIG = (14, 0.2, 0.9)
PHI_BIN_CONFIG = (20, 0, math.pi/2)

Z_BIN_NAME = "phi_z"
PHI_BIN_NAME = "phi_helicity_phi"
PID_WEIGHT_NAME = "phi_weight"


def get_efficiency_hist_from_truth_match(reco_rootFile, truth_rootFile, output_dir):
    eff_rootFile = os.path.join(output_dir, "efficiency_2D.root")

    if not os.path.exists(eff_rootFile):

        rdf_reco = R.RDataFrame("event", reco_rootFile)
        rdf_truth = R.RDataFrame("truth", truth_rootFile)

        hist_model = R.RDF.TH2DModel("hist_eff", ";#phi;z", PHI_BIN_CONFIG[0], PHI_BIN_CONFIG[1], PHI_BIN_CONFIG[2],
                                            Z_BIN_CONFIG[0], Z_BIN_CONFIG[1], Z_BIN_CONFIG[2])

        hist_reco = rdf_reco.Histo2D(hist_model, PHI_BIN_NAME, Z_BIN_NAME, PID_WEIGHT_NAME).GetValue()
        hist_truth = rdf_truth.Histo2D(hist_model, PHI_BIN_NAME, Z_BIN_NAME).GetValue()

        eff = R.TEfficiency(hist_reco, hist_truth)
        eff.SetName("eff_phi_z")
        eff.SetTitle("Efficiency vs #phi and z;#phi;z;Efficiency")
        
        # For visualization, create a histogram from TEfficiency
        hist_eff = hist_reco.Clone("hist_eff")
        for i in range(1, hist_eff.GetNbinsX() + 1):
            for j in range(1, hist_eff.GetNbinsY() + 1):
                global_bin = eff.GetGlobalBin(i, j)
                eff_value = eff.GetEfficiency(global_bin)
                eff_error_low = eff.GetEfficiencyErrorLow(global_bin)
                eff_error_high = eff.GetEfficiencyErrorUp(global_bin)
                
                eff_error = (eff_error_low + eff_error_high) / 2.0
                
                hist_eff.SetBinContent(i, j, eff_value)
                hist_eff.SetBinError(i, j, eff_error)

        output_file = R.TFile(eff_rootFile, "RECREATE")

        hist_eff.Write()
        hist_reco.Write("hist_reco")    
        output_file.Close()
    else : 
        eff_file = R.TFile.Open(eff_rootFile, "READ")
        hist_eff = eff_file.Get("hist_eff")
        hist_eff.SetDirectory(0)  # Detach from file

        hist_reco = eff_file.Get("hist_reco")

    R.gStyle.SetOptStat(0)
    c = R.TCanvas("c_eff_final", "c_eff_final", 800, 600)
    R.gStyle.SetPaintTextFormat(".3f")
    hist_eff.Draw("COLZ TEXT")
    c.SaveAs(os.path.join(output_dir, "efficiency_phi_vs_z.png"))

    return hist_eff


def get_efficiency_hist_from_fit(csv_file, truth_rootFile, output_dir):

    eff_rootFile = os.path.join(output_dir, "efficiency_2D.root")

    z_bin_name = Z_BIN_NAME + "_center"
    phi_bin_name = PHI_BIN_NAME + "_center"

    df = pd.read_csv(csv_file)
    df_valid = df[(df['nsig'] > 0) | (df[z_bin_name] > 0)].reset_index(drop=True)
    df_valid = df_valid[(df_valid[z_bin_name] >= Z_BIN_CONFIG[1]) & (df_valid[z_bin_name] <= Z_BIN_CONFIG[2])].reset_index(drop=True)

    if not os.path.exists(eff_rootFile):

        # Build hist_reco from CSV nsig (fit-based reco yield)
        hist_reco = R.TH2D("hist_reco", ";#phi;z",
                           PHI_BIN_CONFIG[0], PHI_BIN_CONFIG[1], PHI_BIN_CONFIG[2],
                           Z_BIN_CONFIG[0],   Z_BIN_CONFIG[1],   Z_BIN_CONFIG[2])
        for _, row in df_valid.iterrows():
            phi_c = row[phi_bin_name]
            z_c   = row[z_bin_name]
            nsig  = row['nsig']
            nsig_err = row['nsig_err_hi']
            bx = hist_reco.GetXaxis().FindBin(phi_c)
            by = hist_reco.GetYaxis().FindBin(z_c)
            hist_reco.SetBinContent(bx, by, max(nsig, 0))
            hist_reco.SetBinError(bx, by, nsig_err)

        # Build hist_truth from truth ROOT file
        rdf_truth = R.RDataFrame("truth", truth_rootFile)
        hist_truth = rdf_truth.Histo2D(
            R.RDF.TH2DModel("hist_truth", ";#phi;z",
                            PHI_BIN_CONFIG[0], PHI_BIN_CONFIG[1], PHI_BIN_CONFIG[2],
                            Z_BIN_CONFIG[0],   Z_BIN_CONFIG[1],   Z_BIN_CONFIG[2]),
            PHI_BIN_NAME, Z_BIN_NAME).GetValue()

        # Compute efficiency = nsig_fit / n_truth
        hist_eff = hist_reco.Clone("hist_eff")
        hist_eff.Reset()
        for i in range(1, hist_eff.GetNbinsX() + 1):
            for j in range(1, hist_eff.GetNbinsY() + 1):
                n_truth = hist_truth.GetBinContent(i, j)
                n_reco  = hist_reco.GetBinContent(i, j)
                n_reco_err = hist_reco.GetBinError(i, j)
                if n_truth > 0:
                    eff_val = n_reco / n_truth
                    eff_err = n_reco_err / n_truth
                else:
                    eff_val, eff_err = 0.0, 0.0
                hist_eff.SetBinContent(i, j, eff_val)
                hist_eff.SetBinError(i, j, eff_err)

        output_file = R.TFile(eff_rootFile, "RECREATE")
        hist_eff.Write()
        hist_reco.Write("hist_reco")
        output_file.Close()
    else:
        eff_file = R.TFile.Open(eff_rootFile, "READ")
        hist_eff = eff_file.Get("hist_eff")
        hist_eff.SetDirectory(0)
        hist_reco = eff_file.Get("hist_reco")

    R.gStyle.SetOptStat(0)
    c = R.TCanvas("c_eff_final", "c_eff_final", 800, 600)
    R.gStyle.SetPaintTextFormat(".3f")
    hist_eff.Draw("COLZ TEXT")
    c.SaveAs(os.path.join(output_dir, "efficiency_phi_vs_z.png"))

    return hist_eff


def extract_rerho1m1(csv_file:str, output_dir:str, eff_hist:R.TH2):
    df = pd.read_csv(csv_file)

    z_bin_name = Z_BIN_NAME + "_center"
    phi_bin_name = PHI_BIN_NAME + "_center"

    # filter invalid data
    df_valid = df[(df['nsig'] > 0) | (df[z_bin_name] > 0)].reset_index(drop=True)
    
    # the range i wanna use
    df_valid = df_valid[(df_valid[z_bin_name] >= Z_BIN_CONFIG[1]) & (df_valid[z_bin_name] <= Z_BIN_CONFIG[2])].reset_index(drop=True)

    grouped = df_valid.groupby(z_bin_name)

    unique_z = sorted(df_valid[z_bin_name].unique())
    unique_z = [x for x in unique_z if x > 0]

    z_values = []
    rho00_values = []
    rho00_errors = []

    for idx, z_value in enumerate(unique_z):

        group = grouped.get_group(z_value)

        hist = R.TH1F(f"hist_z_{z_value:.3f}", f"z = {z_value:.3f};" + "#phi^{*};N_{sig}/#varepsilon", PHI_BIN_CONFIG[0], PHI_BIN_CONFIG[1], PHI_BIN_CONFIG[2])
        
        for i in range(len(group)):
            phi = group[phi_bin_name].values[i]
            z_center = group[z_bin_name].values[i]
            nsig = group['nsig'].values[i]
            nsig_err = group['nsig_err_hi'].values[i]
            
            bin_x = eff_hist.GetXaxis().FindBin(phi)
            bin_y = eff_hist.GetYaxis().FindBin(z_center)
            efficiency = eff_hist.GetBinContent(bin_x, bin_y)
            eff_error = eff_hist.GetBinError(bin_x, bin_y)
            if efficiency > 1:
                print(f"Warning: Efficiency > 1 at phi={phi}, z={z_center} : {efficiency}")
            if math.isnan(eff_error):
                eff_error = 0.001
            
            if efficiency > 0:
                corrected_nsig = nsig / efficiency
                # Error propagation: delta(N/eff) = sqrt((deltaN/eff)^2 + (N*deltaEff/eff^2)^2)
                corrected_err = corrected_nsig * R.TMath.Sqrt(
                    (nsig_err / nsig)**2 + (eff_error / efficiency)**2
                ) if nsig > 0 else 0
            else:
                corrected_nsig = 0# nsig
                corrected_err = 0# nsig_err
        
            bin_idx = hist.FindBin(phi)
            hist.SetBinContent(bin_idx, corrected_nsig)
            hist.SetBinError(bin_idx, corrected_err)

        value, err = fit_rerho1m1(hist, os.path.join(output_dir, f"fit_{idx:02d}.png"), True, extra_legend=f"{z_value - 0.025:.2f} < z < {z_value + 0.025:.2f}")
        z_values.append(z_value)
        rho00_values.append(value)
        rho00_errors.append(err)
    
    hist_rerho1m1_z = R.TH1F("hist_rerho1m1_z", ";z;Re#rho_{1,-1}", 20, 0, 1)
    for i in range(len(z_values)):
        bin_x = hist_rerho1m1_z.GetXaxis().FindBin(z_values[i])
        hist_rerho1m1_z.SetBinContent(bin_x, rho00_values[i])
        hist_rerho1m1_z.SetBinError(bin_x, rho00_errors[i])

    c = style_draw([hist_rerho1m1_z], "", styles=[HistStyle.error_bars(R.kBlack)], y_min = -0.5, y_max=0.5, use_user_y_range=True, save = False)
    line = R.TLine(0, 0, 1, 0)
    line.SetLineColor(R.kRed)
    line.SetLineStyle(R.kDashed)
    line.SetLineWidth(2)
    line.Draw("same")

    leg = R.TLegend(0.6, 0.3, 0.9, 0.55)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(line, "Re#rho_{1,-1} = 0 (unpolarized)", "l")

    leg.Draw("same")
    c.SaveAs(os.path.join(output_dir, "rerho1m1_vs_z.png"))
        
    with open(os.path.join(output_dir, "rerho1m1_results.csv"), 'w') as f:
        f.write('#z_center,rerho1m1,rerho1m1_error\n')
        for z, rerho1m1, rerho1m1_err in zip(z_values, rho00_values, rho00_errors):
            f.write(f'{z:.4f},{rerho1m1:.6f},{rerho1m1_err:.6f}\n')
    print(f"fitting result has been saved to rerho1m1_results.csv")

    hist_rerho1m1_z.SetDirectory(0)  # Detach from any file

    return hist_rerho1m1_z
        

if __name__ == "__main__":
    truth_rootFile = "/gpfs/home/belle2/wangz/truth.root"
    reco_rootFile = "/gpfs/home/belle2/wangz/reco_truth_match.root"
    nsig_txt = "./fit_results/data_fit/nsig_results.csv"
    MC_nsig_txt = "./fit_results/MC_fit/nsig_results.csv"
    output_dir = "./images/rerho1m1/data_bin10/"
    os.makedirs(output_dir, exist_ok=True)

    #hist_eff  = get_efficiency_hist_from_truth_match(reco_rootFile, truth_rootFile, "./images/") 
    #extract_rerho1m1(nsig_txt, output_dir, hist_eff)

    hist_eff = get_efficiency_hist_from_fit(MC_nsig_txt, truth_rootFile, "./images/")
    extract_rerho1m1(nsig_txt, output_dir, hist_eff)