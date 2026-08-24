import pandas as pd 
import os 
import sys
import ROOT as R
from DRAW import style_draw, HistStyle
from FIT import fit_rho00
import math

#Z_BIN_CONFIG = (16, 0.20, 0.90)
Z_BIN_CONFIG = (15, 0.25, 0.90)
COSTHETA_BIN_CONFIG = (10, -1, 1)
#COSTHETA_BIN_CONFIG = (8, -0.8, 0.8)

HELICITY_BRANCH_NAME = "phi_helicity_angle"
Z_BRANCH_NAME = "phi_z"

def get_efficiency_hist(mix_rootFile, output_dir):
    eff_rootFile = os.path.join(output_dir, "detector_effect_2D.root")

    if not os.path.exists(eff_rootFile):

        rdf_reco = R.RDataFrame("event", mix_rootFile)

        hist_model = R.RDF.TH2DModel("hist_mix", ";cos#theta;z", COSTHETA_BIN_CONFIG[0], COSTHETA_BIN_CONFIG[1], COSTHETA_BIN_CONFIG[2],
                                            Z_BIN_CONFIG[0], Z_BIN_CONFIG[1], Z_BIN_CONFIG[2])

        hist_reco = rdf_reco.Histo2D(hist_model, HELICITY_BRANCH_NAME, Z_BRANCH_NAME).GetValue()

        output_file = R.TFile(eff_rootFile, "RECREATE")

        hist_reco.Write("hist_mix")    
        output_file.Close()
    else : 
        eff_file = R.TFile.Open(eff_rootFile, "READ")
        hist_mix = eff_file.Get("hist_mix")
        hist_mix.SetDirectory(0)  # Detach from file

    R.gStyle.SetOptStat(0)
    c = R.TCanvas("c_eff_final", "c_eff_final", 800, 600)
    R.gStyle.SetPaintTextFormat(".3f")
    hist_mix.Draw("COLZ TEXT")
    c.SaveAs(os.path.join(output_dir, "efficiency_costheta_vs_z.png"))

    return hist_mix
    

def extract_rho00(nsig_txt:str, output_dir:str, eff_hist:R.TH2):
    df = pd.read_csv(nsig_txt)

    z_center_col = Z_BRANCH_NAME + "_center"
    helicity_col = HELICITY_BRANCH_NAME + "_center"

    # filter invalid data
    df_valid = df[(df['nsig'] > 0) | (df[z_center_col] > 0)].reset_index(drop=True)
    
    # the range i wanna use
    df_valid = df_valid[(df_valid[z_center_col] >= Z_BIN_CONFIG[1]) & (df_valid[z_center_col] <= Z_BIN_CONFIG[2])].reset_index(drop=True)

    grouped = df_valid.groupby(z_center_col)

    unique_z = sorted(df_valid[z_center_col].unique())
    unique_z = [x for x in unique_z if x > 0]

    z_values = []
    rho00_values = []
    rho00_errors = []

    for idx, z_value in enumerate(unique_z):

        group = grouped.get_group(z_value)

        hist = R.TH1F(f"hist_z_{z_value:.3f}", f"z = {z_value:.3f};" + "cos#theta^{*};N_{sig}/#varepsilon", COSTHETA_BIN_CONFIG[0], COSTHETA_BIN_CONFIG[1], COSTHETA_BIN_CONFIG[2])
        
        for i in range(len(group)):
            cos_theta = group[helicity_col].values[i]
            z_center = group[z_center_col].values[i]
            nsig = group['nsig'].values[i]
            nsig_err = group['nsig_err_hi'].values[i]
            
            bin_x = eff_hist.GetXaxis().FindBin(cos_theta)
            bin_y = eff_hist.GetYaxis().FindBin(z_center)
            efficiency = eff_hist.GetBinContent(bin_x, bin_y)
            eff_error = eff_hist.GetBinError(bin_x, bin_y)
            if efficiency > 0:
                corrected_nsig = nsig / efficiency
                # Error propagation: delta(N/eff) = sqrt((deltaN/eff)^2 + (N*deltaEff/eff^2)^2)
                corrected_err = corrected_nsig * R.TMath.Sqrt(
                    (nsig_err / nsig)**2 + (eff_error / efficiency)**2
                ) if nsig > 0 else 0
            else:
                corrected_nsig = 0# nsig
                corrected_err = 0# nsig_err
        
            bin_idx = hist.FindBin(cos_theta)
            hist.SetBinContent(bin_idx, corrected_nsig)
            hist.SetBinError(bin_idx, corrected_err)
            #hist.SetBinContent(bin_idx, nsig)
            #hist.SetBinError(bin_idx, nsig_err)
        style_draw([hist],"temp.png")

        value, err = fit_rho00(hist, os.path.join(output_dir, f"fit_{idx:02d}.png"), True, extra_legend=f"{z_value - 0.025:.2f} < z < {z_value + 0.025:.2f}")
        z_values.append(z_value)
        rho00_values.append(value)
        rho00_errors.append(err)
    
    hist_rho00_z = R.TH1F("hist_rho00_z", ";z;#rho_{00}", 20, 0, 1)
    for i in range(len(z_values)):
        bin_x = hist_rho00_z.GetXaxis().FindBin(z_values[i])
        hist_rho00_z.SetBinContent(bin_x, rho00_values[i])
        hist_rho00_z.SetBinError(bin_x, rho00_errors[i])

    c = style_draw([hist_rho00_z], "", styles=[HistStyle.error_bars(R.kBlack)], y_min = 0, y_max=0.6, use_user_y_range=True, save = False)
    line = R.TLine(0, 1/3, 1, 1/3)
    line.SetLineColor(R.kRed)
    line.SetLineStyle(R.kDashed)
    line.SetLineWidth(2)
    line.Draw("same")

    leg = R.TLegend(0.6, 0.3, 0.9, 0.55)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(line, "#rho_{00} = 1/3(unpolarized)", "l")

    leg.Draw("same")
    c.SaveAs(os.path.join(output_dir, "rho00_vs_z.png"))
        
    with open(os.path.join(output_dir, "rho00_results.csv"), 'w') as f:
        f.write('z_center,rho00,rho00_error\n')
        for z, rho00, rho00_err in zip(z_values, rho00_values, rho00_errors):
            f.write(f'{z:.4f},{rho00:.6f},{rho00_err:.6f}\n')
    print(f"fitting result has been saved to rho00_results.csv")

    hist_rho00_z.SetDirectory(0)  # Detach from any file

    return hist_rho00_z
        

if __name__ == "__main__":
    nsig_txt = "../Extract_rho00/fit_results/data_fit/nsig_results.csv"
    output_dir = "./images/rho00/data_bin10/"
    os.makedirs(output_dir, exist_ok=True)
    mix_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/PhiSpinAlignment_v2.1.3_data_svd2/test/reco_data_mixed.root"
    hist_mix = get_efficiency_hist(mix_rootFile, output_dir)
    extract_rho00(nsig_txt, output_dir, hist_mix)