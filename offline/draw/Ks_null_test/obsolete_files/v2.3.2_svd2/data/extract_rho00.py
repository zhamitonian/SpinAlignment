import pandas as pd 
import os 
import sys
import ROOT as R
from DRAW import style_draw, HistStyle
from FIT import fit_rho00
import math

def get_2d_hist_from_txt(nsig_txt):
    df_nsig = pd.read_csv(nsig_txt, sep=r"\s+",
                        names=["z_center", "z_width", "helicity_center", 
                               "helicity_width", "nsig", "nsig_err", "nsig_err2"],
                        skiprows=1)
    
    df_nsig = df_nsig[(df_nsig['nsig'] > 0)]

    hist_nsig = R.TH2D("hist_nsig", "nsig in cos#theta vs Ks z;cos#theta*;Ks z", 10, -1, 1, 14, 0.05, 0.75)

    for index, row in df_nsig.iterrows():
        binx = hist_nsig.GetXaxis().FindBin(row['helicity_center'])
        biny = hist_nsig.GetYaxis().FindBin(row['z_center'])

        hist_nsig.SetBinContent(binx, biny, row['nsig'])
        hist_nsig.SetBinError(binx, biny, row['nsig_err'])

    R.gStyle.SetOptStat(0)
    c2 = R.TCanvas("c_eff_nsig", "c_eff_nsig", 1600, 1080)
    hist_nsig.Draw("COLZ TEXT")
    c2.SaveAs(os.path.join("./", "Ks_nsig_costheta_vs_z.png"))

    return hist_nsig

def calculate_eff_original(truth_rootFile, nsig_txt, output_dir):
    rdf_truth = R.RDataFrame("truth", truth_rootFile )
    hist_model = R.RDF.TH2DModel("", ";cos#theta*;z", 10, -1, 1, 15, 0.05, 0.8)

    #hist_cosTheta_z = rdf_truth.Histo2D(hist_model , "Ks_helicity_angle", "Ks_z")
    hist_cosTheta_z = rdf_truth.Histo2D(hist_model , "Ks_helicity_angle", "Ks_z", "Ks_M")

    R.gStyle.SetOptStat(0)

    c = R.TCanvas("c_eff", "c_eff", 800, 600)
    hist_cosTheta_z.Draw("COLZ TEXT")
    c.SaveAs(os.path.join(output_dir, "Ks_costheta_vs_z.png"))
    
    # read nsig extracted by fitting frin ./MC_fitting/nsig_results.txt
    df_nsig = pd.read_csv(nsig_txt, sep=r"\s+",
                        names=["z_center", "z_width", "helicity_center", 
                               "helicity_width", "nsig", "nsig_err", "nsig_err2"],
                        skiprows=1)
    
    df_nsig = df_nsig[(df_nsig['nsig'] > 0)]

    hist_nsig = R.TH2D("hist_nsig", "nsig in cos#theta vs Ks z;cos#theta*;Ks z", 10, -1, 1, 15, 0.05, 0.8)

    for index, row in df_nsig.iterrows():
        binx = hist_nsig.GetXaxis().FindBin(row['helicity_center'])
        biny = hist_nsig.GetYaxis().FindBin(row['z_center'])

        hist_nsig.SetBinContent(binx, biny, row['nsig'])
        hist_nsig.SetBinError(binx, biny, row['nsig_err'])
    
    c2 = R.TCanvas("c_eff_nsig", "c_eff_nsig", 1600, 1080)
    hist_nsig.Draw("COLZ TEXT")
    c2.SaveAs(os.path.join(output_dir, "Ks_nsig_costheta_vs_z.png"))

    # Use TEfficiency for correct error calculation
    # TEfficiency takes (passed, total) histograms
    eff = R.TEfficiency(hist_nsig, hist_cosTheta_z.GetPtr())
    eff.SetName("eff_costheta_z")
    eff.SetTitle("Efficiency vs cos#theta and z;cos#theta;Ks z;Efficiency")
    
    # For visualization, create a histogram from TEfficiency
    hist_eff = hist_nsig.Clone("hist_eff")
    hist_eff.SetTitle("Efficiency vs cos#theta and z;cos#theta*; z")
    
    # Fill efficiency histogram with correct values and errors
    for i in range(1, hist_eff.GetNbinsX() + 1):
        for j in range(1, hist_eff.GetNbinsY() + 1):
            global_bin = eff.GetGlobalBin(i, j)
            eff_value = eff.GetEfficiency(global_bin)
            eff_error_low = eff.GetEfficiencyErrorLow(global_bin)
            eff_error_high = eff.GetEfficiencyErrorUp(global_bin)
            
            # Use average of asymmetric errors for symmetric error bar
            eff_error = (eff_error_low + eff_error_high) / 2.0
            
            hist_eff.SetBinContent(i, j, eff_value)
            hist_eff.SetBinError(i, j, eff_error)

    c3 = R.TCanvas("c_eff_final", "c_eff_final", 800, 600)
    #R.gStyle.SetPaintTextFormat(".2%")
    R.gStyle.SetPaintTextFormat(".3f")
    hist_eff.Draw("COLZ TEXT")
    c3.SaveAs(os.path.join(output_dir, "Ks_efficiency_costheta_vs_z.png"))
    
    # Also save the TEfficiency object to a ROOT file for later use
    output_file = R.TFile(os.path.join(output_dir ,"efficiency_2D.root"), "RECREATE")
    eff.Write()
    hist_eff.Write()
    hist_cosTheta_z.GetPtr().Write()
    hist_nsig.Write()
    output_file.Close()
    
    return hist_eff, eff

def extract_rho00(nsig_txt:str, output_dir:str, eff_root:str = None):
    df = pd.read_csv(nsig_txt, sep=r"\s+",
                    names=["z_center", "z_width", "helicity_center", 
                            "helicity_width", "nsig", "nsig_err", "nsig_err2"],
                    skiprows=1,
                    comment="#")

    # filter invalid data
    df_valid = df[(df['nsig'] > 0) | (df['z_center'] > 0)].reset_index(drop=True)
    
    # the range i wanna use
    df_valid = df_valid[(df_valid['z_center'] >= 0.05) & (df_valid['z_center'] <= 0.8)].reset_index(drop=True)

    grouped = df_valid.groupby('z_center')

    unique_z = sorted(df_valid['z_center'].unique())
    unique_z = [x for x in unique_z if x > 0]

    # get efficiency from root file
    if eff_root is not None:
        eff_file = R.TFile.Open(eff_root, "READ")
        eff_hist = eff_file.Get("hist_eff")
    else :
        eff_hist = None

    # 存储拟合结果
    z_values = []
    rho00_values = []
    rho00_errors = []

    # 对每个 z_center 值绘图和拟合
    for idx, z_value in enumerate(unique_z):

        group = grouped.get_group(z_value)
        #max = group['nsig'].max()
        #min = group['nsig'].min()

        hist = R.TH1F(f"hist_z_{z_value:.3f}", f"z = {z_value:.3f};" + "cos#theta^{*};N_{sig}/#varepsilon", 10, -1, 1)
        
        for i in range(len(group)):
            cos_theta = group['helicity_center'].values[i]
            z_center = group['z_center'].values[i]
            nsig = group['nsig'].values[i]
            nsig_err = group['nsig_err'].values[i]
            
            # Apply efficiency correction if available
            if eff_hist is not None:
                # For 2D efficiency histogram (cos_theta vs z)
                bin_x = eff_hist.GetXaxis().FindBin(cos_theta)
                bin_y = eff_hist.GetYaxis().FindBin(z_center)
                efficiency = eff_hist.GetBinContent(bin_x, bin_y)
                eff_error = eff_hist.GetBinError(bin_x, bin_y)
                if efficiency > 1:
                    print(f"Warning: Efficiency > 1 at cos_theta={cos_theta}, z={z_center} : {efficiency}")
                if math.isnan(eff_error):
                    eff_error = 0.001
                
                # Avoid division by zero
                if efficiency > 0:
                    # Corrected nsig = nsig / efficiency
                    corrected_nsig = nsig / efficiency
                    # Error propagation: delta(N/eff) = sqrt((deltaN/eff)^2 + (N*deltaEff/eff^2)^2)
                    corrected_err = corrected_nsig * R.TMath.Sqrt(
                        (nsig_err / nsig)**2 + (eff_error / efficiency)**2
                    ) if nsig > 0 else 0
                else:
                    # If efficiency is 0, skip this bin or use uncorrected values
                    corrected_nsig = 0# nsig
                    corrected_err = 0# nsig_err
            else:
                # No efficiency correction
                corrected_nsig = nsig
                corrected_err = nsig_err
            
            bin_idx = hist.FindBin(cos_theta)
            hist.SetBinContent(bin_idx, corrected_nsig)
            hist.SetBinError(bin_idx, corrected_err)

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
    #leg.AddEntry(hist_rho00_z, "#rho_{00} from fits", "ep")  # "ep" 表示误差棒和点

    leg.Draw("same")
    c.SaveAs(os.path.join(output_dir, "rho00_vs_z.png"))
        
    with open(os.path.join(output_dir, "rho00_results.txt"), 'w') as f:
        f.write('# z_center  rho00  rho00_error\n')
        for z, rho00, rho00_err in zip(z_values, rho00_values, rho00_errors):
            f.write(f'{z:.4f}  {rho00:.6f}  {rho00_err:.6f}\n')
    print(f"fitting result has been saved to rho00_results.txt")

    hist_rho00_z.SetDirectory(0)  # Detach from any file

    return hist_rho00_z
        

if __name__ == "__main__":
    truth_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar_svd2/svd2_st0_truth_processed.root"
    nsig_txt = "./output_cbshape/nsig_results.txt"

    nsig_mc_txt = "../qqbar/output_st0/nsig_results.txt"

    eff_dir = "../qqbar/eff"
    os.makedirs(eff_dir, exist_ok=True)
    
    if not os.path.exists(os.path.join(eff_dir, "efficiency_2D.root")):
        calculate_eff_original(truth_rootFile, nsig_mc_txt, eff_dir)

    rho00_output_dir = "./rho00"
    os.makedirs(rho00_output_dir, exist_ok=True)
    
    h1 = extract_rho00(nsig_txt, rho00_output_dir, os.path.join(eff_dir, "efficiency_2D.root"))