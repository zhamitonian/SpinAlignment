import pandas as pd 
import os 
import sys
import ROOT as R
from DRAW import style_draw, HistStyle
from FIT import fit_rho00

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

def truth_level_rho00(truth_rootFile,output_dir):
    os.makedirs(output_dir, exist_ok=True)

    rdf = R.RDataFrame("truth", truth_rootFile )
    z_values = []
    rho00_values = []
    rho00_errors = []

    hist_model = R.RDF.TH1DModel("", ";cos#theta^{*};N_{sig}", 10, -1, 1)
    for i in range(15):
        #hist = rdf.Filter(f"Ks_z >= {0.05 + i*0.05} && Ks_z < {0.05 + (i+1)*0.05}").Histo1D(hist_model , "Ks_helicity_angle").GetValue()
        hist = rdf.Define("temp",f"Ks_helicity_angle[Ks_z >= {0.05 + i*0.05} && Ks_z < {0.05 + (i+1)*0.05}]").Histo1D(hist_model , "temp").GetValue()
        value, err = fit_rho00(hist, os.path.join(output_dir, f"fit_truth_{i:02d}.png"), True, extra_legend=f"{0.05 + i*0.05:.2f} < z < {0.05 + (i+1)*0.05:.2f}")
        z_values.append(0.05 + i*0.05 + 0.025)
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
    

def calculate_eff_original(truth_rootFile, nsig_txt, output_dir):
    rdf_truth = R.RDataFrame("truth", truth_rootFile )
    hist_model = R.RDF.TH2DModel("", ";cos#theta*;z", 10, -1, 1, 15, 0.05, 0.8)
    hist_cosTheta_z = rdf_truth.Histo2D(hist_model , "Ks_helicity_angle", "Ks_z")

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


def extract_eff(reco_truth_matched_rootFile, truth_rootFile, eff_root="efficiency_2D.root", output_dir = "."):
    """
    calculate 2D efficiency of Ks reconstruction in (cos*theta, z) space
    saving efficiency histogram to efficiency_2D.root
    """
    hist_model = R.RDF.TH2DModel("", ";cos#theta*;z", 10, -1, 1, 15, 0.05, 0.8)

    rdf_reco_truth = R.RDataFrame("event", reco_truth_matched_rootFile )
    hist_cosTheta_z_reco = rdf_reco_truth.Histo2D(hist_model, "Ks_helicity_angle", "Ks_z").GetValue()

    rdf_truth = R.RDataFrame("truth", truth_rootFile )
    hist_cosTheta_z_truth = rdf_truth.Histo2D(hist_model, "Ks_helicity_angle", "Ks_z").GetValue()

    R.gStyle.SetOptStat(0)

    c = R.TCanvas("c_eff", "c_eff", 800, 600)
    hist_cosTheta_z_reco.Draw("COLZ TEXT")
    c.SaveAs(os.path.join(output_dir, "Ks_costheta_vs_z_reco_truth_matched.png"))

    c2 = R.TCanvas("c_eff_truth", "c_eff_truth", 800, 600)
    hist_cosTheta_z_truth.Draw("COLZ TEXT")
    c2.SaveAs(os.path.join(output_dir, "Ks_costheta_vs_z_truth.png"))

    eff = R.TEfficiency(hist_cosTheta_z_reco, hist_cosTheta_z_truth)
    eff.SetName("hist_eff")

    # create a histogram from TEfficiency for visualization
    hist_eff = hist_cosTheta_z_reco.Clone("hist_eff")
    hist_eff.SetTitle("Efficiency vs cos#theta and z;cos#theta*; z")

    for i in range(1, hist_eff.GetNbinsX() + 1):
        for j in range(1, hist_eff.GetNbinsY() + 1):
            global_bin = eff.GetGlobalBin(i, j)
            eff_value = eff.GetEfficiency(global_bin)
            eff_error_low = eff.GetEfficiencyErrorLow(global_bin)
            eff_error_high = eff.GetEfficiencyErrorUp(global_bin)
            
            # use average of asymmetric errors for symmetric error bar
            eff_error = (eff_error_low + eff_error_high) / 2.0

            hist_eff.SetBinContent(i, j, eff_value)
            hist_eff.SetBinError(i, j, eff_error)
    c3 = R.TCanvas("c_efficiency", "c_efficiency", 800, 600)
    R.gStyle.SetPaintTextFormat(".3f")
    hist_eff.Draw("COLZ TEXT")
    c3.SaveAs(os.path.join(output_dir, "Ks_efficiency_costheta_vs_z.png"))

    eff_file = R.TFile.Open(eff_root, "RECREATE")
    hist_cosTheta_z_reco.Write()
    hist_cosTheta_z_truth.Write()
    hist_eff.Write()
    eff_file.Close()
    print("Efficiency histogram saved to efficiency_2D.root")



def extract_rho00(nsig_txt:str, output_dir:str, eff_root:str = None):
    df = pd.read_csv(nsig_txt, sep=r"\s+",
                    names=["z_center", "z_width", "helicity_center", 
                            "helicity_width", "nsig", "nsig_err", "nsig_err2"],
                    skiprows=1,
                    comment="#")

    # filter invalid data
    df_valid = df[(df['nsig'] > 0) | (df['z_center'] > 0)].reset_index(drop=True)
    
    # the range i wanna use
    df_valid = df_valid[(df_valid['z_center'] >= 0.05) & (df_valid['z_center'] <= 0.75)].reset_index(drop=True)

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
        """
        # Create fit function: N(cos(theta)) = N0 * (1 - rho00 + (3*rho00 - 1)*cos^2(theta))
        angular_dist = R.TF1(f"angular_dist_{idx}", "[0]*(1-[1]+(3*[1]- 1)*x*x)", -1, 1, 2)
        
        # Set initial parameter values for better convergence
        # [0] = normalization, estimate from histogram integral
        # [1] = rho00, start from 1/3 (unpolarized case)
        hist_integral = hist.Integral()
        angular_dist.SetParameter(0, hist_integral / 2.0)  # rough normalization
        angular_dist.SetParameter(1, 0.33)  # start from unpolarized hypothesis
        
        # Set parameter limits
        angular_dist.SetParLimits(0, 0, hist.GetMaximum()*2.0)  # allow more flexibility
        angular_dist.SetParLimits(1, 0, 1)  # physical range for rho00
        
        # Set parameter names for better output
        angular_dist.SetParName(0, "N0")
        angular_dist.SetParName(1, "rho00")
        
        # Perform fit with improved options
        # S - return fit status, M - improve fit, E - better error calc, R - use function range, Q - quiet
        fit_result = hist.Fit(angular_dist, "RSMEQ")
        
        # Check fit quality and retry if needed
        if fit_result.Status() != 0 or not fit_result.IsValid():
            # Try with different initial value for rho00
            for init_rho in [0.1, 0.5, 0.7]:
                angular_dist.SetParameter(1, init_rho)
                fit_result = hist.Fit(angular_dist, "RSMEQ")
                if fit_result.Status() == 0 and fit_result.IsValid():
                    break
        
        # fit will clear hist's draw options, so draw after fit

        max = hist.GetMaximum()
        min = hist.GetMinimum()
        c = style_draw([hist], "", styles = [HistStyle.error_bars(1)], save= False, y_min = min*0.9, y_max = max*1.1 ,use_user_y_range=True)
        #c.SetCanvasSize(1000, 1000)
        c.SetCanvasSize(800, 600)
        c.SetWindowSize(800, 600)
        c.Update()
        
        rho00_value = angular_dist.GetParameter(1)
        rho00_error = angular_dist.GetParError(1)
        chi2 = angular_dist.GetChisquare()
        ndf = angular_dist.GetNDF()
        chi2_ndf = chi2 / ndf if ndf > 0 else 0
        
        z_values.append(z_value)
        rho00_values.append(rho00_value)
        rho00_errors.append(rho00_error)
        
        angular_dist.Draw("same")
        
        leg = R.TLegend(0.6, 0.7, 0.9, 0.90)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetFillColor(0)
        leg.SetTextSize(0.04)
        leg.AddEntry(0, "#rho_{00}"+ f" = {rho00_value:.3f} #pm {rho00_error:.3f} ", "")
        leg.AddEntry(0, "#chi^{2}/ndf" + f" = {chi2:.2f}/{ndf:.0f}", "")
        leg.AddEntry(0, f"{z_value - 0.025:.2f} < z < {z_value + 0.025:.2f}", "")
        leg.Draw("same")

        c.SaveAs(os.path.join(output_dir, f"fit_{idx:02d}.png"))
    """
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
        

if __name__ == "__main__":
    eff_rootFile = "./efficiency_2D.root"
    eff_original_rootFile = "./images/eff_v2.1.0/efficiency_2D.root"
    reco_truth_matched_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.0.0_qqbar/exp55_reco_truth_matched.root"
    truth_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.1.0_qqbar/exp55_truth_processed.root"
    nsig_txt = "images/fit_dt_e55_v2.1.0/nsig_results.txt"

    if not os.path.exists(eff_rootFile):
        extract_eff(reco_truth_matched_rootFile, truth_rootFile, eff_root=eff_rootFile, output_dir="./images/eff")
    if not os.path.exists(eff_original_rootFile):
        nsig_mc_txt = "images/fit_mc_e55_v2.1.0/nsig_results.txt"
        calculate_eff_original(truth_rootFile, nsig_mc_txt, output_dir="images/eff_v2.1.0/")
    
    output_dir = "./images/fit_rho00/"
    os.makedirs(output_dir, exist_ok=True)
    #sextract_rho00(nsig_txt, output_dir, eff_rootFile)
    #extract_rho00(nsig_txt, output_dir)

    output_dir_original = "./images/fit_rho00_v2.1.0/"
    os.makedirs(output_dir_original, exist_ok=True)
    #extract_rho00(nsig_txt, output_dir_original, eff_root=eff_original_rootFile)    
    #extract_rho00(nsig_txt, output_dir_original)    
    
    #get_2d_hist_from_txt(nsig_txt)

    truth_level_rho00(truth_rootFile, output_dir="./test/fit_rho00_truth/")

