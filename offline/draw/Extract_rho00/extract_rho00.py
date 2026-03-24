import pandas as pd 
import os 
import ROOT as R
from DRAW import style_draw, HistStyle

def extract_rho00(nsig_txt:str, output_dir:str, eff_root:str = None):
    df = pd.read_csv(nsig_txt, sep=r"\s+",
                    names=["z_center", "z_width", "helicity_center", 
                            "helicity_width", "nsig", "nsig_err", "nsig_err2"],
                    skiprows=1)

    # filter invalid data
    df_valid = df[(df['nsig'] > 0) | (df['z_center'] > 0)]

    grouped = df_valid.groupby('z_center')

    unique_z = sorted(df_valid['z_center'].unique())
    unique_z = [x for x in unique_z if x > 0]

    # get efficiency from root file
    if eff_root is not None:
        eff_file = R.TFile.Open(eff_root, "READ")
        eff_hist = eff_file.Get("hist_eff")

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
                    corrected_nsig = nsig
                    corrected_err = nsig_err
            else:
                # No efficiency correction
                corrected_nsig = nsig
                corrected_err = nsig_err
            
            bin_idx = hist.FindBin(cos_theta)
            hist.SetBinContent(bin_idx, corrected_nsig)
            hist.SetBinError(bin_idx, corrected_err)

        angular_dist = R.TF1(f"angular_dist_{idx}", "[0]*(1-[1]+(3*[1]- 1)*x*x)", -1, 1, 2)
        angular_dist.SetParLimits(0, 0, hist.GetMaximum()*1.5) 
        angular_dist.SetParLimits(1, 0, 1)  
        hist.Fit(angular_dist, "RQ")
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

        c.SaveAs(os.path.join(output_dir, f"fit_{idx}.png"))

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
    nsig_txt = "./images_test/fit_phi/nsig_results.txt"
    output_dir = "./images_test/fit_rho00/"
    eff_root = "../Costheta_z_efficiency/efficiency_2D.root"
    os.makedirs(output_dir, exist_ok=True)
    extract_rho00(nsig_txt, output_dir, eff_root)



