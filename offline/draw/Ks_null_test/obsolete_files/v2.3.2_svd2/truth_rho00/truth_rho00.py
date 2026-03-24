import ROOT as R
from FIT import fit_rho00
import os 
from DRAW import style_draw, HistStyle


def truth_level_rho00(truth_rootFile,output_dir):
    os.makedirs(output_dir, exist_ok=True)

    rdf = R.RDataFrame("truth", truth_rootFile )
    z_values = []
    rho00_values = []
    rho00_errors = []

    hist_model = R.RDF.TH1DModel("", ";cos#theta^{*};N_{sig}", 10, -1, 1)
    for i in range(15):
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

if __name__ == "__main__":
    truth_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar_svd2/svd2_truth_processed.root"
    output_dir = "./output/"
    truth_level_rho00(truth_rootFile, output_dir)

