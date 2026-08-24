import pandas as pd
import ROOT as R
from DRAW import style_draw, HistStyle

orignial_csv = "../Extract_rho00/images/rho00/data_bin10/rho00_results.csv"
ME_correction_csv = "images/rho00/data_bin10/rho00_results.csv"
no_correction_csv = "images/rho00/data_bin10_no_correction/rho00_results.csv"

def get_hist_from_csv(csv_file):
    df = pd.read_csv(csv_file)
    hist = R.TH1D("hist", ";z;#rho_{00}", 20, 0, 1)
    for idx, row in df.iterrows():
        z_value = row["z_center"]
        rho00 = row["rho00"]
        rho00_err = row["rho00_error"]
        bin_idx = hist.FindBin(z_value)
        hist.SetBinContent(bin_idx, rho00)
        hist.SetBinError(bin_idx, rho00_err)
    return hist 

orignial_hist = get_hist_from_csv(orignial_csv)
ME_correction_hist = get_hist_from_csv(ME_correction_csv)
no_correction_hist = get_hist_from_csv(no_correction_csv)

c = style_draw([orignial_hist, ME_correction_hist, no_correction_hist], "rho00_comparison.png",
           ["MC eff corrected", "Mixed Event corrected", "No correction"],
           [HistStyle.error_bars(4), HistStyle.error_bars(8), HistStyle.error_bars(6)], save=False, y_min=0, y_max=0.8, use_user_y_range=True)

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
c.SaveAs("./rho00_comparison.png")