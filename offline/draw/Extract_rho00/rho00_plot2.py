import ROOT as R
from DRAW import style_draw, HistStyle

import pandas as pd

def main():
    Kstar_rho00_csv = "../Kstar892_rho00extraction/images/rho00/data/rho00_results.csv"
    theory_csv = "../Kstar892_rho00extraction/rho00_1052MeV.csv"

    df_Kstar = pd.read_csv(Kstar_rho00_csv)
    df_theory = pd.read_csv(theory_csv)

    def get_hist(df):
        hist = R.TH1D("", ";z;#rho_{00}", 20, 0, 1)
        for i, row in df.iterrows():
            z = row["z_center"]
            rho00 = row["rho00"]
            rho00_err = row["rho00_error"]

            bin_idx = hist.FindBin(z)
            hist.SetBinContent(bin_idx, rho00)
            hist.SetBinError(bin_idx, rho00_err)
        return hist

    hist_Kstar = get_hist(df_Kstar)


    def get_theory_hist(df):
        hist = R.TH1D("", ";z;#rho_{00}", 100, 0.005, 1.005)
        for i, row in df.iterrows():
            z = row["z"]
            rho00 = row["rho00"]

            bin_idx = hist.FindBin(z)
            hist.SetBinContent(bin_idx, rho00)
            hist.SetBinError(bin_idx, 0)  # No error for theory points
        return hist

    hist_theory = get_theory_hist(df_theory)

    c = style_draw([hist_theory, hist_Kstar], "", ["Theory calculation", "K* #rho_{00}"],
                   styles=[HistStyle.error_bars(R.TColor.GetColor("#D55E00")), HistStyle.error_bars(R.TColor.GetColor("#0072B2"))],
                   legend_position= 2 , 
                   save=False, y_min=0, y_max=0.8, use_user_y_range=True)

    line = R.TLine(0, 1/3, 1, 1/3)
    #line.SetLineColor(R.TColor.GetColor("#CC79A7"))
    line.SetLineColor(2)
    line.SetLineStyle(R.kDashed)
    line.SetLineWidth(4)
    line.Draw("same")

    leg = R.TLegend(0.6, 0.2, 0.9, 0.45)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(line, "#rho_{00} = 1/3(unpolarized)", "l")

    leg.Draw("same")
    c.SaveAs("images/rho00_theory_vs_kstar.png")

if __name__ == "__main__":
    main()


            

    
    