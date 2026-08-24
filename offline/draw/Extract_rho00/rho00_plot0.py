import ROOT as R
from DRAW import style_draw, HistStyle

import pandas as pd

def main():
    phi_rho00_csv = "/home/belle2/wangz/Work/SpinAlignment/offline/draw/Extract_rho00/images/rho00/data_bin10/rho00_results.csv"
    Ks_rho00_csv = "/home/belle2/wangz/Work/SpinAlignment/offline/draw/Ks_null_test/images/rho00/data_bin10/rho00_results.csv"

    df_phi = pd.read_csv(phi_rho00_csv)
    df_Ks = pd.read_csv(Ks_rho00_csv)

    def get_hist(df):
        hist = R.TH1D("", ";z;#rho_{00}", 20, 0, 1)
        for i, row in df.iterrows():
            z = row["#z_center"]
            rho00 = row["rho00"]
            rho00_err = row["rho00_error"]

            bin_idx = hist.FindBin(z)
            hist.SetBinContent(bin_idx, rho00)
            hist.SetBinError(bin_idx, rho00_err)
        return hist

    hist_phi = get_hist(df_phi)
    hist_Ks = get_hist(df_Ks)
    
    c = style_draw([hist_phi, hist_Ks], "", ["#phi #rho_{00}", "K_{S}^{0} #rho_{00}"],
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
    c.SaveAs("images/rho00_vs_z.png")

if __name__ == "__main__":
    main()


            

    
    