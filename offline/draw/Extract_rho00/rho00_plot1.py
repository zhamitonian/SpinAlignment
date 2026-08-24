import ROOT as R
from DRAW import style_draw, HistStyle

import pandas as pd
import numpy as np

M_KSTAR = 0.896 # GeV 
# @ LEP , opal measurement 10.1016/0031-9163(66)90451-3
E_CMS = 91.2 # GeV, Z^0 mass 
P_BEAM = E_CMS / 2 # GeV


def main():
    phi_bin8_rho00_csv = "./images/rho00/data_bin8/rho00_results.csv"
    Kstar_rho00_csv = "../Kstar892_rho00extraction/images/rho00/data/rho00_results.csv"
    OPAL_Kstar_rho00_csv = "./OPAL_Kstar.csv"

    df_phi_bin8 = pd.read_csv(phi_bin8_rho00_csv)
    df_Kstar = pd.read_csv(Kstar_rho00_csv)

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

    hist_phi_bin8 = get_hist(df_phi_bin8)
    hist_Kstar = get_hist(df_Kstar)


    # ---- OPAL Kstar data ----
    df_OPAL_Kstar = pd.read_csv(OPAL_Kstar_rho00_csv)

    df_OPAL_Kstar['z'] = np.sqrt((df_OPAL_Kstar['xp'] * P_BEAM) **2 + M_KSTAR**2) / P_BEAM
    df_OPAL_Kstar['z_low'] = np.sqrt((df_OPAL_Kstar['xp_low'] * P_BEAM) **2 + M_KSTAR**2) / P_BEAM
    df_OPAL_Kstar['z_high'] = np.sqrt((df_OPAL_Kstar['xp_high'] * P_BEAM) **2 + M_KSTAR**2) / P_BEAM
    df_OPAL_Kstar['rho00_error'] = np.sqrt(df_OPAL_Kstar['stat_plus']**2 + df_OPAL_Kstar['sys_plus']**2)

    bin_edges = list(df_OPAL_Kstar['z_low']) + [df_OPAL_Kstar['z_high'].iloc[-1]]
    hist_opal = R.TH1D("", ";z;#rho_{00}", len(bin_edges)-1, np.array(bin_edges, dtype=np.float64))

    for i, row in df_OPAL_Kstar.reset_index().iterrows():
        z = row['z']
        rho00 = row['rho00']
        rho00_err = row['rho00_error']
        bin_idx = hist_opal.FindBin(z)
        hist_opal.SetBinContent(bin_idx, rho00)
        hist_opal.SetBinError(bin_idx, rho00_err)

    
    c = style_draw([hist_phi_bin8, hist_Kstar, hist_opal], "", ["#phi", "K^{*0}", "K^{*0} @ OPAL"],
                   styles=[HistStyle.error_bars(R.TColor.GetColor("#E69F00")), HistStyle.error_bars(R.TColor.GetColor("#0072B2")), HistStyle.error_bars(R.TColor.GetColor("#009E73"))],
                   legend_position= 0 , 
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


    c.SaveAs("images/rho00_vs_z_1.png")

if __name__ == "__main__":
    main()


            

    
    