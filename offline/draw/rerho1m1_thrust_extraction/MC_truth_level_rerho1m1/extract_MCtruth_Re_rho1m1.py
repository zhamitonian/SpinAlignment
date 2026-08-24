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
PHI_BIN_NAME = "phi_thrust_helicity_phi"

def extract_rerho1m1(MCtruth_path:str, output_dir:str):
    df = R.RDataFrame("truth", MCtruth_path)

    z_values = []
    rho00_values = []
    rho00_errors = []

    for z in range(Z_BIN_CONFIG[0]):
        z_width = (Z_BIN_CONFIG[2] - Z_BIN_CONFIG[1]) / Z_BIN_CONFIG[0]
        zmin = Z_BIN_CONFIG[1] + z * z_width
        zmax = zmin + z_width

        hist = df.Redefine(PHI_BIN_NAME, f"{PHI_BIN_NAME}[{Z_BIN_NAME} >= {zmin} && {Z_BIN_NAME} < {zmax}]") \
                 .Histo1D(("", ";#phi^{#ast}(rad);candidate", 20, 0, math.pi/2), PHI_BIN_NAME).GetValue()

        value, err = fit_rerho1m1(hist, os.path.join(output_dir, f"fit_{z:02d}.png"), True, extra_legend=f"{zmin:.2f} < z < {zmax:.2f}")
        z_values.append((zmin + zmax) / 2)
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
    truth_rootFile = "../truth_output_qqbar.root"
    output_dir = "../images/truth_rerho1m1/"
    os.makedirs(output_dir, exist_ok=True)

    extract_rerho1m1(truth_rootFile, output_dir)


