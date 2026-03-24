## ------------------------------
## perform fitting to primary vertex 's Z position distribution
## to extraction it's nominal position for each experiment
## ------------------------------

import sys, os
import ROOT as R
from DRAW import style_draw, HistStyle


def perform_fitting(hist, output_dir , tag):
    """
    perform fitting to primary vertex 's Z position distribution
    to extraction it's nominal position
    """
    gaus = R.TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", -3.5, 3.5)
    gaus.SetParameters(hist.GetMaximum(), 0, 0.5)
    hist.Fit(gaus, "RQ")
    mean = gaus.GetParameter(1)
    sigma = gaus.GetParameter(2)

    ## plotting
    c = style_draw([hist], "./test.png", styles=[HistStyle.line_hist(1)], save=False)
    gaus.Draw("SAME")
    leg = R.TLegend(0.73, 0.75, 0.93, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.AddEntry(hist, "Primary Vertex Z Position", "l")
    leg.AddEntry(gaus, "Gaussian Fit", "l")
    leg.AddEntry(0, f"Mean = {mean:.4f} cm", "")
    leg.AddEntry(0, f"Sigma = {sigma:.4f} cm", "")
    leg.Draw("same")

    c.SaveAs(os.path.join(output_dir, f"primeVz_fitting_e{tag}.png"))
    return mean, sigma

def forEach_exp():
    """
    loop over each experiment's root file to extract center prime Vz
    """
    output_dir = "./images/primeVz_fitting"
    os.makedirs(output_dir, exist_ok=True)

    # record results in a txt file
    txt_file = os.path.join(output_dir, "center_primeVz_byExp.txt")  
    with open(txt_file, "w") as f_out:
        f_out.write("# exp_id\tmean\t sigma\n")

    for rootFile_name in os.listdir("./hist_roots"):
        exp_id = rootFile_name.replace("exp", "").replace(".root", "")

        rootFile_path = os.path.join("./hist_roots", rootFile_name)
        rootFile = R.TFile.Open(rootFile_path)
        hist_vz = rootFile.Get("PrimeVz_data").Clone()
        hist_vz.SetDirectory(0)  # Detach from file
        rootFile.Close()

        mean, sigma = perform_fitting(hist_vz, output_dir, exp_id)

        with open(txt_file, "a") as f_out:
            f_out.write(f"{exp_id}\t{mean:.6f}\t{sigma:.6f}\n")

forEach_exp()

    

    

