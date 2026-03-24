import ROOT as R
import pandas as pd
import os

def calculate_eff(truth_rootFile, nsig_txt, output_dir):
    """
    calculate 2D efficiency of phi reconstruction in (cos*theta, z) space
    saving efficiency histogram to efficiency_2D.root
    """
    # truth info stored in ./truth_processed.root , processed by function phi_truth in ../../Process.py
    rdf_truth = R.RDataFrame("truth", truth_rootFile )
    hist_cosTheta_z = rdf_truth.Histo2D(("hist_cosTheta_z", "cos#theta vs phi z;cos#theta*;phi z", 10, -1, 1, 17, 0.15, 1), "phi_helicity_angle", "phi_z")

    R.gStyle.SetOptStat(0)

    c = R.TCanvas("c_eff", "c_eff", 800, 600)
    hist_cosTheta_z.Draw("COLZ TEXT")
    c.SaveAs(os.path.join(output_dir, "phi_costheta_vs_z.png"))
    
    # read nsig extracted by fitting frin ./MC_fitting/nsig_results.txt
    df_nsig = pd.read_csv(nsig_txt, sep=r"\s+",
                        names=["z_center", "z_width", "helicity_center", 
                               "helicity_width", "nsig", "nsig_err", "nsig_err2"],
                        skiprows=1)
    
    df_nsig = df_nsig[(df_nsig['nsig'] > 0)]

    hist_nsig = R.TH2D("hist_nsig", "nsig in cos#theta vs phi z;cos#theta*;phi z", 10, -1, 1, 17, 0.15, 1)

    for index, row in df_nsig.iterrows():
        binx = hist_nsig.GetXaxis().FindBin(row['helicity_center'])
        biny = hist_nsig.GetYaxis().FindBin(row['z_center'])

        hist_nsig.SetBinContent(binx, biny, row['nsig'])
        hist_nsig.SetBinError(binx, biny, row['nsig_err'])
    
    c2 = R.TCanvas("c_eff_nsig", "c_eff_nsig", 1600, 1080)
    hist_nsig.Draw("COLZ TEXT")
    c2.SaveAs(os.path.join(output_dir, "phi_nsig_costheta_vs_z.png"))

    # Use TEfficiency for correct error calculation
    # TEfficiency takes (passed, total) histograms
    eff = R.TEfficiency(hist_nsig, hist_cosTheta_z.GetPtr())
    eff.SetName("eff_costheta_z")
    eff.SetTitle("Efficiency vs cos#theta and z;cos#theta;phi z;Efficiency")
    
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
    c3.SaveAs(os.path.join(output_dir, "phi_efficiency_costheta_vs_z.png"))
    
    # Also save the TEfficiency object to a ROOT file for later use
    output_file = R.TFile(os.path.join(output_dir ,"efficiency_2D.root"), "RECREATE")
    eff.Write()
    hist_eff.Write()
    hist_cosTheta_z.GetPtr().Write()
    hist_nsig.Write()
    output_file.Close()
    
    return hist_eff, eff


if __name__ == "__main__":
    output_dir = "./new_images/"
    truth_rootFile = "/gpfs/group/belle2/users2022/luruihua/for_wangz/data_gMC_belle1/2025-11-12_SpinAlignment_gMC/continuum_truth_processed.root"
    nsig_txt = "./new_images/fitting/nsig_results.txt"
    calculate_eff(truth_rootFile , nsig_txt, output_dir)

