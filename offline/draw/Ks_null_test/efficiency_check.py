import ROOT as R
from DRAW import style_draw, HistStyle
import pandas as pd
import os

def get_hist_from_csv(csv_file):
    df = pd.read_csv(csv_file)
    hist = R.TH2D("", ";cos#theta;z", 10, -1, 1, 15, 0.05, 0.8)

    for index, row in df.iterrows():
        binx = hist.GetXaxis().FindBin(row["Ks_helicity_angle_center"])
        biny = hist.GetYaxis().FindBin(row["Ks_z_center"])

        hist.SetBinContent(binx, biny, row['nsig'])
        hist.SetBinError(binx, biny, row['nsig_err_hi'])

    return hist

def get_2d_hist_from_txt(nsig_txt):
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

    R.gStyle.SetOptStat(0)
    c2 = R.TCanvas("c_eff_nsig", "c_eff_nsig", 1600, 1080)
    hist_nsig.Draw("COLZ TEXT")
    c2.SaveAs(os.path.join("./", "Ks_nsig_costheta_vs_z.png"))

    return hist_nsig

def plot_2d_hist(hist, output_path):
    R.gStyle.SetOptStat(0)
    c = R.TCanvas("c_eff", "c_eff", 1600, 1080)
    R.gStyle.SetPaintTextFormat(".1f")
    hist.Draw("COLZ TEXT")
    c.SaveAs(output_path)

def plot_TEfficiency(eff, output_path):
    c = R.TCanvas("c_efficiency", "c_efficiency", 1600, 1080)
    R.gStyle.SetPaintTextFormat(".3f")
    eff.SetTitle("Efficiency vs cos#theta and z;cos#theta;Ks z;Efficiency")
    eff.Draw("COLZ TEXT")
    c.SaveAs(output_path)

def main():
    """
    truth_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.3_qqbar_svd2/svd2_only/svd2_st0_truth_processed.root"
    reco_rootFile = "./rootFiles/svd2_only_v2.3.3/sig_isSignal_v2.3.3.root"
    fit_csv = "./fit_results/MC_fit/nsig_results.csv"

    hist_model = R.RDF.TH2DModel("hist_eff", ";cos#theta;z", 10, -1, 1, 15, 0.05, 0.8)

    rdf_truth = R.RDataFrame("truth", truth_rootFile)
    hist_truth = rdf_truth.Histo2D(hist_model, "Ks_helicity_angle", "Ks_z").GetValue()
    rdf_reco = R.RDataFrame("event", reco_rootFile)
    hist_reco = rdf_reco.Histo2D(hist_model, "Ks_helicity_angle", "Ks_z", "Ks_weight").GetValue()

    hist_fit = get_hist_from_csv(fit_csv)

    hist_diff = hist_fit.Clone("hist_diff")
    hist_diff.Add(hist_reco, -1)

    plot_2d_hist(hist_truth, "./images/nsig_truth.png")
    plot_2d_hist(hist_fit, "./images/nsig_fit.png")
    plot_2d_hist(hist_reco, "./images/nsig_reco.png")
    plot_2d_hist(hist_diff, "./images/nsig_diff.png")

    eff_MCtruth = R.TEfficiency(hist_reco, hist_truth)
    eff_fit = R.TEfficiency(hist_fit, hist_truth)

    plot_TEfficiency(eff_MCtruth, "./images/efficiency_MC_truth.png")
    plot_TEfficiency(eff_fit, "./images/efficiency_fit.png")
    """

    hist_model = R.RDF.TH2DModel("hist_eff", ";cos#theta;z", 10, -1, 1, 15, 0.05, 0.8)

    #------------- comparing to previous results -------------
    truth_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar_svd2/svd2_st0_truth_processed.root"
    fit_txt = "/gpfs/group/belle/users/wangz/Work/SpinAlignment/offline/draw/Ks_null_test/v2.3.2_svd2/qqbar/output_st0/nsig_results.txt"
    #hist_fit_prev = get_2d_hist_from_txt(fit_txt)
    rdf_truth = R.RDataFrame("truth", truth_rootFile)
    """
    hist_truth = rdf_truth.Histo2D(hist_model, "Ks_helicity_angle", "Ks_z").GetValue()
    eff_fit_prev = R.TEfficiency(hist_fit_prev, hist_truth)
    plot_TEfficiency(eff_fit_prev, "./images/efficiency_fit_prev.png")
    plot_2d_hist(hist_fit_prev, "./images/nsig_fit_prev.png")
    plot_2d_hist(hist_truth, "./images/nsig_truth_prev.png")

    hist_eff = hist_fit_prev.Clone("hist_eff")
    hist_eff.SetTitle("Efficiency vs cos#theta and z;cos#theta*; z")
    
    # Fill efficiency histogram with correct values and errors
    for i in range(1, hist_eff.GetNbinsX() + 1):
        for j in range(1, hist_eff.GetNbinsY() + 1):
            global_bin = eff_fit_prev.GetGlobalBin(i, j)
            eff_value = eff_fit_prev.GetEfficiency(global_bin)
            eff_error_low = eff_fit_prev.GetEfficiencyErrorLow(global_bin)
            eff_error_high = eff_fit_prev.GetEfficiencyErrorUp(global_bin)
            
            # Use average of asymmetric errors for symmetric error bar
            eff_error = (eff_error_low + eff_error_high) / 2.0
            
            hist_eff.SetBinContent(i, j, eff_value)
            hist_eff.SetBinError(i, j, eff_error)
    plot_TEfficiency(eff_fit_prev, "./images/efficiency_fit_prev_2.png")
    """
    rdf_truth = rdf_truth.Redefine("Ks_M", "0.5")
    hist_truth = rdf_truth.Histo2D(hist_model, "Ks_helicity_angle", "Ks_z").GetValue()
    hist_truth2 = rdf_truth.Histo2D(hist_model, "Ks_helicity_angle", "Ks_z", "Ks_M").GetValue() #less 1/2 ?
    hist_diff = hist_truth2.Clone("hist_diff")
    hist_diff.Add(hist_truth, -1)
    plot_2d_hist(hist_truth, "./images/nsig_truth_prev.png")
    plot_2d_hist(hist_truth2, "./images/nsig_truth_prev_M.png") 
    plot_2d_hist(hist_diff, "./images/nsig_truth_diff.png")


    


if __name__ == "__main__":
    main()
    