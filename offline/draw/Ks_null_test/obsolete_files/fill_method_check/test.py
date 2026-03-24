import os
import ROOT as R
from DRAW import style_draw, HistStyle
from FIT import fit_rho00


def extract_eff(reco_truth_matched_rootFile, truth_rootFile, eff_root="efficiency_2D.root", output_dir = "."):
    """
    calculate 2D efficiency of Ks reconstruction in (cos*theta, z) space
    saving efficiency histogram to efficiency_2D.root
    """
    hist_model = R.RDF.TH2DModel("", ";cos#theta*;z", 10, -1, 1, 15, 0.05, 0.8)

    rdf_reco_truth = R.RDataFrame("event", reco_truth_matched_rootFile )
    hist_cosTheta_z_reco = rdf_reco_truth.Histo2D(hist_model, "Ks_helicity_angle", "Ks_z").GetValue()

    rdf_truth = R.RDataFrame("truth", truth_rootFile )
    hist_cosTheta_z_truth = rdf_truth.Histo2D(hist_model, "Ks_helicity_angle", "Ks_z").GetValue()

    R.gStyle.SetOptStat(0)

    c = R.TCanvas("c_eff", "c_eff", 800, 600)
    hist_cosTheta_z_reco.Draw("COLZ TEXT")
    c.SaveAs(os.path.join(output_dir, "Ks_costheta_vs_z_reco_truth_matched.png"))

    c2 = R.TCanvas("c_eff_truth", "c_eff_truth", 800, 600)
    hist_cosTheta_z_truth.Draw("COLZ TEXT")
    c2.SaveAs(os.path.join(output_dir, "Ks_costheta_vs_z_truth.png"))

    eff = R.TEfficiency(hist_cosTheta_z_reco, hist_cosTheta_z_truth)
    eff.SetName("hist_eff")

    # create a histogram from TEfficiency for visualization
    hist_eff = hist_cosTheta_z_reco.Clone("hist_eff")
    hist_eff.SetTitle("Efficiency vs cos#theta and z;cos#theta*; z")

    for i in range(1, hist_eff.GetNbinsX() + 1):
        for j in range(1, hist_eff.GetNbinsY() + 1):
            global_bin = eff.GetGlobalBin(i, j)
            eff_value = eff.GetEfficiency(global_bin)
            eff_error_low = eff.GetEfficiencyErrorLow(global_bin)
            eff_error_high = eff.GetEfficiencyErrorUp(global_bin)
            
            # use average of asymmetric errors for symmetric error bar
            eff_error = (eff_error_low + eff_error_high) / 2.0

            hist_eff.SetBinContent(i, j, eff_value)
            hist_eff.SetBinError(i, j, eff_error)
    c3 = R.TCanvas("c_efficiency", "c_efficiency", 800, 600)
    R.gStyle.SetPaintTextFormat(".3f")
    hist_eff.Draw("COLZ TEXT")
    c3.SaveAs(os.path.join(output_dir, "Ks_efficiency_costheta_vs_z.png"))

    eff_file = R.TFile.Open(eff_root, "RECREATE")
    hist_cosTheta_z_reco.Write()
    hist_cosTheta_z_truth.Write()
    hist_eff.Write()
    eff_file.Close()
    print("Efficiency histogram saved to efficiency_2D.root")


def test_fill_method():
    truth_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar_svd2/svd2_truth_processed.root"
    rdf = R.RDataFrame("truth", truth_rootFile)

    rdf = rdf.Filter("Ks_M.size() > 1", "has more than 1 Ks")
    rdf.Report().Print()

    hist_model = R.RDF.TH2DModel("", ";cos#theta*;z", 10, -1, 1, 15, 0.05, 0.8)
    hist_cosTheta_z = rdf.Histo2D(hist_model , "Ks_helicity_angle", "Ks_z")
    c = R.TCanvas("c1", "c1", 800, 600)
    R.gStyle.SetPaintTextFormat(".0f")
    hist_cosTheta_z.Draw("COLZ TEXT")
    c.SaveAs("test_fill_method_before.png")

    return 
    
    rdf = rdf.Define("null_test", "Ks_M.size()/Ks_M.size()")
    rdf = rdf.Define("test2", "rdfentry_/rdfentry_")

    total_entry = rdf.Count().GetValue()
    total_Ks = rdf.Define("n", "Ks_M.size()").Sum("n").GetValue()

    hist_model = R.RDF.TH2DModel("", ";cos#theta*;z", 10, -1, 1, 15, 0.05, 0.8)
    hist_cosTheta_z = rdf.Histo2D(hist_model, "Ks_helicity_angle", "Ks_z", "test2").GetValue() # if it is a list, it will count the number of element not event number
    h2 = rdf.Histo2D(hist_model, "Ks_helicity_angle", "Ks_z", "null_test").GetValue()

    sum1 = 0
    sum2 = 0
    for i in range(1, hist_cosTheta_z.GetNbinsX() + 1):
        for j in range(1, hist_cosTheta_z.GetNbinsY() + 1):
            sum1 += hist_cosTheta_z.GetBinContent(i, j)
            sum2 += h2.GetBinContent(i, j)
            print("Bin ({0}, {1}): hist_cosTheta_z = {2}, h2 = {3}".format(
                i, j,
                hist_cosTheta_z.GetBinContent(i, j),
                h2.GetBinContent(i, j)
            ))


    print("Total entries:", total_entry)
    print("Total Ks:", total_Ks)
    print("Sum of hist_cosTheta_z bin contents:", sum1)
    print("Sum of h2 bin contents:", sum2)

    return 
    # fill each bin with Ks_M.size()
    hist_cosTheta_z_filled = hist_cosTheta_z.Clone("hist_cosTheta_z_filled")
    hist_cosTheta_z_filled.Reset()

    for i in range(1, hist_cosTheta_z.GetNbinsX() + 1):
        for j in range(1, hist_cosTheta_z.GetNbinsY() + 1):
            bin_condition =  "(Ks_helicity_angle > {0} && Ks_helicity_angle <= {1} && Ks_z > {2} && Ks_z <= {3})".format(
                hist_cosTheta_z.GetXaxis().GetBinLowEdge(i),
                hist_cosTheta_z.GetXaxis().GetBinUpEdge(i),
                hist_cosTheta_z.GetYaxis().GetBinLowEdge(j), 
                hist_cosTheta_z.GetYaxis().GetBinUpEdge(j)
            )
            print(i, j, bin_condition)
            
            rdf_ij = rdf
            for var in ["Ks_M", "Ks_weight"]:
                if var in rdf_ij.GetColumnNames():
                    rdf_ij = rdf_ij.Redefine(var, f"{var}[{bin_condition}]")
            rdf_ij.Filter(f"Any({bin_condition})")
            
            #nKs =  rdf_ij.Sum("Ks_M * Ks_weight").GetVaule()
            rdf_ij = rdf_ij.Define("Ks_M_size", "Ks_M.size()")
            nKs =  rdf_ij.Sum("Ks_M_size").GetValue()
            nEvents = hist_cosTheta_z.GetBinContent(i, j)
            print("Bin ({0}, {1}): Filtered value = {2}".format(i, j, nKs), "Original value = {0}".format(nEvents))
            hist_cosTheta_z_filled.SetBinContent(i, j, nKs)

if __name__ == "__main__":
    test_fill_method()