import ROOT as R
from DRAW import style_draw, HistStyle

def check1():
    R.EnableImplicitMT()
    reco_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.0.0_qqbar/exp55_reco_processed.root"
    reco_truth_matched_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.0.0_qqbar/exp55_reco_truth_matched.root"

    df_reco = R.RDataFrame("event", reco_rootFile )
    df_truth_matched = R.RDataFrame("event", reco_truth_matched_rootFile )


    condition = "Ks_z > 0.05 && Ks_z < 0.1"
    df_reco = df_reco.Redefine("Ks_helicity_angle", f"Ks_helicity_angle[{condition}]")
    df_truth_matched = df_truth_matched.Redefine("Ks_helicity_angle", f"Ks_helicity_angle[{condition}]")

    hist_model = R.RDF.TH1DModel("", ";cos#theta;[1]", 10, -1, 1)
    hist_reco = df_reco.Histo1D(hist_model , "Ks_helicity_angle")
    hist_truth_matched = df_truth_matched.Histo1D(hist_model , "Ks_helicity_angle")
    print(hist_reco.GetBinContent(1))

    style_draw([hist_reco, hist_truth_matched], "test.png", ["reco", "truth matched"],
           [HistStyle.line_hist(4), HistStyle.error_bars(1)])
    
def check2():
    df = R.RDataFrame("event", "./images/MC_fit/temp_bin_10.root")
    entries =  df.Count().GetValue()
    print("Entries:", entries)
    hist = df.Histo1D(("hist", "test;Ks_M;Entries", 110, 0.47, 0.525), "Ks_M").GetValue()
    hist2 = df.Histo1D(("hist", "test;Ks_M;Entries", 1000, 0., 10), "Ks_M").GetValue()
    print("Hist entries:", hist.GetEntries())
    print("Hist entries:", hist2.GetEntries())


def check3():
    df = R.RDataFrame("event", "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.0.0_qqbar/exp55_reco_processed.root")
    hist_model = R.RDF.TH2DModel("", ";cos#theta*;z", 10, -1, 1, 15, 0.05, 0.8)

    hist = df.Histo2D(hist_model , "Ks_helicity_angle", "Ks_z").GetValue()
    c = R.TCanvas("c_eff", "c_eff", 800, 600)
    #R.gStyle.SetPaintTextFormat(".3f")
    hist.Draw("COLZ TEXT")
    c.SaveAs("test2D.png")


#check1()
#check2()
check3()
     