import ROOT as R
from DRAW import style_draw, HistStyle

def main():
    rootFile =  "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar/exp55_reco_truth_matched.root"
    df = R.RDataFrame("event", rootFile)
    condition = "Ks_z >= 0.6 && Ks_z < 0.65 && Ks_helicity_angle >= -0.2 && Ks_helicity_angle < 0"
    
    #variables = [""]
    df = df.Redefine("Ks_M", f"Ks_M[{condition}]")
    df = df.Filter(f"Any({condition})")

    model = R.RDF.TH1DModel("", ";p(#pi^{+});[MeV]" ,90, 0, 4.5)
    h1 = df.Histo1D(model, "pip_p")
    h2 = df.Redefine("pip_p", f"pip_p[{condition}]").Histo1D(model, "pip_p")
    print(type(h1), type(h2))

    style_draw([h1, h2], "test.png", leg_texts=["before", "after"], styles=[HistStyle.line_hist(4), HistStyle.error_bars(1)])

    df.Report().Print()

if __name__ == "__main__":
    main()