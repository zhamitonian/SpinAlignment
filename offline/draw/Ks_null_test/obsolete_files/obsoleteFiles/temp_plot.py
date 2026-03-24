from DRAW import style_draw, HistStyle
import ROOT

path1 = "temp/temp_bin_39.root"
path2 = "images/fit/temp_bin_39.root"

def quick_plot(path, output):
    df = ROOT.RDataFrame("event", path)
    hist = df.Histo1D(("", ";M_{#pi^{+}#pi^{-}};[MeV]", 110, 0.47, 0.525), "Ks_M").GetValue()

    double_gaus = ROOT.TF1("double_gaus", "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [1], [4])", hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
    double_gaus.SetParLimits(0, 0, hist.GetMaximum())
    double_gaus.SetParLimits(3, 0, hist.GetMaximum())
    double_gaus.SetParLimits(1, 0.4, 0.6)
    double_gaus.SetParLimits(2, 0.000001, 0.0001)  
    double_gaus.SetParLimits(4, 0.0001, 0.01)    
    double_gaus.SetParLimits(2, 0.01, 0.05)  

    hist.Fit(double_gaus, "RQ")
    c = style_draw([hist], output, styles=[HistStyle.line_hist(1)], save=False)
    double_gaus.Draw("SAME")
    c.SaveAs(output)
    para = [double_gaus.GetParameter(i) for i in range(double_gaus.GetNpar())]
    for i, p in enumerate(para):
        print(f"par[{i}] = {p}")

quick_plot(path1, "temp/quick_plot_39.png")
quick_plot(path2, "images/fit/quick_plot_39.png")