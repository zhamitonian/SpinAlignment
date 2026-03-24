import ROOT as R
from DRAW import style_draw, HistStyle

def main(b2bii_path, basf_path):
    file_b2bii = R.TFile(b2bii_path, "READ")
    track_tree = file_b2bii.Get("track")
    track_tree.AddFriend("photon", b2bii_path)
    df_b2bii = R.RDataFrame(track_tree)
    df_basf = R.RDataFrame("event", basf_path)

    df_b2bii = df_b2bii.Define("sum_chg", "sum(sqrt(p_CMS*p_CMS + 0.13957*0.13957))")
    df_b2bii = df_b2bii.Define("sum_gam", "sum(photon.p_CMS)")
    df_b2bii = df_b2bii.Define("Evis", "sum_chg + sum_gam")

    df_basf = df_basf.Define("Evis", "sum(sqrt(trk_p_CMS*trk_p_CMS + 0.13957*0.13957) + gam_p_CMS*)")

    hist_b2bii = df_b2bii.Histo1D(("", ";E_{vis};[MeV]", 50, 0, 15), "Evis")
    hist_basf = df_basf.Histo1D(("", ";E_{vis};[MeV]", 50, 0, 15), "Evis")
    style_draw([hist_b2bii, hist_basf], "Evis_comparison.png", ["b2bii", "basf"], styles = [HistStyle.error_bars(1), HistStyle.line_hist(4)], log_y= True)

if __name__ == "__main__":
