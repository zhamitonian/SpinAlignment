import ROOT as R
from DRAW import style_draw, HistStyle


def main(rootFile):
    rdf = R.RDataFrame("event", rootFile)

    def plot(var, nbin, xmin, xmax, label):
        model =  R.RDF.TH1DModel("", label, nbin, xmin, xmax)
        h1 = rdf.Histo1D(model, f"gam_{var}")
        h2 = rdf.Histo1D(model, f"gam2_{var}")
        style_draw([h1, h2], f"{var}.png", ["evtcls_hadron", "Mdst_gamma"],
                   [HistStyle.error_bars(1), HistStyle.line_hist(4)])

    plot("p", 60, 0.103, 0.106, ";p_{#gamma};[MeV]")
    plot("costheta", 25, -1, 1, ";cos#theta_{#gamma};[]")
    plot("phi", 40, -R.TMath.Pi(), R.TMath.Pi(), ";#phi_{#gamma};[rad]")

if __name__ == "__main__":
    rootFile = "/gpfs/home/belle2/wangz/Work/SpinAlignment/steeringFile/belle1_steeringFile/hadronic_selection_test/hadronic_belle.root"
    main(rootFile)