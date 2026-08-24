import ROOT as R
from DRAW import style_draw, HistStyle
from OFFLINE_PROCESS import RDF_process
from math import pi


def test1():
    df = R.RDataFrame("event", "reco_HadronBorJ_e69r1344_1344.root")

    tools = RDF_process()

    df = df.Define("z_costheta", "1").Define("z_phi", "0")

    df = tools.cal_pola_angles(df, ("kp_px_cms", "kp_py_cms", "kp_pz_cms"),
                            ("km_px_cms", "km_py_cms", "km_pz_cms"),
                            (0.493677, 0.493677), ("z_costheta", "z_phi"), ("kp_index", "km_index"),
                            ("cos_theta_beam", "phi_beam"))

    df = df.Redefine("phi_beam", "fmod(phi_beam + 2 * TMath::Pi(), 2 * TMath::Pi())")
    df = df.Redefine("phi_beam", "TMath::Pi()/2 - abs(fmod(phi_beam, TMath::Pi()) - TMath::Pi()/2)")

    # ---

    phi_var = ["kp_index", "km_index"]
    for var in phi_var:
        df = df.Redefine(var, f"{var}[phi_thrust_costheta > 0]")
    df = df.Filter("kp_index.size()>0", "at least one pair with positive thrust costheta")

    df = tools.cal_pola_angles(df, ("kp_px_cms", "kp_py_cms", "kp_pz_cms"),
                            ("km_px_cms", "km_py_cms", "km_pz_cms"),
                            (0.493677, 0.493677), ("thrust[1]", "thrust[2]"), ("kp_index", "km_index"),
                            ("phi_pola_thrust_costheta", "phi_pola_thrust_phi"))

    df.Snapshot("event", "output.root")


def test2():
    df = R.RDataFrame("event", "output.root")

    hist_model = R.RDF.TH1DModel("", ";cos#theta;[]", 100, -1, 1)
    
    hist1 = df.Histo1D(hist_model, "phi_helicity_angle")
    hist2 = df.Histo1D(hist_model, "cos_theta_beam")

    style_draw([hist1, hist2], "comp.png", ["original", "calculated"], [HistStyle.filled_hist(2), HistStyle.line_hist(4)])

    df.Filter("cos_theta_beam.size() != phi_helicity_angle.size()")

    df = df.Filter("cos_theta_beam.size() == phi_helicity_angle.size()", "same cos_theta size")
    df = df.Define("Delta_cos_theta", "cos_theta_beam - phi_helicity_angle")
    hist_delta = df.Histo1D(("", ";#Delta cos#theta;[]", 100, -0.01, 0.01), "Delta_cos_theta")
    style_draw([hist_delta], "delta_cos_theta.png", styles =[HistStyle.line_hist(2)])

    
    # ---    
    df = df.Filter("phi_beam.size() == phi_helicity_phi.size()", "same phi size")
    hist_model_2 = R.RDF.TH1DModel("", ";#phi;[rad]", 100, 0, pi/2)
    hist_phi1 = df.Histo1D(hist_model_2, "phi_beam")
    hist_phi2 = df.Histo1D(hist_model_2, "phi_helicity_phi")
    style_draw([hist_phi1, hist_phi2], "comp_phi.png", ["original", "calculated"], [HistStyle.filled_hist(2), HistStyle.line_hist(4)])

    df = df.Define("Delta_phi", "phi_beam - phi_helicity_phi")
    hist_delta_phi = df.Histo1D(("", ";#Delta #phi;[rad]", 100, -0.01, 0.01), "Delta_phi")
    style_draw([hist_delta_phi], "delta_phi.png", styles =[HistStyle.line_hist(2)])


    df.Report().Print()

if __name__ == "__main__":
    test1()
    #test2()


