import ROOT as R
from OFFLINE_PROCESS import RDF_process
import sys
sys.path.append("/home/belle2/wangz/Work/SpinAlignment/offline/draw/belle_kid_eff/")
from get_pid_eff import get_pid_eff
from DRAW import style_draw, HistStyle

def Ks_reco(df:R.RDataFrame):
    df = df.Define("pip_pt", "pip_p* sqrt(1 - pip_costheta*pip_costheta)")
    df = df.Define("pim_pt", "pim_p* sqrt(1 - pim_costheta*pim_costheta)")

    #track_selection = "{0}_p >0.5 && {0}_costheta > -0.511 && {0}_costheta < 0.842 && {0}_atcKPi < 0.4"
    track_selection = "{0}_costheta > -0.511 && {0}_costheta < 0.842 && {0}_atcKPi < 0.4 && {0}_pt > 0.1"
    Ks_daughter_tracks_cut = track_selection.format("pip") + " && " + track_selection.format("pim")
    df = df.Define("Ks_daughter_tracks_cut", Ks_daughter_tracks_cut)

    for var in ["E_cms", "px_cms", "pz_cms", "p", "costheta", "atcKPi"]:
        for particle in ["pip", "pim"]:
            df = df.Redefine(f"{particle}_{var}", f"{particle}_{var}[Ks_daughter_tracks_cut]")

    tools = RDF_process()
    df = tools.save_all_pairs(df,
                              p1_branches=("pip_px_cms","pip_py_cms","pip_pz_cms"),
                              p2_branches=("pim_px_cms","pim_py_cms","pim_pz_cms"),
                              mass=(0.13957061, 0.13957061),
                              particle_name=("Ks","pip","pim"),
                              cross_mode= False) 

    df = df.Define("Ks_M", "sqrt(Ks_E*Ks_E - Ks_px*Ks_px - Ks_py*Ks_py - Ks_pz*Ks_pz)")
    df = df.Define("Ks_z", "2*Ks_E/sqrts")

    df = get_pid_eff(df, period="svd2", is_pion=True, cut_value=6, particle_names = ("pip", "pim"), eff_source="data")

    df = tools.convert_cartesian_to_spherical(df, particles=["Ks"])
    df = tools.convert_cartesian_to_spherical(df, particles=["pip", "pim"],
                                               px_branch = "px_cms", py_branch = "py_cms", pz_branch = "pz_cms")

    df = df.Define("Ks_eff", "pip_pid_eff * pim_pid_eff")

    df = df.Filter("Ks_M.size() > 0", "has Ks")
    df.Report().Print()

    df = df.Define("Ks_eff_correction", "1.0/Ks_eff")
    hist = df.Histo1D(("", ";M(K_{s});[MeV]",100 , 0.47, 0.52), "Ks_M")
    hist2 = df.Histo1D(("", ";M(K_{s});[MeV]",100 , 0.47, 0.52), "Ks_M", "Ks_eff_correction")

    style_draw([hist, hist2], "Ks_M.png", ["", "with eff correction"], [HistStyle.error_bars(2), HistStyle.error_bars(4)])

    Br2Save = ["Ks_M", "Ks_z","Ks_helicity_angle", "Ks_eff", "Ks_eff_correction"]
    Br2Save += ["Ks_costheta", "Ks_phi", "Ks_p", "pip_p", "pip_costheta", "pip_phi", "pim_p", "pim_costheta", "pim_phi", "pip_pid_eff", "pim_pid_eff"]
    Br2Save += ["pip_index", "pim_index"]

    df.Snapshot("event", "./Ks_processed_with_pt_cut.root", Br2Save)

if __name__ == "__main__":
    rootFile = "/gpfs/group/belle/users/dues/data_gMC_belle1/KsSpinAlignment_only_goodKs_pid_v1.0.0/steered_svd2.root"
    df = R.RDataFrame("event", rootFile)
    Ks_reco(df)
