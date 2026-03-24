import ROOT as R
from DRAW import style_draw , HistStyle

def phi_multi(rootFile:str, treeName:str="event"):
    df = R.RDataFrame(treeName, rootFile)
    hist = df.Define("phi_multplicity", "phi_M.size()").Histo1D(("", ";#phi multiplicity;[]", 6, 1, 7), "phi_multplicity").GetValue()
    #style_draw([hist], "./phi_multiplicity.png", styles=[HistStyle.error_bars(1)], log_y= True)
    return hist

if __name__ == "__main__":
    gMC_root = "../Bining_var_resolution/reso_check_gMC.root"
    truth_root = "../Costheta_z_efficiency/truth_processed.root"
    h_gMC = phi_multi(gMC_root)#.Scale(0.25)
    h_truth = phi_multi(truth_root, treeName="truth")#.Scale(0.25)

    style_draw([h_gMC, h_truth], "./Phi_multiplicity_comparison.png", 
               styles=[HistStyle.error_bars(2), HistStyle.error_bars(4)], 
               leg_texts=["udsc MC", "Generated"], log_y = True)

    gMC_1 = h_gMC.GetBinContent(1)
    truth_1 = h_truth.GetBinContent(1)
    gMC_2 = h_gMC.GetBinContent(2)
    truth_2 = h_truth.GetBinContent(2)
    print(f"gMC 1 phi: {gMC_1}, gMC 2 phi: {gMC_2}, ratio: {gMC_2/gMC_1:.4f}")
    print(f"truth 1 phi: {truth_1}, truth 2 phi: {truth_2}, ratio: {truth_2/truth_1:.4f}")

"""
gMC 1 phi: 1752245.0, gMC 2 phi: 23621.0, ratio: 0.0135
truth 1 phi: 17720778.0, truth 2 phi: 522828.0, ratio: 0.0295
"""