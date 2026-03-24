import ROOT as R
from PHY_CALCULATOR import PhysicsCalculator
from DRAW import style_draw , HistStyle
import os

def efficiency_calculation(input_rootFile ,output_dir, nsig_txt):
    
    tools = PhysicsCalculator(os.path.join(output_dir,"container.root"))
    tools.set_bins([0,1], bin_num=[40])

    h_nsig =  tools.getNsigHist(nsig_txt)
    style_draw([h_nsig], os.path.join(output_dir, "nsig_xp.png"), styles=[HistStyle.error_bars(1)])

    df_truth = R.RDataFrame("truth", input_rootFile)
    h_truth = df_truth.Histo1D(("h_truth", ";phi_xp;Counts", 40, 0, 1), "phi_xp").GetValue()
    style_draw([h_truth], os.path.join(output_dir, "truth_xp.png"), styles=[HistStyle.error_bars(1)])

    h_nsig.Divide(h_truth)
    style_draw([h_nsig], os.path.join(output_dir, "efficiency.png"), styles=[HistStyle.error_bars(1)], y_min= 0, y_max = 0.75, use_user_y_range=True)

    """
    eff = R.TEfficiency(h_nsig, h_truth)
    eff.SetName("efficiency_phi_xp")
    eff.SetTitle("Efficiency vs phi_xp;phi_xp;Efficiency")
    c_eff = R.TCanvas("c_eff", "c_eff", 800, 600)
    R.gStyle.SetOptStat(0)
    eff.Draw()
    c_eff.SaveAs("./efficiency_phi_xp.png")
    """

if __name__ == "__main__":
    #input_rootFile = "/gpfs/group/belle2/users2022/luruihua/for_wangz/data_gMC_belle1/2025-11-13_SpinAlignment_gMC/continuum_truth_processed.root"
    #output_dir = "./images/"
    #nsig_txt = "./images/fitting/nsig_results.txt"
    
    #input_rootFile = "/gpfs/group/belle2/users2022/luruihua/for_wangz/data_gMC_belle1/2025-11-25_SpinAlignment_gMC/continuum_truth_processed.root"
    #output_dir = "./eff_Nov26/"
    #nsig_txt = "./eff_Nov26/fitting/nsig_results.txt"
    
    #input_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/2025-12-02_SpinAlignment_qqbarMC/continuum_truth_processed_test.root"
    #utput_dir = "./eff_Dec02/"
    #nsig_txt = "./eff_Dec02/fitting/nsig_results.txt"
    #efficiency_calculation(input_rootFile, output_dir, nsig_txt)

    #input_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/2025-12-02_SpinAlignment_qqbarMC2/continuum_truth_processed.root"
    #output_dir = "./eff_Dec02_2/"
    #nsig_txt = "./eff_Dec02_2/fitting/nsig_results.txt"
    #efficiency_calculation(input_rootFile, output_dir, nsig_txt)

    input_rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/2025-12-03_SpinAlignment_qqbarMC/continuum_truth_processed_cutPt.root"
    output_dir = "./eff_Dec03/"
    nsig_txt = "./eff_Dec03/fitting/nsig_results.txt"
    efficiency_calculation(input_rootFile, output_dir, nsig_txt)