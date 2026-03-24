import ROOT as R
from DRAW import HistStyle, style_draw
from OFFLINE_PROCESS import RDF_process

b2bii_path = "/gpfs/group/belle2/users2022/wangz/other_rootFiles/SpinAlignment/b2bii_basf_comparing_rootProd/b2bii.root"
basf_path = "/gpfs/group/belle2/users2022/wangz/other_rootFiles/SpinAlignment/b2bii_basf_comparing_rootProd/basf/exp71_rs2249_re2348_evtgen-uds_0_tree.root"
basf_path_rotated = "/gpfs/group/belle2/users2022/wangz/other_rootFiles/SpinAlignment/b2bii_basf_comparing_rootProd/basf/after_rotate.root"  ## after coordinates rotation

def check_cosTheta_thrust():
    df_basf = R.RDataFrame("event", basf_path)
    df_b2bii = R.RDataFrame("event_trk", b2bii_path)
    count_b2bii = int(df_b2bii.Count().GetValue())
    count_basf = int(df_basf.Count().GetValue())
    print(f"event entries,\n b2bii: {count_b2bii}\nbasf: {count_basf}")

    df_basf = df_basf.Define("temp", "thrust[1]").Define("abs_cosTheta_thrust", "abs(temp)")
    df_b2bii = df_b2bii.Define("abs_cosTheta_thrust", "abs(thrustAxisCosTheta)")

    h_basf = df_basf.Histo1D(("h_basf", "cos#theta_{thrust} from BASF;|cos#theta_{thrust}|;[1]", 50, 0, 1), "abs_cosTheta_thrust").GetValue()
    h_b2bii = df_b2bii.Histo1D(("h_b2bii", "cos#theta_{thrust} from B2BII;|cos#theta_{thrust}|;[1]", 50, 0, 1), "abs_cosTheta_thrust").GetValue()

    # plot
    c_all = R.TCanvas("c_all", "Combined", 1600, 1080)

    pad1 = R.TPad("pad1", "pad1", 0, 0.3, 1, 1)  
    pad1.SetBottomMargin(0.01)  
    pad1.SetTopMargin(0.1)
    pad1.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad1.Draw()
    
    pad2 = R.TPad("pad2", "pad2", 0, 0, 1, 0.3)  
    pad2.SetTopMargin(0.02)     
    pad2.SetBottomMargin(0.3)
    pad2.SetLeftMargin(0.15)
    pad2.SetRightMargin(0.05)
    
    pad2.Draw()

    c_all.Update()

    style_draw([h_basf, h_b2bii], "", ["basf","b2bii"], [HistStyle.error_bars(1), HistStyle.line_hist(4)], legend_position= 0, pad =pad1, save= False)

    h_ratio = h_b2bii.Clone("h_ratio")
    h_ratio.Sumw2()
    h_ratio.Divide(h_basf, h_b2bii, 1, 1, "B")
    h_ratio.GetYaxis().SetTitle("ratio")

    h_ratio.GetXaxis().SetLabelSize(0.10)   # X轴标签大小，因为pad2较小，需要更大
    h_ratio.GetYaxis().SetLabelSize(0.10)   # Y轴标签大小
    h_ratio.GetXaxis().SetTitleSize(0.16)   # X轴标题大小
    h_ratio.GetYaxis().SetTitleSize(0.14)   # Y轴标题大小
    h_ratio.GetXaxis().SetLabelFont(22)     # 字体 (42=Helvetica)
    h_ratio.GetYaxis().SetLabelFont(22)
    h_ratio.GetXaxis().SetTitleFont(22)
    h_ratio.GetYaxis().SetTitleFont(22)
    h_ratio.GetXaxis().SetTitleOffset(0.8)  # 标题偏移
    h_ratio.GetYaxis().SetTitleOffset(0.4)  # Y轴标题偏移

    style_draw([h_ratio], "abs_cosTheta_thrust.png", styles= [HistStyle.error_bars(1)] ,pad =pad2, y_min= 0.9, y_max= 1.1, use_user_y_range= True)

def calculate_trk_pt():
    """
    """
    df_basf = R.RDataFrame("event", basf_path)
    df_b2bii = R.RDataFrame("event_trk", b2bii_path)

    tools = RDF_process()

    df_basf.Define("thrust_theta","TMath::ACos(thrust[1])")
    df_basf = tools.calculate_pt_toAxis(df_basf, ("trk_px_CMS", "trk_py_CMS", "trk_pz_CMS"),
                                        ("TMath::ACos(thrust[1])", "thrust[2]"), "trk", "thrust")
    hist = df_basf.Histo1D(("", ";pt;[MeV]",40, 0, 2),"trk_pt_toAxis_thrust")
    style_draw([hist], "test.png", styles= [HistStyle.line_hist(4)] )
        
    hist_b2bii = df_b2bii.Histo1D(("thrustFrame_pt", ";pt;[MeV]",40, 0, 2), "thrustFrame_pt")

    style_draw([hist_b2bii, hist], "test_b2bii.png", styles= [HistStyle.line_hist(4), HistStyle.error_bars(1)] )





if __name__ == "__main__":
    #check_cosTheta_thrust()
    calculate_trk_pt()