#!/usr/bin/env python3

import ROOT as R
from DRAW import style_draw, HistStyle
from math import sqrt
import os

def extract_reso(hist_reso, output_name = None):
    """
    extract resolution of variables by fitting to double gaussian 
    and plotting the fit result
    """
    title =  hist_reso.GetYaxis().GetTitle()
    title = title.replace("[","").replace("]","")
    """
    double_gaus = R.TF1("double_gaus", "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5])", hist_reso.GetXaxis().GetXmin(), hist_reso.GetXaxis().GetXmax())
    double_gaus.SetParLimits(0, 0, hist_reso.GetMaximum())
    double_gaus.SetParLimits(3, 0, hist_reso.GetMaximum())
    double_gaus.SetParLimits(2, 0.0001, 0.01)  
    double_gaus.SetParLimits(5, 0.0001, 0.01)    
    #double_gaus.SetParLimits(2, 0.0001, 0.01)
    #double_gaus.SetParLimits(5, 0.01, 0.1)    
    """
    double_gaus = R.TF1("double_gaus", "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [1], [4])", hist_reso.GetXaxis().GetXmin(), hist_reso.GetXaxis().GetXmax())
    double_gaus.SetParLimits(0, 0, hist_reso.GetMaximum())
    double_gaus.SetParLimits(3, 0, hist_reso.GetMaximum())
    double_gaus.SetParLimits(1, -0.02, 0.02)
    double_gaus.SetParLimits(2, 0.000001, 0.0001)  
    double_gaus.SetParLimits(4, 0.0001, 0.01)    
    double_gaus.SetParLimits(2, 0.01, 0.05)    
    #double_gaus.SetParLimits(2, 0.05, 0.4)  # pt
    #double_gaus.SetParLimits(4, 0.0001, 0.05)    

    hist_reso.Fit(double_gaus, "RQ")
    c = style_draw((hist_reso,), f"./images_reso/reso_{hist_reso.GetName()}.png", leg_texts=[], styles=[HistStyle.line_hist(1)], save = False)

    para = [double_gaus.GetParameter(i) for i in range(double_gaus.GetNpar())]
    #sigma = sqrt(para[0]/(para[0] +para[3]) *para[2]**2 + para[3]/(para[0] + para[3]) *para[5]**2 ) 
    sigma = sqrt(para[0]/(para[0] +para[3]) *para[2]**2 + para[3]/(para[0] + para[3]) *para[4]**2 ) 

    double_gaus.Draw("SAME")
    leg = R.TLegend(0.73, 0.7, 0.93, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetTextFont(22)
    leg.SetTextSize(0.03)
    leg.AddEntry(hist_reso, f"#delta {hist_reso.GetName()}", "l")
    leg.AddEntry(double_gaus, "Double Gaussian Fit", "l")

    if "MeV" in title:
        leg.AddEntry(0, f"#sigma = {sigma*1000:.2f} MeV", "")
    else:
        leg.AddEntry(0, f"#sigma = {sigma*1000:.2f} #times 10^{{-3}} {title}", "")
    leg.Draw("same")

    if output_name is not None:
        c.SaveAs(f"./images_reso/{output_name}.png")
    else:
        c.SaveAs(f"./images_reso/reso_{hist_reso.GetName()}.png")

    return sigma


def brush(h_gen, h_reco, output_name = None):
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

    style_draw([h_gen, h_reco], "", ["generate","reconstruct"],styles = [HistStyle.line_hist(4, 1, 2), HistStyle.error_bars(1)], pad = pad1, save=False, log_y= True)

    h_ratio = h_gen.Clone("h_ratio")
    h_ratio.Sumw2()
    h_ratio.Divide(h_gen, h_reco, 1, 1, "B")
    h_ratio.GetYaxis().SetTitle("ratio")
    
    h_ratio.GetXaxis().SetTitleOffset(0.8)  
    h_ratio.GetYaxis().SetTitleOffset(0.4)  

    png_file = os.path.join("./images_reso/", f"{h_gen.GetName()}.png")
    if output_name is not None:
        png_file = os.path.join("./images_reso/", f"{output_name}.png")
    style_draw([h_ratio], png_file, styles = [HistStyle.error_bars(1)], y_min= 0, y_max=2, use_user_y_range= True,pad= pad2)
    #style_draw([h_ratio], png_file, styles = [HistStyle.line_hist(1)], y_min= 0, y_max=2, use_user_y_range= True,pad= pad2)


def calculate_angle_between_axes(df: R.RDataFrame, 
                                  axis1_costheta: str, axis1_phi: str,
                                  axis2_costheta: str, axis2_phi: str,
                                  output_name: str = "angle_between_axes") -> R.RDataFrame:
    """
    Calculate angle between two axes given their spherical coordinates (costheta, phi)
    
    Args:
        df: Input RDataFrame
        axis1_costheta: costheta of first axis
        axis1_phi: phi of first axis
        axis2_costheta: costheta of second axis
        axis2_phi: phi of second axis
        output_name: name of output variable
    
    Returns:
        RDataFrame with angle between axes added
    """
    
    # Convert to Cartesian coordinates and calculate angle using dot product
    code = f"""
    double costheta1 = {axis1_costheta};
    double phi1 = {axis1_phi};
    double costheta2 = {axis2_costheta};
    double phi2 = {axis2_phi};
    
    double sintheta1 = sqrt(1 - costheta1*costheta1);
    double sintheta2 = sqrt(1 - costheta2*costheta2);
    
    // Cartesian coordinates
    double x1 = sintheta1 * cos(phi1);
    double y1 = sintheta1 * sin(phi1);
    double z1 = costheta1;
    
    double x2 = sintheta2 * cos(phi2);
    double y2 = sintheta2 * sin(phi2);
    double z2 = costheta2;
    
    // Dot product (both vectors are unit vectors)
    double cos_angle = x1*x2 + y1*y2 + z1*z2;
    cos_angle = std::abs(cos_angle);
    
    // Clamp to [0, 1] to avoid numerical issues
    cos_angle = std::min(1.0, cos_angle);
    
    return acos(cos_angle);
    """
    
    df = df.Define(output_name, code)
    
    return df


def phi_var(rootFile_path, cutThrust = False):
    var_config = {"M": (-0.01, 0.01, 100, "MeV", 0.99, 1.06, 200),
              "z": (-0.01, 0.01, 100, "", 0.2, 1, 100),
              "xp": (-0.01, 0.01, 100, "", 0, 1, 100),
              "helicity_angle": (-0.05, 0.05, 100, "", -1, 1, 100),
              "thrust_pt": (-1.5, 1.5, 300, "MeV", 0, 3.5, 100),
              "pt_zQ": (-0.25, 0.25, 100, "", 0, 0.5, 100)}
    
    var_latex_dict = {"M":"M_{#phi}", "z":"z", "xp":"x_{p}", "helicity_angle":"cos#theta^{#star}", "thrust_pt":"p_{t}", "pt_zQ":"p_{t}/zQ"}
    
    df = R.RDataFrame("event", rootFile_path)
    sqrts = 10.516469955444336 
    df = df.Define("phi_pt_zQ", f"phi_thrust_pt / phi_z/{sqrts}").Define("phi_gen_pt_zQ", f"phi_gen_thrust_pt / phi_gen_z/{sqrts}")
    if cutThrust:
        df = df.Filter("thrust > 0.8", "thrust cut")
    #for var in ["M", "z", "helicity_angle", "thrust_pt", "xp"]:
    for var in ["pt_zQ"]:
        varname = f"phi_{var}"

        df = df.Filter(f"{varname}.size() == phi_gen_{var}.size()", "match count")
        df = df.Filter(f"{varname}.size() > 0", "has phi after match")

        reso_hist = df.Define("diff", f"{varname} - phi_gen_{var}").Histo1D((f"{var_latex_dict[var]}", f";#Delta {var_latex_dict[var]};[{var_config[var][3]}]", 
                                                    var_config[var][2], var_config[var][0], var_config[var][1]), "diff").GetValue()

        suffix = "_cutted" if cutThrust else ""

        sigma = extract_reso(reso_hist, f"reso_{var}{suffix}")

        h_gen = df.Histo1D((f"h_gen_{var}", f"{var} distribution;{var_latex_dict[var]};[{var_config[var][3]}]", 
                                    var_config[var][6], var_config[var][4], var_config[var][5]), f"phi_gen_{var}").GetValue()
        h_reco = df.Histo1D((f"h_reco_{var}", f"{var} distribution;{var_latex_dict[var]};[{var_config[var][3]}]", 
                                    var_config[var][6], var_config[var][4], var_config[var][5]), varname).GetValue()

        brush(h_gen, h_reco, f"{var}{suffix}")

        print(f"{var} resolution: {sigma*1000:.2f} {var_config[var][3]}")

        df.Snapshot("event", "temp.root")

def event_var(rootFile_path, cutThrust = False):

    df = R.RDataFrame("event", rootFile_path.replace("_reco_truth_matched", ""))
    if cutThrust:
        df = df.Filter("thrust[0] > 0.8", "thrust cut")
    df = calculate_angle_between_axes(df, "thrust[1]", "thrust[2]", 
                                        "qqbar_axis[0]", "qqbar_axis[1]",
                                        "angle_thrust_qqbar")
    df = calculate_angle_between_axes(df, "thrust_truth[1]", "thrust_truth[2]", 
                                        "qqbar_axis[0]", "qqbar_axis[1]",
                                        "angle_thrustgen_qqbar")

    h_angle = df.Histo1D(("", ";angle;[rad]", 100, 0, 3.14/2), "angle_thrust_qqbar").GetValue()
    h_angle_gen = df.Histo1D(("", ";angle;[rad]", 100, 0, 3.14/2), "angle_thrustgen_qqbar").GetValue()
    suffix = "_cutted" if cutThrust else ""
    outfile = f"./images_reso/angle_thrust_qqbar{suffix}.png"
    style_draw([h_angle, h_angle_gen], outfile, leg_texts=["reconstruct" ,"generate"], styles=[HistStyle.line_hist(2), HistStyle.line_hist(4)], save = True)

    reso_hist = df.Define("diff", "thrust[0] - thrust_truth[0]").Histo1D(("thrust", ";#Delta thrust;[]", 100, -0.5, 0.5), "diff").GetValue()
    sigma = extract_reso(reso_hist, f"reso_thrust{suffix}")
    h_gen = df.Define("temp","thrust_truth[0]").Histo1D(("h_gen_thrust", "thrust distribution;thrust;[]", 100, 0.5, 1.01), "temp").GetValue()
    h_reco = df.Define("temp", "thrust[0]").Histo1D(("h_reco_thrust", "thrust distribution;thrust;[]", 100, 0.5, 1.01), "temp").GetValue()
    brush(h_gen, h_reco, output_name = None if not cutThrust else "thrust_cutted")    

    print(f"thrust resolution: {sigma:.4f} ")


if __name__ == "__main__":
    rootFile_path = "/gpfs/group/belle2/users2022/luruihua/for_wangz/data_gMC_belle1/2025-11-25_SpinAlignment_gMC/continuum_reco_truth_matched.root" 
    rootFile_path_evt = "/gpfs/group/belle2/users2022/luruihua/for_wangz/data_gMC_belle1/2025-11-25_SpinAlignment_gMC/continuum.root" 
    #event_var(rootFile_path_evt)
    #event_var(rootFile_path_evt, cutThrust= True)
    phi_var(rootFile_path, cutThrust= False)
    #phi_var(rootFile_path, cutThrust= True)
    

