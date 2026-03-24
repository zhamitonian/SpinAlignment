#!/usr/bin/env python3

import ROOT as R
from DRAW import style_draw, HistStyle
from math import sqrt
import os

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

    style_draw([h_gen, h_reco], "", ["generate","reconstruct"],styles = [HistStyle.line_hist(4, 1, 2), HistStyle.error_bars(1)], pad = pad1, save=False)

    h_ratio = h_gen.Clone("h_ratio")
    h_ratio.Sumw2()
    h_ratio.Divide(h_gen, h_reco, 1, 1, "B")
    h_ratio.GetYaxis().SetTitle("ratio")
    
    # 调整 pad2 (ratio plot) 的坐标轴标签大小和字体
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

    png_file = os.path.join("./images/", f"{h_gen.GetName()}.png")
    if output_name is not None:
        png_file = os.path.join("./images/", f"{output_name}.png")
    style_draw([h_ratio], png_file, styles = [HistStyle.error_bars(1)], y_min= 0, y_max=2, use_user_y_range= True,pad= pad2)

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

def comparing():
    labels = ["all", "no_nu", "no_kl"]

    for label in labels:
        df = R.RDataFrame("event", f"./rootFiles/{label}.root")
        df = calculate_angle_between_axes(df, "thrust[1]", "thrust[2]", 
                                            "qqbar_axis[0]", "qqbar_axis[1]",
                                            "angle_thrust_qqbar")
        df = calculate_angle_between_axes(df, "thrust_truth[1]", "thrust_truth[2]", 
                                            "qqbar_axis[0]", "qqbar_axis[1]",
                                            "angle_thrustgen_qqbar")
        h_angle = df.Histo1D(("", ";angle;[rad]", 100, 0, 3.14/2), "angle_thrust_qqbar").GetValue()
        h_angle_gen = df.Histo1D(("", ";angle;[rad]", 100, 0, 3.14/2), "angle_thrustgen_qqbar").GetValue()
        style_draw([h_angle, h_angle_gen], f"./images/{label}_angle.png", leg_texts=["reconstruct" ,"generate"], styles=[HistStyle.line_hist(2), HistStyle.line_hist(4)], save = True)

        h_gen = df.Define("temp","thrust_truth[0]").Histo1D(("h_gen_thrust", "thrust distribution;thrust;[]", 100, 0.5, 1.01), "temp").GetValue()
        h_reco = df.Define("temp", "thrust[0]").Histo1D(("h_reco_thrust", "thrust distribution;thrust;[]", 100, 0.5, 1.01), "temp").GetValue()
        brush(h_gen, h_reco, output_name = f"{label}_thrust")    
    

if __name__ == "__main__":
    comparing()