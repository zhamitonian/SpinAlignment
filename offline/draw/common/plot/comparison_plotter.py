"""Plot helpers for data vs MC overlays."""

import os
import ROOT as R
from DRAW import style_draw, HistStyle


def plot_data_mc(h_data, h_mc, output_path, leg_position=2, normalize=True):
    """Overlay data and MC histograms and save."""
    if normalize and h_mc.Integral() != 0:
        scale = h_data.Integral() / h_mc.Integral()
        h_mc.Scale(scale)

    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    style_draw(
        [h_data, h_mc],
        output_path,
        ["sPlotted Data", "Truth Matched MC"],
        [HistStyle.error_bars(R.kBlack), HistStyle.line_hist(4)],
        legend_position=leg_position,
    )

    """ version with pull distribution
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

        style_draw([h_data ,h_mc], os.path.join(self.output_dir, f"{var}.png"),
                ["sPlotted Data", "Truth Matched MC"],
                [HistStyle.error_bars(R.kBlack), HistStyle.line_hist(4)], 
                legend_position = leg_position,pad = pad1, save = False)
        
        h_pull = h_mc.Clone("h_pull")
        for i in range(1, h_pull.GetNbinsX() + 1):
            data_val = h_data.GetBinContent(i)
            data_err = h_data.GetBinError(i)
            mc_val = h_mc.GetBinContent(i)
            mc_err = h_mc.GetBinError(i)
            
            total_err = (data_err**2 + mc_err**2)**0.5
            if total_err > 0:
                pull = (data_val - mc_val) / total_err
                h_pull.SetBinContent(i, pull)
                h_pull.SetBinError(i, 1.0)  # Pull 的误差为 1
            else:
                h_pull.SetBinContent(i, 0)
                h_pull.SetBinError(i, 0)

    
        h_pull.GetYaxis().SetTitle("Pull")
        h_pull.GetXaxis().SetTitleOffset(0.75)
        h_pull.GetYaxis().SetTitleOffset(0.35)
        style_draw([h_pull], os.path.join(self.output_dir, f"{var}.png"), styles=[HistStyle.filled_line_hist(R.kGray+1, 1001)], 
                   y_min=-4, y_max=4, use_user_y_range=True, pad=pad2, save=True)
    """