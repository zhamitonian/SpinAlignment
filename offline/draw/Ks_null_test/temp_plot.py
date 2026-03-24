import ROOT as R
from DRAW import style_draw, HistStyle
import array


def load_rho00_to_th1(txt_path: str, hist_name: str = "hist_rho00") -> R.TH1F:
	z_vals = []
	rho_vals = []
	rho_errs = []

	with open(txt_path, "r") as f:
		for line in f:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			z, rho00, rho00_err = map(float, line.split())
			z_vals.append(z)
			rho_vals.append(rho00)
			rho_errs.append(rho00_err)

	#hist = R.TH1F(hist_name, ";z;#rho_{00}", len(z_vals), array.array("d", bin_edges))
	hist = R.TH1F(hist_name, ";z;#rho_{00}", 20, 0, 1)
	for i, (rho00, rho00_err) in enumerate(zip(rho_vals, rho_errs), start=1):
		bin_index = hist.GetXaxis().FindBin(z_vals[i-1])
		
		hist.SetBinContent(bin_index, rho00)
		hist.SetBinError(bin_index, rho00_err)

	hist.SetDirectory(0)
	return hist


txt_path_bin10 = "./images/rho00/data_bin10/rho00_results.txt"
txt_path_bin8 = "./images/rho00/data_bin8/rho00_results.txt"
hist_bin10 = load_rho00_to_th1(txt_path_bin10)
hist_bin8 = load_rho00_to_th1(txt_path_bin8)
c = style_draw([hist_bin10, hist_bin8], "./images/rho00_different_bin.png", ["fit with 10 bins", "fit with 8 bins"], [HistStyle.error_bars(2), HistStyle.error_bars(4)], save= False , y_min= 0, y_max= 0.6, use_user_y_range=True)
line = R.TLine(0., 1/3, 1, 1/3)
line.SetLineColor(R.kRed)
line.SetLineStyle(R.kDashed)
line.SetLineWidth(2)
line.Draw("same")

leg = R.TLegend(0.6, 0.3, 0.9, 0.55)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.AddEntry(line, "#rho_{00} = 1/3(unpolarized)", "l")

leg.Draw("same")
c.SaveAs("./images/rho00_different_bin.png")




