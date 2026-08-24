import ROOT as R
import os
import sys
sys.path.append("/home/belle2/wangz/Work/SpinAlignment/offline")
from KsSA_processor import Ks_reco, Ks_reco_truth_matched
from FIT.generic_fit import *
from DRAW import style_draw, HistStyle

def process_rootFile(rootFile, sig_path, bkg_path):
    if not os.path.exists(sig_path) or not os.path.exists(bkg_path):
        rdf = R.RDataFrame("event", rootFile)
        rdf_sig = rdf.Filter("n_ks_truth > 0", "signal region")
        rdf_sig, Br2Save = Ks_reco(rdf_sig)
        rdf_sig.Snapshot("event", sig_path, Br2Save)
        rdf_bkg = rdf.Filter("n_ks_truth == 0", "background region")
        rdf_bkg, Br2Save = Ks_reco(rdf_bkg)
        rdf_bkg.Snapshot("event", bkg_path)
    else :
        rdf_sig = R.RDataFrame("event", sig_path)
        rdf_bkg = R.RDataFrame("event", bkg_path)

    hist_model = R.RDF.TH1DModel("", ";M(K_{s});[MeV]", 200, 0.47, 0.52)
    h_mKs_sig = rdf_sig.Histo1D(hist_model, "Ks_M")
    h_mKs_bkg = rdf_bkg.Histo1D(hist_model, "Ks_M")
    style_draw([h_mKs_sig, h_mKs_bkg], "./images/mKs_dist.png", ["signal", "background"], 
               [HistStyle.filled_line_hist(4,3011), HistStyle.line_hist(2)], log_y=True)
    style_draw([h_mKs_bkg], "./images/mKs_bkg.png", ["background"],
               [HistStyle.filled_line_hist(4,3011)])

    hist_model = R.RDF.TH2DModel("", ";cos#theta*;z", 10, -1, 1, 15, 0.05, 0.8)
    hist_cosTheta_z = rdf_bkg.Histo2D(hist_model , "Ks_helicity_angle", "Ks_z", "Ks_M")
    hist_cosTheta_z.SetMarkerSize(1.5) 
    c = R.TCanvas("c_eff", "c_eff", 1600, 1080)
    R.gStyle.SetPaintTextFormat(".1f")
    hist_cosTheta_z.Draw("COLZ TEXT")
    c.SaveAs(os.path.join("images", "bkg_Ks_costheta_vs_z.png"))

def check_isSignal(rootFile):  # get truth-matched MC
    rootFile_bkg = "./rootFiles/bkg_isSignal_v2.3.3_with_MCGen.root"
    rootFile_sig = "./rootFiles/sig_isSignal_v2.3.3.root"
    if not os.path.exists(rootFile_bkg) or not os.path.exists(rootFile_sig):
        rdf = R.RDataFrame("event", rootFile)
        #rdf = rdf.Filter("n_ks_truth == 0", "background region")

        from OFFLINE_PROCESS import RDF_process
        sys.path.append("/home/belle2/wangz/Work/SpinAlignment/offline/draw/belle_kid_eff/")
        from get_pid_eff_weight import get_weight
        tools = RDF_process()
        rdf = tools.save_all_pairs(rdf,
                                p1_branches=["pip_px_cms","pip_py_cms","pip_pz_cms"],
                                p2_branches=["pim_px_cms","pim_py_cms","pim_pz_cms"],
                                mass=(0.13957061, 0.13957061),
                                particle_name=("Ks","pip","pim"),
                                cross_mode= False)
        rdf = get_weight(rdf, period = "svd2", is_pion=True, cut_value=6, particle_names = ("pip", "pim"))
        rdf = rdf.Define("Ks_weight", "pip_pid_weight * pim_pid_weight")
        rdf = tools.convert_cartesian_to_spherical(rdf, particles=["Ks"])
        rdf = tools.convert_cartesian_to_spherical(rdf, particles=["pip", "pim"],
                                                px_branch = "px_cms", py_branch = "py_cms", pz_branch = "pz_cms")
        
        rdf = rdf.Define("Ks_M", "sqrt(Ks_E*Ks_E - Ks_px*Ks_px - Ks_py*Ks_py - Ks_pz*Ks_pz)")
        rdf = rdf.Define("Ks_z", "2*Ks_E/sqrts")

        Ks_variables = ["M", "z", "p", "costheta", "phi", "helicity_angle", "weight"]
        pion_variables = ["p", "costheta", "phi", "pid_weight", "index"]
        rdfs = {}
        for label, condition, output in [
            ("bkg", "!pip_isSignal || !pim_isSignal", rootFile_bkg),
            ("sig", "pip_isSignal && pim_isSignal", rootFile_sig),
        ]:
            rdf_out = rdf
            for var in Ks_variables:
                rdf_out = rdf_out.Redefine(f"Ks_{var}", f"Ks_{var}[{condition}]")
            for var in pion_variables:
                rdf_out = rdf_out.Redefine(f"pip_{var}", f"pip_{var}[{condition}]")
                rdf_out = rdf_out.Redefine(f"pim_{var}", f"pim_{var}[{condition}]")
            rdf_out = rdf_out.Filter("Ks_M.size() > 0", "Ks list not empty")
            rdf_out.Report().Print()
            mcgen_branches = [str(c) for c in rdf_out.GetColumnNames() if str(c).startswith("MCGen")] + ["nMCGen"]
            output_branches = [f"pip_{var}" for var in pion_variables] + [f"pim_{var}" for var in pion_variables] + [f"Ks_{var}" for var in Ks_variables]
            if label == "bkg":
                output_branches += mcgen_branches
                
            rdf_out.Snapshot("event", output, output_branches)
            rdfs[label] = rdf_out

        rdf_bkg = rdfs["bkg"]
 
    else:
        rdf_bkg = R.RDataFrame("event", rootFile_bkg)

    hist_model = R.RDF.TH1DModel("", ";M(K_{s});[MeV]", 200, 0.47, 0.52)
    h_bkg_isSignal = rdf_bkg.Histo1D(hist_model, "Ks_M")

    rdf_bkg_truth = R.RDataFrame("event", "rootFiles/obsolete/svd2_4Soffres_st0_v2.3.3_bkg.root")
    h_bkg_truth = rdf_bkg_truth.Histo1D(hist_model, "Ks_M")
    style_draw([h_bkg_isSignal, h_bkg_truth], "./images/bkg_isSignal_check.png", ["bkg with isSignal", "bkg with truth match"], 
                [HistStyle.line_hist(2), HistStyle.line_hist(4)], log_y=True)  

def bkg_fit(rootFile, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    file = R.TFile(rootFile)
    tree = file.Get("event")
    branches_name = ["Ks_M"]
    log_file = os.path.join(output_dir, "fit_log.log")

    Ks_z_conf = (0.05, 0.80, 15)
    Ks_helicity_angle_conf = (-1, 1, 10)
    Ks_z_width = (Ks_z_conf[1] - Ks_z_conf[0]) / Ks_z_conf[2]
    Ks_helicity_angle_width = (Ks_helicity_angle_conf[1] - Ks_helicity_angle_conf[0]) / Ks_helicity_angle_conf[2]
    pad_width = len(str(Ks_helicity_angle_conf[2] * Ks_z_conf[2]))

    result, nbkg, nbkg_err = None, None, None

    #for i in range(Ks_z_conf[2]):
    for i in range(2):
        for j in range(Ks_helicity_angle_conf[2]):
            z_low = Ks_z_conf[0] + i * Ks_z_width
            z_high = Ks_z_conf[0] + (i + 1) * Ks_z_width
            h_low = Ks_helicity_angle_conf[0] + j * Ks_helicity_angle_width
            h_high = Ks_helicity_angle_conf[0] + (j + 1) * Ks_helicity_angle_width

            cut_expr = f"Ks_z > {z_low} && Ks_z < {z_high}"
            cut_expr += f" && Ks_helicity_angle > {h_low} && Ks_helicity_angle < {h_high}"

            index = i * Ks_helicity_angle_conf[2] + j

            # FIX: avoid "Input and output lists are the same!"
            if index > 0:
                tree.SetEventList(0)  # clear current input list
                R.gDirectory.Delete(f"{list_name};*")  # remove old object with same name if any
            list_name = f"event_list_{i}_{j}"

            tree.Draw(f">>{list_name}", cut_expr, "goff")
            event_list = R.gDirectory.Get(list_name)
            if not event_list:
                continue

            tree.SetEventList(event_list)

            for order in range(3):
                pdf_config = FitDefinition(
                    [Variable("Ks_M", 0.47, 0.52, 100)],
                    #[Variable("Ks_M", 0.48, 0.51, 120)],
                    [PDFSpec("bkg", "Ks_M", "chebychev", {"order": order, "coef1": (0, -10, 10), "coef2": (0, -10, 10)})],
                    model="bkg",
                )
                plot_config = PlotConfiguration(plot_config={
                    "xlabel": {"Ks_M": "M_{#pi^{+}#pi^{-}} (GeV/c^{2})"},
                    "components": {"model": {"label": "Total Fit", "color": 4}},
                    "legend": {"extra_text": [f"bkg order {order}", f"z: [{z_low:.3f}, {z_high:.3f}]", f"helicity angle: [{h_low:.3f}, {h_high:.3f}]"]},
                })
                dataset_config = DatasetConfig(
                    binned_fit=True, branches_name=branches_name, perform_splot=False, weight_branch="Ks_weight"
                )
                fit_config = FitterConfig(two_step_fit=True, use_minos=False)

                fitter = GenericFit(
                    tree, output_dir, log_file=log_file, fit_definition=pdf_config,
                    dataset_config=dataset_config, plot_config=plot_config, fitter_config=fit_config
                )
                result, fit_results = fitter.run()
                #nbkg, nbkg_err = fit_results["nsig"], fit_results["nsig_err"]

                src = os.path.join(output_dir, "_Ks_M.png")
                #dst = os.path.join(output_dir, f"bin{index:0{pad_width}d}_order{order}_Ks_M_narrower.png")
                dst = os.path.join(output_dir, f"bin{index:0{pad_width}d}_order{order}_Ks_M.png")
                if os.path.exists(src):
                    os.rename(src, dst)

    file.Close()
    #return result, nbkg, nbkg_err
    return result, 0, 0

def sig_fit(rootFile, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    file = R.TFile(rootFile)
    tree = file.Get("event")
    branches_name = ["Ks_M"]

    Ks_z_conf = (0.05, 0.80, 15)
    Ks_helicity_angle_conf = (-1, 1, 10)
    Ks_z_width = (Ks_z_conf[1] - Ks_z_conf[0]) / Ks_z_conf[2]
    Ks_helicity_angle_width = (Ks_helicity_angle_conf[1] - Ks_helicity_angle_conf[0]) / Ks_helicity_angle_conf[2]
    pad_width = len(str(Ks_helicity_angle_conf[2] * Ks_z_conf[2]))

    result, nbkg, nbkg_err = None, None, None

    #for i in range(Ks_z_conf[2]):
    #    for j in range(Ks_helicity_angle_conf[2]):
    for i in range(1):
        for j in range(10):
            z_low = Ks_z_conf[0] + i * Ks_z_width
            z_high = Ks_z_conf[0] + (i + 1) * Ks_z_width
            h_low = Ks_helicity_angle_conf[0] + j * Ks_helicity_angle_width
            h_high = Ks_helicity_angle_conf[0] + (j + 1) * Ks_helicity_angle_width

            cut_expr = f"Ks_z > {z_low} && Ks_z < {z_high}"
            cut_expr += f" && Ks_helicity_angle > {h_low} && Ks_helicity_angle < {h_high}"

            index = i * Ks_helicity_angle_conf[2] + j
            log_file = os.path.join(output_dir, f"bin{index:0{pad_width}d}.log")

            # FIX: avoid "Input and output lists are the same!"
            if index > 0:
                tree.SetEventList(0)  # clear current input list
                R.gDirectory.Delete(f"{list_name};*")  # remove old object with same name if any
            list_name = f"event_list_{i}_{j}"

            tree.Draw(f">>{list_name}", cut_expr, "goff")
            event_list = R.gDirectory.Get(list_name)
            if not event_list:
                continue

            tree.SetEventList(event_list)

            pdf_config = FitDefinition([Variable("Ks_M", 0.47, 0.52, 200)],
                            [PDFSpec("sig", "Ks_M", "crystal_ball",
                                     {"mean":(0.497,0.496,0.498),
                                      "sigma":(0.001,0.0005,0.01),
                                      "alpha":(1.5,0.01,5),
                                      "n":(2.0,0.01,100),
                                      "n_right":(2.0,0.01,100),
                                      "alpha_right":(1.5,0.01,5),}),
                                      #"sigma_right":(0.001,0.0005,0.01)}),
                            PDFSpec("gauss", "Ks_M", "gaussian",
                                    {"mean":(0.497, 0.496, 0.498),
                                     "sigma":(0.0006, 0.0005, 0.005)}),],
                            model = "nsig[7000,0,3000000]* (frac[0.8,0.75,1]*sig + gauss)")
            plot_config = PlotConfiguration(plot_config={
                "xlabel" : {"Ks_M" : "M_{#pi^{+}#pi^{-}} (GeV/c^{2})"},
                "components" : {"model" : {"label" : "Total Fit", "color" : 4},
                                "sig": {"label" : "Signal", "color" : 2, "style" : 4, "width" : 3}, # plot will draw MCshape conv gauss
                                "gauss": {"label" : "Gaussian", "color" : R.kGreen + 2, "style" : 7, "width" : 3},},
                "legend": {"extra_text": [f"z: [{z_low:.3f}, {z_high:.3f}]", f"helicity angle: [{h_low:.3f}, {h_high:.3f}]"]}},)
            dataset_config = DatasetConfig(binned_fit= True, branches_name=branches_name, perform_splot=False,weight_branch="Ks_weight")
            fit_config = FitterConfig(two_step_fit=True, use_minos=False)

            fitter = GenericFit(tree, output_dir, log_file=log_file,fit_definition=pdf_config,
                            dataset_config=dataset_config, plot_config=plot_config, fitter_config=fit_config)
            result, fit_results = fitter.run() 
            nbkg, nbkg_err = fit_results["nsig"], fit_results["nsig_err"]

            src = os.path.join(output_dir, "_Ks_M.png")
            dst = os.path.join(output_dir, f"bin{index:0{pad_width}d}_Ks_M.png")
            if os.path.exists(src):
                os.rename(src, dst)

    return result, nbkg, nbkg_err

def isSignal_before_after():
    rootFile1 = "./rootFiles/bkg_isSignal_check.root"
    rootFiles = "./rootFiles/bkg_isSignal_v2.3.3_with_MCGen.root"
    rootFile3 = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.3_qqbar_svd2/svd2_st0_reco_processed.root"
    rdf1 = R.RDataFrame("event", rootFile1)
    rdf2 = R.RDataFrame("event", rootFiles)
    rdf3 = R.RDataFrame("event", rootFile3)
    hist_model = R.RDF.TH1DModel("", ";M(K_{s});[MeV]", 200, 0.47, 0.52)
    hist_model2 = R.RDF.TH1DModel("", ";M(K_{s});[MeV]", 200, 0.49, 0.505)
    h1 = rdf1.Histo1D(hist_model, "Ks_M")
    h2 = rdf2.Histo1D(hist_model, "Ks_M")
    h2a = rdf2.Histo1D(hist_model2, "Ks_M")
    h3 = rdf3.Histo1D(hist_model, "Ks_M")
    style_draw([h1, h2], "./images/isSignal_before_after.png",["before", "after"], [HistStyle.line_hist(2), HistStyle.line_hist(4)])
    style_draw([h2, h3], "./images/isSignal_bkg_total.png", ["bkg", "total"], [HistStyle.line_hist(4), HistStyle.line_hist(2)], log_y=True)
    print("bkg: ", h2a.Integral()," total: ", h3.Integral()," fraction: ", h2a.Integral()/h3.Integral())

def main():
    rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.3_qqbar_svd2/svd2_only/svd2_4S_offres_st0.root"
    sig_path = "./rootFiles/svd2_4Soffres_st0_v2.3.3_sig.root"
    bkg_path = "./rootFiles/svd2_4Soffres_st0_v2.3.3_bkg.root"
    isSignal_bkg_path = "rootFiles/bkg_isSignal_v2.3.3_with_MCGen.root"
    #process_rootFile(rootFile, sig_path, bkg_path) 
    #isSignal_before_after()
    #bkg_fit(isSignal_bkg_path, "./fit_results/bkg_fit2/")
    #sig_fit(sig_path, "./fit_results/sig_fit/")

    #----- get truth matched MC    
    check_isSignal(rootFile)

if __name__ == "__main__":
    main()



