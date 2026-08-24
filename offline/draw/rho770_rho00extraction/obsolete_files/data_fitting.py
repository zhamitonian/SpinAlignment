#!/usr/bin/env python3
from FIT import QUICK_FIT
import ROOT as R
import os
from math import pi

from FIT.generic_fit import *

def wrapper_func(tree, output_dir, log_file, bin_fit_range, binned_fit):

    ranges = bin_fit_range.split(';')
    range_txt = []
    for range_str in ranges:
        range_str = range_str.strip()
        range_txt.append(range_str)
    bin_index = int(log_file.split("_")[-1].split(".")[0])

    # data fit
    n_entry = tree.GetEntries()
    #----------------------------------------------
    m_pi = 0.13957
    m_rho = 0.77526
    m_omega = 0.78265
    Gamma_rho_0 = 0.1491
    Gamma_omega_0 = 0.00849

    R.gInterpreter.Declare(f"""
#ifndef RHO_OMEGA_SIGNAL_DEFINED
#define RHO_OMEGA_SIGNAL_DEFINED
#include <complex>
Double_t rho_omega_signal(Double_t m, Double_t A, Double_t B, Double_t C, Double_t phase, 
                           Double_t Gamma_rho = {Gamma_rho_0}, Double_t Gamma_omega = {Gamma_omega_0},
                           Double_t m_rho = {m_rho}, Double_t m_omega = {m_omega}) {{
    using cd = std::complex<double>;
    const double m2 = m * m;

    // mass-dependent widths
    double Gr = Gamma_rho * (m_rho/m) * std::pow((m2 - {4*m_pi**2}) / (m_rho*m_rho - {4*m_pi**2}), 1.5);
    double Gw = Gamma_omega * (m_omega/m) * std::pow((m2 - {4*m_pi**2}) / (m_omega*m_omega - {4*m_pi**2}), 1.5);

    // Breit-Wigner amplitudes: sqrt(m*m0*Gamma) / (m0^2 -m^2 - i*m0*Gamma)
    cd bw_rho = std::sqrt(m * m_rho * Gr) / cd(m_rho*m_rho - m2, -m_rho * Gr);
    cd bw_omega = std::sqrt(m * m_omega * Gw) / cd(m_omega*m_omega - m2, -m_omega * Gw);

    cd amp = A * bw_rho + B + C * std::exp(cd(0, 1) * phase) * bw_omega;
    return std::norm(amp);   // |amp|^2
}}
#endif
""")

    pdf_config = FitDefinition([Variable("rho_M", 0.770 - 0.150, 0.770 + 0.150, 120)],
                [PDFSpec("Sig", "rho_M", "generic_pdf",
                        config = {"formula": f"rho_omega_signal(rho_M, A[1,0,10], B[0.1,-5,5], C[1,0,10], phase[{pi},0,{2*pi}])"} ) ,
                                  #Gamma_omega = Gamma_omega[{Gamma_omega_0}, 0.85*{Gamma_omega_0}, 1.15*{Gamma_omega_0}], \
                                  #m_omega = m_omega[{m_omega}, 0.85*{m_omega}, 1.15*{m_omega}] )"}) ,
                                  #Gamma_rho = Gamma_rho[{Gamma_rho_0}, 0.85*{Gamma_rho_0}, 1.15*{Gamma_rho_0}], \
                                  #m_rho = m_rho[{m_rho}, 0.85*{m_rho}, 1.15*{m_rho}] )"}) ,
                                  #Gamma_rho[{Gamma_rho_0}, 0.85*{Gamma_rho_0}, 1.15*{Gamma_rho_0}], Gamma_omega[{Gamma_omega_0}, 0.85*{Gamma_omega_0}, 1.15*{Gamma_omega_0}], \
                                  #m_rho[{m_rho}, 0.85*{m_rho}, 1.15*{m_rho}], m_omega[{m_omega}, 0.85*{m_omega}, 1.15*{m_omega}] )"}) ,
                #PDFSpec("bkg", "rho_M", "chebychev", {"order":3, "coef1": (-0.2, -0.3, 0)}),],
                PDFSpec("bkg", "rho_M", "chebychev", {"order":3,}),],
                model = f"SUM(nsig[{0.15 * n_entry}, 0, {0.3*n_entry}] * Sig, nbkg[{0.85 * n_entry},{0.5* n_entry}, {n_entry}]*bkg)")

    plot_config = PlotConfiguration(plot_config={
        "xlabel" : {"rho_M" : "M_{#pi^{+}#pi^{-}} (GeV/c^{2})"},
        "components" : {"model" : {"label" : "Total Fit", "color" : 4},
                        "sig": {"label" : "Signal", "color" : 2, "style" : 4, "width" : 3},
                        "bkg": {"label" : "Background", "color" : R.kAzure + 1, "style" : 3, "width" : 3},},
        "legend" : {"extra_text" : range_txt,}, 
        "show_pull" :False, "logy":False})
    dataset_config = DatasetConfig(binned_fit= binned_fit, target_branch=["rho_M"], perform_splot=False)
    fit_config = FitterConfig(two_step_fit=True, use_minos=False)

    fitter = GenericFit(tree, output_dir, log_file=log_file,fit_definition=pdf_config,
                    dataset_config=dataset_config, plot_config=plot_config, fitter_config=fit_config)
    fitter.run()
    result_data, yield_results = fitter.run()
    nsig, nsig_err = yield_results["nsig"], yield_results["nsig_err"]

    return result_data , nsig, nsig_err 


def fit(): 
    quick_fit = QUICK_FIT(wrapper_func, {"rho_z": (20, 0, 1), "rho_helicity_angle": (10, -1, 1)})
    quick_fit.parse_arguments()

if __name__ == "__main__":
    fit()