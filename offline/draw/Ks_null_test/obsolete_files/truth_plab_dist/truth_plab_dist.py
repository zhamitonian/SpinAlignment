import ROOT as R
from DRAW import style_draw, HistStyle

def def_boost_function():
# Add a version that takes px, py, pz, E in CMS
    R.gInterpreter.Declare("""
#include <tuple>
#include <cmath>
#include "TLorentzVector.h"
std::tuple<double, double, double> boostToLab(double p_cms, double costheta_cms, double phi_cms) {
    // Create TLorentzVector in CM frame
    TLorentzVector p4_cms;
    double mass = 0.13957; // pion mass in GeV/c^2
    p4_cms.SetPxPyPzE(
        p_cms * std::sqrt(1 - costheta_cms * costheta_cms) * std::cos(phi_cms),
        p_cms * std::sqrt(1 - costheta_cms * costheta_cms) * std::sin(phi_cms),
        p_cms * costheta_cms,
        std::sqrt(p_cms * p_cms + mass * mass)
    );

    double eH = 7.99638; // HER beam energy
    double eL = 3.49841; // LER beam energy
    double me = 0.51099895069e-3; // electron mass
    double pH = std::sqrt(eH * eH - me * me);
    double pL = std::sqrt(eL * eL - me * me);

    double phiH = 0.022; // HER angle in LAB
    double phiL = 0.0;   // LER angle in LAB
    TLorentzVector p4H(pH * std::sin(phiH), 0, pH * std::cos(phiH), eH);
    TLorentzVector p4L(-pL * std::sin(phiL), 0, -pL * std::cos(phiL), eL);
    TLorentzVector p4CM = p4H + p4L;
    TVector3 betaCM = p4CM.BoostVector();

    // Boost to LAB frame
    p4_cms.Boost(-betaCM);

    double p_lab = p4_cms.P();
    double costheta_lab = p4_cms.Pz() / p_lab;
    double phi_lab = std::atan2(p4_cms.Py(), p4_cms.Px());
    return std::make_tuple(p_lab, costheta_lab, phi_lab);
}

// New version: input px, py, pz, E in CMS
std::tuple<double, double, double> boostToLabVec(double px_cms, double py_cms, double pz_cms, double E_cms) {
    TLorentzVector p4_cms(px_cms, py_cms, pz_cms, E_cms);
    double eH = 7.99638; // HER beam energy
    double eL = 3.49841; // LER beam energy
    double me = 0.51099895069e-3; // electron mass
    double pH = std::sqrt(eH * eH - me * me);
    double pL = std::sqrt(eL * eL - me * me);
    double phiH = 0.022; // HER angle in LAB
    double phiL = 0.0;   // LER angle in LAB
    TLorentzVector p4H(pH * std::sin(phiH), 0, pH * std::cos(phiH), eH);
    TLorentzVector p4L(-pL * std::sin(phiL), 0, -pL * std::cos(phiL), eL);
    TLorentzVector p4CM = p4H + p4L;
    TVector3 betaCM = p4CM.BoostVector();
    p4_cms.Boost(-betaCM);
    double p_lab = p4_cms.P();
    double costheta_lab = p4_cms.Pz() / p_lab;
    double phi_lab = std::atan2(p4_cms.Py(), p4_cms.Px());
    return std::make_tuple(p_lab, costheta_lab, phi_lab);
}

    """)
    a =R.boostToLab(1.0, 0.5, 0.0)
    print(a)

def plot_plab(rootFile, output):
    df = R.RDataFrame("truth", rootFile)
    #df = df.Define("pip_p_lab", "boostToLab(pip_p_truth, pip_costheta_truth, pip_phi_truth).get<0>()")
    #df = df.Define("pip_p_lab", "boostToLabVec(pip_px_cms_truth, pip_py_cms_truth, pip_pz_cms_truth, pip_E_cms_truth).get<0>()") // they are vec branch
    hist = df.Histo1D(("", ";momentum of #pi^{#pm};[MeV]", 120, 0, 3), "pip_p_truth")
    hist2 = df.Histo1D(("", ";momentum of #pi^{#pm};[MeV]", 120, 0, 3), "pim_p_truth")
    hist.GetXaxis().CenterTitle() # not work
    style_draw([hist], output, styles=[HistStyle.filled_line_hist(4, 3011)])
    



if __name__ == "__main__":
    def_boost_function()

    rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.1_qqbar/exp55_4S_offres.root"
    plot_plab(rootFile, "test.png")