import ROOT as R
from DRAW import style_draw , HistStyle
from OFFLINE_PROCESS import RDF_process

"""
get phi multiplicity distribution (in phi mass window: nominal mass +- 3*width) 
for event with at least one generate level phi meson

for detected level: not use gen level info (specifically not use isSignal for Kaon)
"""

PHI_MASS = 1.0195
PHI_WIDTH = 0.004249 

def event_process(df:R.RDataFrame, phi_mass_window=(0.99, 1.15)):
    """
    not same to the phi_reco_truth_match_func in offline/Process.py 
    """
    tools = RDF_process()

    df = df.Filter("n_phi_truth > 0", "has phi") # pre-filter

    # ----------- generation info --------------------
    df = tools.save_all_pairs(df, 
                        p1_branches=["kp_px_cms_gen","kp_py_cms_gen","kp_pz_cms_gen"],
                        p2_branches=["km_px_cms_gen","km_py_cms_gen","km_pz_cms_gen"],
                        mass=(0.493677, 0.493677),
                        particle_name=("phi_gen","kp_gen","km_gen"),
                        cross_mode= False) 
    
    df = df.Define("phi_gen_M", "sqrt(phi_gen_E*phi_gen_E - phi_gen_px*phi_gen_px - phi_gen_py*phi_gen_py - phi_gen_pz*phi_gen_pz)")
    """
    mass_condition = f"(phi_gen_M > {phi_mass_window[0]}) && (phi_gen_M < {phi_mass_window[1]})"
    df = df.Redefine("phi_gen_E", f"phi_gen_E[{mass_condition}]")
    df = df.Redefine("phi_gen_px", f"phi_gen_px[{mass_condition}]")
    df = df.Redefine("phi_gen_py", f"phi_gen_py[{mass_condition}]")
    df = df.Redefine("phi_gen_pz", f"phi_gen_pz[{mass_condition}]")
    df = df.Redefine("phi_gen_helicity_angle", f"phi_gen_helicity_angle[{mass_condition}]")
    df = df.Redefine("phi_gen_M", f"phi_gen_M[{mass_condition}]")
    """
    df = tools.calculate_pt_toAxis(df, 
                        particle=["phi_gen_px","phi_gen_py","phi_gen_pz"],
                        axis = ["TMath::ACos(thrust_truth[1])", "thrust_truth[2]"],
                        particle_name="phi_gen",axis_name="thrust")

    df = df.Redefine("thrust_truth", "thrust_truth[0]")
    df = df.Define("phi_gen_z", f"2*phi_gen_E/sqrts")
    df = df.Define("phi_gen_xp", "sqrt(phi_gen_E*phi_gen_E - phi_gen_M*phi_gen_M)/sqrt(sqrts*sqrts*0.25 - 1.0195*1.0195)")

    df = df.Filter("phi_gen_M.size() > 0", "has phi in mass window at gen level")    
    # has phi in mass window at gen level: pass=3196823    all=4783173    -- eff=66.83 % 

    # ----------- reconstruction info --------------------
    df = tools.save_all_pairs(df, 
                         p1_branches=["kp_px_cms","kp_py_cms","kp_pz_cms"],
                         p2_branches=["km_px_cms","km_py_cms","km_pz_cms"],
                         mass=(0.493677, 0.493677),
                         particle_name=("phi","kp","km"))

    df = df.Define("phi_M", "sqrt(phi_E*phi_E - phi_px*phi_px - phi_py*phi_py - phi_pz*phi_pz)")
    
    
    mass_condition = f"(phi_M > {phi_mass_window[0]}) && (phi_M < {phi_mass_window[1]})"
    # Apply mass window filter to all phi-related variables
    df = df.Redefine("phi_E", f"phi_E[{mass_condition}]")
    df = df.Redefine("phi_px", f"phi_px[{mass_condition}]")
    df = df.Redefine("phi_py", f"phi_py[{mass_condition}]")
    df = df.Redefine("phi_pz", f"phi_pz[{mass_condition}]")
    df = df.Redefine("phi_helicity_angle", f"phi_helicity_angle[{mass_condition}]")
    df = df.Redefine("kp_index", f"kp_index[{mass_condition}]")
    df = df.Redefine("km_index", f"km_index[{mass_condition}]")
    df = df.Redefine("phi_M", f"phi_M[{mass_condition}]")
    
    df = df.Define("phi_z", "2*phi_E/sqrts")
    df = df.Define("phi_xp", "sqrt(phi_E*phi_E - phi_M*phi_M)/sqrt(sqrts*sqrts*0.25 - 1.0195*1.0195)")
    df = tools.calculate_pt_toAxis(df, 
                                  particle=["phi_px","phi_py","phi_pz"],
                                  axis = ["TMath::ACos(thrust[1])", "thrust[2]"],
                                  particle_name="phi",axis_name="thrust")
    for kaon in ["kp", "km"]:
        df = tools.calculate_pt_toAxis(df,
                                    particle=[f"{kaon}_px_cms", f"{kaon}_py_cms",f"{kaon}_pz_cms"],
                                    axis = ["TMath::ACos(thrust[1])", "thrust[2]"],
                                    particle_name=kaon,axis_name="thrust")

    df = df.Redefine("thrust", "thrust[0]")


    Br2Save = ["phi_M", "phi_gen_M", "phi_z", "phi_gen_z", "phi_xp", "phi_gen_xp", "phi_helicity_angle", "phi_gen_helicity_angle", "phi_thrust_pt"]
    #Br2Save += ["phi_costheta", "phi_phi", "phi_p", "kp_p", "kp_costheta", "kp_phi", "km_p", "km_costheta", "km_phi"]

    df.Report().Print()

    return df

def phi_multi():
    """
    test the mc match truth efficiency
    """
    gMC = "/gpfs/group/belle/users/wangz/data_gMC_belle1/PhiSpinAlignment_v2.1.2_qqbar/continuum.root"
    df = R.RDataFrame("event", gMC)
    df = event_process(df, (PHI_MASS - 3*PHI_WIDTH, PHI_MASS + 3*PHI_WIDTH))

    hist_model = R.RDF.TH1DModel("", "#phi multiplicity", 6, -0.5, 5.5)

    hist = df.Define("phi_multplicity", "phi_M.size()").Histo1D(hist_model, "phi_multplicity").GetValue()
    hist_gen = df.Define("phi_gen_multplicity", "phi_gen_M.size()").Histo1D(hist_model, "phi_gen_multplicity").GetValue()
    style_draw([hist, hist_gen],"./images/phi_multiplicity.png", ["detect level", "gen level"], [HistStyle.error_bars(2), HistStyle.error_bars(4)], log_y=True)

    h_phi_M_gen = df.Histo1D(("", "#phi mass gen level;Mass (GeV/c^{2})", 100, 0.99, 1.15), "phi_gen_M").GetValue()
    style_draw([h_phi_M_gen], "test.png", ["gen level"], [HistStyle.line_hist(4)])

    df = df.Filter("thrust > 0.8", "thrust > 0.8")
    
    # with additional requirement at detected level : K+ K- in same semisphere
    # Define a helper function to filter phi_M based on same-semisphere condition
    R.gInterpreter.Declare("""
        #include <ROOT/RVec.hxx>
        using namespace ROOT::VecOps;
        
        RVec<double> filter_phi_M_same_semisphere(
            const RVec<double>& phi_M,
            const RVec<int>& kp_index,
            const RVec<int>& km_index,
            const RVec<double>& kp_thrust_costheta,
            const RVec<double>& km_thrust_costheta) 
        {
            RVec<double> result;
            for (size_t i = 0; i < phi_M.size(); ++i) {
                if (kp_thrust_costheta[kp_index[i]] * km_thrust_costheta[km_index[i]] > 0) {
                    result.push_back(phi_M[i]);
                }
            }
            return result;
        }
    """)

    R.gInterpreter.Declare("""
        #include <ROOT/RVec.hxx>
        using namespace ROOT::VecOps;
        
        RVec<double> filter_phi_M_unsame_semisphere(
            const RVec<double>& phi_M,
            const RVec<int>& kp_index,
            const RVec<int>& km_index,
            const RVec<double>& kp_thrust_costheta,
            const RVec<double>& km_thrust_costheta) 
        {
            RVec<double> result;
            for (size_t i = 0; i < phi_M.size(); ++i) {
                if (kp_thrust_costheta[kp_index[i]] * km_thrust_costheta[km_index[i]] < 0) {
                    result.push_back(phi_M[i]);
                }
            }
            return result;
        }
    """)
    
    df_ss = df.Define("phi_M_ss", "filter_phi_M_same_semisphere(phi_M, kp_index, km_index, kp_thrust_costheta, km_thrust_costheta)")
    hist_detect_ss = df_ss.Define("phi_multplicity_ss", "phi_M_ss.size()").Histo1D(hist_model, "phi_multplicity_ss").GetValue()

    df_nss = df.Define("phi_M_nss", "filter_phi_M_unsame_semisphere(phi_M, kp_index, km_index, kp_thrust_costheta, km_thrust_costheta)")
    hist_detect_nss = df_nss.Define("phi_multplicity_nss", "phi_M_nss.size()").Histo1D(hist_model, "phi_multplicity_nss").GetValue()

    style_draw([hist_detect_ss, hist_detect_nss, hist_gen], 
               "./images/phi_multiplicity_semi_sphere.png", 
               ["Same hemisphere", "Opposite hemisphere", "Gen level"], 
               [HistStyle.error_bars(2), HistStyle.error_bars(4), HistStyle.error_bars(1)], 
               log_y=True)


if __name__ == "__main__":
    phi_multi()