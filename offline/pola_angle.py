import ROOT
from typing import Tuple, Optional

def get_polarization_angles(df: ROOT.RDataFrame, 
                            p1_branches : Tuple[str, str, str],
                            p2_branches : Tuple[str, str, str],
                            mass: Tuple[float, float], axis: Tuple[str, str],
                            output_names: Tuple[str, str] = ("cos_theta", "phi")) -> ROOT.RDataFrame:
    """
    Calculate the polarization angles for a given DataFrame.
    """

    new_df = df

    ROOT.gInterpretation.Declare("""
        #include <ROOT/RVec.hxx>
        #include <TLorentzVector.h>
        using namespace ROOT::VecOps;
        
        std::tuple<double, double> 
        calculate_angles(
            double p1_px, double p1_py, double p1_pz,
            double p2_px, double p2_py, double p2_pz,
            double mass_p1, double mass_p2,
            double axis_costheta, double axis_phi)
        {
            TLorentzVector p1,p2;
            p1.SetXYZM(p1_px, p1_py, p1_pz, mass_p1);
            p2.SetXYZM(p2_px, p2_py, p2_pz, mass_p2);
            TLorentzVector parent = p1 + p2;

            // --- define coordinate system in parent rest frame ---
            // z: parent flight direction
            TVector3 z_hat = parent.Vect().Unit();

            // reference axis from spherical coordinates (axis_costheta, axis_phi)
            double axis_sintheta = sqrt(1.0 - axis_costheta * axis_costheta);
            TVector3 axis_vec(axis_sintheta * cos(axis_phi),
                              axis_sintheta * sin(axis_phi),
                              axis_costheta);

            // y = z x axis  (normal to the plane spanned by z and axis)
            TVector3 y_hat = z_hat.Cross(axis_vec).Unit();

            // x = y x z  (right-handed: x x y = z)
            TVector3 x_hat = y_hat.Cross(z_hat);

            // --- boost p1 into parent rest frame ---
            TLorentzVector p1_rest = p1;
            p1_rest.Boost(-parent.BoostVector());
            TVector3 p1_vec = p1_rest.Vect();

            // --- project onto new axes ---
            double cos_theta = p1_vec.Dot(z_hat) / p1_vec.Mag();
            double phi       = atan2(p1_vec.Dot(y_hat), p1_vec.Dot(x_hat));

            return std::make_tuple(cos_theta, phi);
        }

        std::tuple<RVec<double>, RVec<double>>
        calculate_angles(
            RVec<double> p1_px, RVec<double> p1_py, RVec<double> p1_pz,
            RVec<double> p2_px, RVec<double> p2_py, RVec<double> p2_pz,
            RVec<double> mass_p1, RVec<double> mass_p2,
            RVec<double> axis_costheta, RVec<double> axis_phi)
        {
            size_t n = p1_px.size();
            RVec<double> cos_theta(n);
            RVec<double> phi(n);
            for (size_t i = 0; i < n; ++i) {
                std::tie(cos_theta[i], phi[i]) = calculate_angles(
                    p1_px[i], p1_py[i], p1_pz[i],
                    p2_px[i], p2_py[i], p2_pz[i],
                    mass_p1[i], mass_p2[i],
                    axis_costheta[i], axis_phi[i]
                );
            }
            return std::make_tuple(cos_theta, phi);
        }
        """)

    if "Pairs" not in df.GetColumnNames():
        new_df = new_df.Define("Pairs", f"calculate_angles({p1_branches[0]}, {p1_branches[1]}, {p1_branches[2]}, {p2_branches[0]}," 
                                        f"{p2_branches[1]}, {p2_branches[2]}, {mass[0]}, {mass[1]}, {axis[0]}, {axis[1]})")
    else: 
        print("Warning: 'Pairs' column already exists. Overwriting with new angles.")
        new_df = new_df.Redefine("Pairs", f"calculate_angles({p1_branches[0]}, {p1_branches[1]}, {p1_branches[2]}, {p2_branches[0]}," 
                                        f"{p2_branches[1]}, {p2_branches[2]}, {mass[0]}, {mass[1]}, {axis[0]}, {axis[1]})")
    
    if output_names[0] in df.GetColumnNames() or output_names[1] in df.GetColumnNames():
        print(f"Warning: Output columns {output_names} already exist. The calculation will not be written.")
        return df
    else :
        new_df = new_df.Define(output_names[0], "Pairs.get<0>()")
        new_df = new_df.Define(output_names[1], "Pairs.get<1>()")

    return new_df