#!/usr/bin/env python3

"""
Generic data-MC comparison tool for rho00 measurements
Compares sWeighted data with truth-matched MC in N-dimensional bins
Supports multiple particle types (phi, Ks, etc.) with flexible binning
"""

import ROOT as R
import os
from DRAW import style_draw, HistStyle
import numpy as np

class comparing:
    _cpp_functions_defined = False  # Class variable to track if C++ functions are defined
    
    def __init__(self, splot_rootFile, df_data, df_mc, output_dir):
        self.splot_rootFile = splot_rootFile
        self.df_data = df_data
        self.df_mc = df_mc
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        
    def get_splot_map(self):
        sweight_map = {}
        
        file = R.TFile.Open(self.splot_rootFile, "READ")
        tree = file.Get("RooTreeDataStore_sWeighted_data_Combined_dataset")
        
        for i in range(tree.GetEntries()):
            tree.GetEntry(i)
            evt_idx = int(tree.event_idx)
            cnd_idx = int(tree.cand_idx)
            nsig_sw = tree.nsig_sw
            sweight_map[(evt_idx, cnd_idx)] = nsig_sw
            if i < 5:
                print(f"sWeighted evt_idx: {evt_idx}, cnd_idx: {cnd_idx}, nsig_sw: {nsig_sw}")

        self.sweight_map = sweight_map
        file.Close()

    def get_data_parent_hist(self, df:R.RDataFrame, label, var_name="phi_M", nbins=60, xmin=1.0, xmax=1.06):
        """
        Get histogram of parent particle data with sWeights applied
        
        Parameters:
        -----------
        df : R.RDataFrame
            Input dataframe with event_idx and cand_idx columns (from vector branches)
        label : str
            Histogram label for axis
        var_name : str
            Parent particle variable name to plot (e.g., "phi_M", "Ks_p")
        nbins, xmin, xmax : int, float, float
            Histogram binning parameters
            
        Returns:
        --------
        R.TH1D : Weighted histogram of parent particle variable
        """
        if not hasattr(self, 'sweight_map'):
            raise RuntimeError("sweight_map not initialized. Call get_splot_map() first.")
        
        # Define C++ functions only once
        if not comparing._cpp_functions_defined:
            R.gInterpreter.Declare("""
                #include <map>
                #include <tuple>
                std::map<std::pair<int, int>, double> g_sweight_map;
                
                void set_sweight(int event_idx, int cand_idx, double weight) {
                    g_sweight_map[std::make_pair(event_idx, cand_idx)] = weight;
                }
                
                void clear_sweight_map() {
                    g_sweight_map.clear();
                }
                
                double get_sweight(int event_idx, int cand_idx) {
                    auto key = std::make_pair(event_idx, cand_idx);
                    auto it = g_sweight_map.find(key);
                    if (it != g_sweight_map.end()) {
                        return it->second;
                    }
                    return 0.0;  // Return 0 if not found
                }
            """)
            comparing._cpp_functions_defined = True
        
        # Clear and refill the map (in case sweight_map has changed)
        R.clear_sweight_map()
        
        # Fill the C++ map with Python dictionary
        for (evt_idx, cnd_idx), weight in self.sweight_map.items():
            R.set_sweight(evt_idx, cnd_idx, weight)
        
        # Check if dataframe has event_idx column
        columns = df.GetColumnNames()
        has_event_idx = "event_idx" in columns
        has_cand_idx = "cand_idx" in columns
        
        # If no event_idx, use rdfentry_ as event index
        if not has_event_idx:
            print("Warning: event_idx not found in dataframe. Using rdfentry_ as event index.")
            df = df.Define("event_idx", "static_cast<int>(rdfentry_)")
        
        # For vector branches, we need cand_idx
        # If phi_M is a vector, we need to handle it differently
        if not has_cand_idx:
            # Check if phi_M is a vector
            phi_M_type = df.GetColumnType(var_name)
            if "vector" in phi_M_type or "RVec" in phi_M_type:
                print(f"Warning: {var_name} is a vector but no cand_idx found.")
                print("Assuming each event has candidates indexed 0, 1, 2, ...")
                
                # Define vector version only once
                if not hasattr(comparing, '_vec_functions_defined'):
                    R.gInterpreter.Declare("""
                        #include <ROOT/RVec.hxx>
                        using namespace ROOT::VecOps;
                        RVec<double> get_sweight_vec(int event_idx, const RVec<double>& var) {
                            RVec<double> weights;
                            weights.reserve(var.size());
                            for (size_t i = 0; i < var.size(); ++i) {
                                weights.push_back(get_sweight(event_idx, static_cast<int>(i)));
                            }
                            return weights;
                        }
                    """)
                    comparing._vec_functions_defined = True
                    
                df = df.Define("nsig_sw", f"get_sweight_vec(static_cast<int>(event_idx), {var_name})")
            else:
                # Scalar case: assume cand_idx = 0
                df = df.Define("nsig_sw", "get_sweight(static_cast<int>(event_idx), 0)")
        else:
            # Both event_idx and cand_idx exist
            df = df.Define("nsig_sw", "get_sweight(static_cast<int>(event_idx), static_cast<int>(cand_idx))")
        
        # Create weighted histogram
        model = (var_name, label, nbins, xmin, xmax)
        h = df.Histo1D(model, var_name, "nsig_sw")
        
        # Must call GetValue() or trigger execution
        h_result = h.GetValue()
        h_result.SetDirectory(0)  # Detach from file
        
        return h_result

    def get_data_daughter_hist(self, df:R.RDataFrame, label, var_name="kp_p", 
                     nbins=60, xmin=0.0, xmax=3.0):
        """
        Get histogram of daughter particle variables with sWeights applied through parent index
        
        For daughter particle variables (kp_p, km_p, pip_p, pim_p, etc.), we:
        1. Use daughter_index to map daughter -> parent candidate
        2. Apply the sweight of the corresponding parent candidate
        
        Parameters:
        -----------
        df : R.RDataFrame
            Input dataframe with daughter variables and corresponding index branches
        label : str
            Histogram label for axis
        var_name : str
            Daughter particle variable name (e.g., "kp_p", "pip_costheta")
        nbins, xmin, xmax : int, float, float
            Histogram binning parameters
            
        Returns:
        --------
        R.TH1D : Weighted histogram of daughter particle variable
        """
        if not hasattr(self, 'sweight_map'):
            raise RuntimeError("sweight_map not initialized. Call get_splot_map() first.")
        
        # Define C++ functions only once (same as get_data_hist)
        if not comparing._cpp_functions_defined:
            R.gInterpreter.Declare("""
                #include <map>
                #include <tuple>
                std::map<std::pair<int, int>, double> g_sweight_map;
                
                void set_sweight(int event_idx, int cand_idx, double weight) {
                    g_sweight_map[std::make_pair(event_idx, cand_idx)] = weight;
                }
                
                void clear_sweight_map() {
                    g_sweight_map.clear();
                }
                
                double get_sweight(int event_idx, int cand_idx) {
                    auto key = std::make_pair(event_idx, cand_idx);
                    auto it = g_sweight_map.find(key);
                    if (it != g_sweight_map.end()) {
                        return it->second;
                    }
                    return 0.0;
                }
            """)
            comparing._cpp_functions_defined = True
        
        # Define daughter-specific weight function only once
        if not hasattr(comparing, '_daughter_functions_defined'):
            R.gInterpreter.Declare("""
                #include <ROOT/RVec.hxx>
                using namespace ROOT::VecOps;
                
                // Apply sweight to daughter variables using parent candidate index
                RVec<double> get_daughter_sweight(int event_idx, 
                                                  const RVec<double>& daughter_var,
                                                  const RVec<int>& daughter_index) {
                    RVec<double> weights;
                    weights.reserve(daughter_index.size());
                    
                    // For each parent candidate (cand_idx)
                    for (size_t cand_idx = 0; cand_idx < daughter_index.size(); ++cand_idx) {
                        // Get the sweight for this parent candidate
                        double sw = get_sweight(event_idx, static_cast<int>(cand_idx));
                        weights.push_back(sw);
                    }
                    
                    return weights;
                }
                
                // Extract daughter values using index mapping
                RVec<double> extract_daughter_by_index(const RVec<double>& daughter_var,
                                                       const RVec<int>& daughter_index) {
                    RVec<double> result;
                    result.reserve(daughter_index.size());
                    
                    for (size_t i = 0; i < daughter_index.size(); ++i) {
                        int idx = daughter_index[i];
                        if (idx >= 0 && idx < static_cast<int>(daughter_var.size())) {
                            result.push_back(daughter_var[idx]);
                        } else {
                            result.push_back(0.0);  // Invalid index
                        }
                    }
                    
                    return result;
                }
            """)
            comparing._daughter_functions_defined = True
        
        # Clear and refill the map
        R.clear_sweight_map()
        for (evt_idx, cnd_idx), weight in self.sweight_map.items():
            R.set_sweight(evt_idx, cnd_idx, weight)
        
        # Determine index variable name
        daughter_type = var_name.split("_")[0]  # e.g., "kp" from "kp_p"
        index_var = f"{daughter_type}_index"
        
        # Check if we need to use event_idx
        columns = df.GetColumnNames()
        if "event_idx" not in columns:
            df = df.Define("event_idx", "static_cast<int>(rdfentry_)")
        
        # Extract daughter values by index and apply weights
        df = df.Define(f"__{var_name}_indexed__", 
                      f"extract_daughter_by_index({var_name}, {index_var})")
        df = df.Define("__daughter_nsig_sw__", 
                      f"get_daughter_sweight(static_cast<int>(event_idx), {var_name}, {index_var})")
        
        # Create weighted histogram
        model = (var_name, label, nbins, xmin, xmax)
        h = df.Histo1D(model, f"__{var_name}_indexed__", "__daughter_nsig_sw__")
        
        h_result = h.GetValue()
        h_result.SetDirectory(0)
        
        return h_result


    def get_MC_hist(self, df:R.RDataFrame, label, var_name="phi_M", nbins=60, xmin=1.0, xmax=1.06):      
        """
        Get histogram from MC data (unweighted)
        
        Parameters:
        -----------
        df : R.RDataFrame
            MC dataframe (already filtered for bin)
        label : str
            Histogram label for axis
        var_name : str
            Variable name to plot
        nbins, xmin, xmax : int, float, float
            Histogram binning parameters
            
        Returns:
        --------
        R.TH1D : MC histogram
        """
        h_MC = df.Histo1D((var_name+"_mc", label, nbins, xmin, xmax), var_name)
        h_result = h_MC.GetValue()
        return h_result

    def brush(self, hist_list, var, leg_position=2, normailize=False):
        h_data = hist_list[0]
        h_mc = hist_list[1]
        if normailize:
            h_mc.Scale(h_data.Integral() / h_mc.Integral())

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

    def start_plotting(self, var_config):
        if var_config is None:
            var_config = {"phi_M": (60, 1.0, 1.06, ";M_{#phi};[MeV]"), 
                        "phi_p": (100, 0, 5, ";p(#phi);[GeV/c]"),
                        "phi_costheta": (20, -1, 1, ";cos#theta(#phi);[]"),
                        "phi_phi" : (20, -3.14, 3.14, ";#varphi(#phi);[rad]"),
                        "kp_p": (90, 0, 4.5, ";p(K^{+});[MeV]"),
                        "km_p": (90, 0, 4.5, ";p(K^{-});[MeV]"),
                        "kp_costheta": (20, -1, 1, ";cos#theta({K^{+});[]"),
                        "km_costheta": (20, -1, 1, ";cos#theta(K^{-});[]"),
                        "kp_phi": (20, -3.14, 3.14, ";#phi(K^{+});[rad]"),
                        "km_phi": (20, -3.14, 3.14, ";#phi({K^{-});[rad]"),
                        }

        for var, (nbins, xmin, xmax, label) in var_config.items():
            #if var.startswith("phi_"):
            if var.split("_")[0] in ["phi", "Ks", "Kstar"]:
                h_data = self.get_data_parent_hist(self.df_data, label, var_name=var, nbins=nbins, xmin=xmin, xmax=xmax)
            else:
                h_data = self.get_data_daughter_hist(self.df_data, label, var_name=var, nbins=nbins, xmin=xmin, xmax=xmax)
            
            h_mc = self.get_MC_hist(self.df_mc, label, var_name=var, nbins=nbins, xmin=xmin, xmax=xmax)

            self.brush([h_data, h_mc], var, leg_position=2, normailize=True)

        

def process_single_bin(flat_bin_idx, binning_config, var_config, data_rootFile_dir, splot_rootFile_dir, mc_rootFile, batch_output_dir):
    """
    Process a single bin for data/MC comparison
    Designed to be called by batch jobs (LSF, condor, etc.)
    
    This function loads data and MC for a specific N-dimensional bin,
    applies sWeights from sPlot, and generates comparison plots.
    
    Parameters:
    -----------
    flat_bin_idx : int
        Flat bin index (0 to total_bins-1) for N-dimensional binning
    binning_config : dict
        Dictionary mapping binning variable names to bin edges
        e.g., {"phi_z": [0.2, 0.3, ...], "phi_thrust_pt": [...]}
    data_rootFile_dir : str
        Directory containing data ROOT files
    splot_rootFile_dir : str
        Directory containing sPlot ROOT files
    mc_rootFile : str
        Path to MC ROOT file
    batch_output_dir : str
        Output directory for comparison plots
    """
    import sys
    
    # Calculate bin dimensions and total bins
    bins = [len(binning_list) - 1 for binning_list in binning_config.values()]
    total_bins = int(np.prod(bins))  # Works for any number of dimensions
    pad_width = len(str(total_bins - 1))
    os.makedirs(batch_output_dir, exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"Processing bin {flat_bin_idx}/{total_bins-1}")
    print(f"{'='*60}")
    
    # File paths for this bin
    data_rootFile = os.path.join(data_rootFile_dir, f"temp_bin_{flat_bin_idx}.root")
    splot_rootFile = os.path.join(splot_rootFile_dir, f"bin_{flat_bin_idx:0{pad_width}d}_splot_output.root")
    
    # Check if files exist
    if not os.path.exists(data_rootFile):
        print(f"Error: Data file not found: {data_rootFile}")
        sys.exit(1)
    if not os.path.exists(splot_rootFile):
        print(f"Error: sPlot file not found: {splot_rootFile}")
        sys.exit(1)
    
    # Check if data ROOT file is valid (not zombie and has entries)
    test_file = R.TFile.Open(data_rootFile, "READ")
    if test_file.IsZombie():
        print(f"Error: Data file is zombie (corrupted): {data_rootFile}")
        test_file.Close()
        sys.exit(1)
    
    test_tree = test_file.Get("event")
    if not test_tree:
        print(f"Error: 'event' tree not found in {data_rootFile}")
        test_file.Close()
        sys.exit(1)
    
    n_entries = test_tree.GetEntries()
    test_file.Close()
    
    if n_entries == 0:
        print(f"Error: Data file has 0 entries: {data_rootFile}")
        sys.exit(1)
    
    print(f"Data file valid: {n_entries} entries")
    
    # Calculate N-dimensional bin indices generically
    bin_idx = []
    remaining = flat_bin_idx
    for i in range(len(bins)):  # Go forward through dimensions
        divisor = int(np.prod(bins[i+1:])) if i+1 < len(bins) else 1
        bin_idx.append(remaining // divisor)
        remaining = remaining % divisor
    
    # Print bin indices with variable names
    bin_names = list(binning_config.keys())
    bin_str = ", ".join([f"{name.split('_')[-1]}_bin={idx}" for name, idx in zip(bin_names, bin_idx)])
    print(f"Bin indices: {bin_str}")
    
    # Construct bin condition dynamically from binning_config
    bin_conditions = []
    for i, (var_name, bin_edges) in enumerate(binning_config.items()):
        idx = bin_idx[i]
        condition = f"{var_name} >= {bin_edges[idx]} && {var_name} < {bin_edges[idx+1]}"
        bin_conditions.append(condition)
    bin_condition = " && ".join(bin_conditions)

    print(f"Bin condition: {bin_condition}")

    # Load and filter MC for this bin
    print(f"Loading MC file: {mc_rootFile}")
    df_mc_full = R.RDataFrame("event", mc_rootFile)
    
    # Extract particle name from binning_config (e.g., "phi" from "phi_z", "Ks" from "Ks_z")
    first_var = list(binning_config.keys())[0]
    particle_name = "_".join(first_var.split("_")[:-1])  # Remove last component (e.g., "z", "thrust_pt")
    
    # Dynamically redefine particle variables based on particle name
    particle_vars = [f"{particle_name}_p", f"{particle_name}_costheta", f"{particle_name}_phi"]
    df_mc_bin = df_mc_full
    for var in particle_vars:
        # Check if variable exists in dataframe before redefining
        if var in df_mc_full.GetColumnNames():
            df_mc_bin = df_mc_bin.Redefine(var, f"{var}[{bin_condition}]")
    
    df_mc_bin = df_mc_bin.Filter(f"Any({bin_condition})")
    ''' same to below
    df_mc_bin = df_mc_full.Redefine("phi_p", f"phi_p[{bin_condition}]") \
                          .Redefine("phi_costheta", f"phi_costheta[{bin_condition}]") \
                          .Redefine("phi_phi", f"phi_phi[{bin_condition}]") \
                          .Filter(f"Any({bin_condition})")
    '''

    # Load data for this bin
    df_data = R.RDataFrame("event", data_rootFile)
    
    # Create output directory for this bin
    bin_output_dir = os.path.join(batch_output_dir, f"bin_{flat_bin_idx:0{pad_width}d}")
    os.makedirs(bin_output_dir, exist_ok=True)
    
    # Create comparing instance
    comparer = comparing(splot_rootFile, df_data, df_mc_bin, bin_output_dir)
    
    # Load sPlot weights
    comparer.get_splot_map()
    print(f"Loaded {len(comparer.sweight_map)} sWeights")
    
    # Perform comparison and generate plots
    print(f"Generating comparison plots...")
    comparer.start_plotting(var_config)
    
    print(f"Bin {flat_bin_idx} completed. Plots saved to {bin_output_dir}")
    print(f"{'='*60}\n")

