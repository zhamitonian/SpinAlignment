"""Data/MC comparison helpers with sWeights and bin slicing."""

import os
import sys
import numpy as np
import ROOT as R

from common.plot.comparison_plotter import plot_data_mc
from common.config.comparison import DEFAULT_KS_VAR_CONFIG, DEFAULT_KS_WEIGHT_CONFIG


class DataMcComparer:
    _cpp_functions_defined = False
    _vec_functions_defined = False
    _daughter_functions_defined = False

    def __init__(self, splot_rootfile, df_data, df_mc, output_dir):
        self.splot_rootfile = splot_rootfile
        self.df_data = df_data
        self.df_mc = df_mc
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def get_splot_map(self):
        sweight_map = {}
        file = R.TFile.Open(self.splot_rootfile, "READ")
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

    def _declare_cpp_functions(self):
        if not DataMcComparer._cpp_functions_defined:
            R.gInterpreter.Declare(
                r"""
                #include <map>
                #include <tuple>
                #include <ROOT/RVec.hxx>
                std::map<std::pair<int, int>, double> g_sweight_map;

                void set_sweight(int event_idx, int cand_idx, double weight) {
                    g_sweight_map[std::make_pair(event_idx, cand_idx)] = weight;
                }

                void clear_sweight_map() { g_sweight_map.clear(); }

                double get_sweight(int event_idx, int cand_idx) {
                    auto key = std::make_pair(event_idx, cand_idx);
                    auto it = g_sweight_map.find(key);
                    if (it != g_sweight_map.end()) return it->second;
                    return 0.0;
                }

                using namespace ROOT::VecOps;
                RVec<double> get_sweight_vec(int event_idx, const RVec<double>& var) {
                    RVec<double> weights;
                    weights.reserve(var.size());
                    for (size_t i = 0; i < var.size(); ++i) {
                        weights.push_back(get_sweight(event_idx, static_cast<int>(i)));
                    }
                    return weights;
                }
                """
            )
            DataMcComparer._cpp_functions_defined = True

    def _declare_daughter_functions(self):
        if not DataMcComparer._daughter_functions_defined:
            R.gInterpreter.Declare(
                r"""
                #include <ROOT/RVec.hxx>
                using namespace ROOT::VecOps;

                RVec<double> extract_daughter_by_index(const RVec<double>& daughter_var,
                                                       const RVec<int>& daughter_index) {
                    RVec<double> result;
                    result.reserve(daughter_index.size());
                    for (size_t i = 0; i < daughter_index.size(); ++i) {
                        int idx = daughter_index[i];
                        if (idx >= 0 && idx < static_cast<int>(daughter_var.size())) {
                            result.push_back(daughter_var[idx]);
                        } else {
                            result.push_back(0.0);
                        }
                    }
                    return result;
                }
                
                RVec<double> get_daughter_sweight(int event_idx,
                                                  const RVec<int>& daughter_index) {
                    RVec<double> weights;
                    weights.reserve(daughter_index.size());
                    for (size_t cand_idx = 0; cand_idx < daughter_index.size(); ++cand_idx) {
                        double sw = get_sweight(event_idx, static_cast<int>(cand_idx));
                        weights.push_back(sw);
                    }
                    return weights;
                }

                """
            )
            DataMcComparer._daughter_functions_defined = True

    def _refresh_sweight_map(self):
        self._declare_cpp_functions()
        R.clear_sweight_map()
        for (evt_idx, cnd_idx), weight in self.sweight_map.items():
            R.set_sweight(evt_idx, cnd_idx, weight)

    def get_data_parent_hist(self, df: R.RDataFrame, label, var_name="phi_M", nbins=60, xmin=1.0, xmax=1.06):
        if not hasattr(self, "sweight_map"):
            raise RuntimeError("sweight_map not initialized. Call get_splot_map() first.")

        self._refresh_sweight_map()

        columns = df.GetColumnNames()
        has_event_idx = "event_idx" in columns
        has_cand_idx = "cand_idx" in columns

        if not has_event_idx:
            print("Warning: event_idx not found in dataframe. Using rdfentry_ as event index.")
            df = df.Define("event_idx", "static_cast<int>(rdfentry_)")

        if not has_cand_idx:
            phi_M_type = df.GetColumnType(var_name)
            if "vector" in phi_M_type or "RVec" in phi_M_type:
                print(f"Warning: {var_name} is a vector but no cand_idx found.")
                print("Assuming each event has candidates indexed 0, 1, 2, ...")
                df = df.Define("nsig_sw", f"get_sweight_vec(static_cast<int>(event_idx), {var_name})")
            else:
                df = df.Define("nsig_sw", "get_sweight(static_cast<int>(event_idx), 0)")
        else:
            df = df.Define("nsig_sw", "get_sweight(static_cast<int>(event_idx), static_cast<int>(cand_idx))")

        model = (var_name, label, nbins, xmin, xmax)
        h = df.Histo1D(model, var_name, "nsig_sw")
        h_result = h.GetValue()
        h_result.SetDirectory(0)
        return h_result

    def get_data_daughter_hist(self, df: R.RDataFrame, label, var_name="kp_p", nbins=60, xmin=0.0, xmax=3.0):
        if not hasattr(self, "sweight_map"):
            raise RuntimeError("sweight_map not initialized. Call get_splot_map() first.")

        self._refresh_sweight_map()
        self._declare_daughter_functions()

        daughter_type = var_name.split("_")[0]
        index_var = f"{daughter_type}_index"

        columns = df.GetColumnNames()
        if "event_idx" not in columns:
            df = df.Define("event_idx", "static_cast<int>(rdfentry_)")

        df = df.Define(f"__{var_name}_indexed__", f"extract_daughter_by_index({var_name}, {index_var})")
        df = df.Define("__daughter_nsig_sw__", f"get_daughter_sweight(static_cast<int>(event_idx), {index_var})")

        model = (var_name, label, nbins, xmin, xmax)
        h = df.Histo1D(model, f"__{var_name}_indexed__", "__daughter_nsig_sw__")
        h_result = h.GetValue()
        h_result.SetDirectory(0)
        return h_result

    def get_MC_hist(self, df: R.RDataFrame, label, var_name="phi_M", nbins=60, xmin=1.0, xmax=1.06, weight=None):
        if weight is not None:
            h_mc = df.Histo1D((var_name + "_mc", label, nbins, xmin, xmax), var_name, weight)
        else:
            h_mc = df.Histo1D((var_name + "_mc", label, nbins, xmin, xmax), var_name)
        h_result = h_mc.GetValue()
        return h_result

    def _plot_pair(self, h_data, h_mc, var, leg_position=2, normalize=True):
        output_path = os.path.join(self.output_dir, f"{var}.png")
        plot_data_mc(h_data, h_mc, output_path, leg_position=leg_position, normalize=normalize)

    def start_plotting(self, var_config=None, weight_config=None):
        if var_config is None:
            var_config = DEFAULT_KS_VAR_CONFIG

        for var, (nbins, xmin, xmax, label) in var_config.items():
            particle_name = var.split("_")[0]
            if particle_name in ["phi", "Ks", "Kstar"]:
                h_data = self.get_data_parent_hist(self.df_data, label, var_name=var, nbins=nbins, xmin=xmin, xmax=xmax)
            else:
                h_data = self.get_data_daughter_hist(self.df_data, label, var_name=var, nbins=nbins, xmin=xmin, xmax=xmax)

            if weight_config is not None and particle_name in weight_config.keys():
                print(f"Applying weight for {particle_name}: {weight_config[particle_name]}")
                h_mc = self.get_MC_hist(self.df_mc, label, var_name=var, nbins=nbins, xmin=xmin, xmax=xmax, weight=weight_config[particle_name])
            else:
                h_mc = self.get_MC_hist(self.df_mc, label, var_name=var, nbins=nbins, xmin=xmin, xmax=xmax)

            self._plot_pair(h_data, h_mc, var, leg_position=2, normalize=True)


def process_single_bin(flat_bin_idx, binning_config, var_config, data_root_dir, splot_root_dir, mc_rootfile, batch_output_dir, mc_weight_config=DEFAULT_KS_WEIGHT_CONFIG):
    bins = [len(binning_list) - 1 for binning_list in binning_config.values()]
    total_bins = int(np.prod(bins))
    pad_width = len(str(total_bins - 1))
    os.makedirs(batch_output_dir, exist_ok=True)

    print("\n" + "=" * 60)
    print(f"Processing bin {flat_bin_idx}/{total_bins - 1}")
    print("=" * 60)

    data_rootfile = os.path.join(data_root_dir, f"temp_bin_{flat_bin_idx}.root")
    splot_rootfile = os.path.join(splot_root_dir, f"bin_{flat_bin_idx:0{pad_width}d}_splot_output.root")

    if not os.path.exists(data_rootfile):
        print(f"Error: Data file not found: {data_rootfile}")
        sys.exit(1)
    if not os.path.exists(splot_rootfile):
        print(f"Error: sPlot file not found: {splot_rootfile}")
        sys.exit(1)

    test_file = R.TFile.Open(data_rootfile, "READ")
    if test_file.IsZombie():
        print(f"Error: Data file is zombie (corrupted): {data_rootfile}")
        test_file.Close()
        sys.exit(1)

    test_tree = test_file.Get("event")
    if not test_tree:
        print(f"Error: 'event' tree not found in {data_rootfile}")
        test_file.Close()
        sys.exit(1)

    n_entries = test_tree.GetEntries()
    test_file.Close()

    if n_entries == 0:
        print(f"Error: Data file has 0 entries: {data_rootfile}")
        sys.exit(1)

    print(f"Data file valid: {n_entries} entries")

    bin_idx = []
    remaining = flat_bin_idx
    for i in range(len(bins)):
        divisor = int(np.prod(bins[i + 1:])) if i + 1 < len(bins) else 1
        bin_idx.append(remaining // divisor)
        remaining = remaining % divisor

    bin_names = list(binning_config.keys())
    bin_str = ", ".join([f"{name.split('_')[-1]}_bin={idx}" for name, idx in zip(bin_names, bin_idx)])
    print(f"Bin indices: {bin_str}")

    bin_conditions = []
    for i, (var_name, bin_edges) in enumerate(binning_config.items()):
        idx = bin_idx[i]
        condition = f"{var_name} >= {bin_edges[idx]} && {var_name} < {bin_edges[idx + 1]}"
        bin_conditions.append(condition)
    bin_condition = " && ".join(bin_conditions)

    print(f"Bin condition: {bin_condition}")

    print(f"Loading MC file: {mc_rootfile}")
    df_mc_full = R.RDataFrame("event", mc_rootfile)

    particle_vars = [var for var in var_config.keys()]
    if mc_weight_config is not None:
        for var in mc_weight_config.values():
            particle_vars.append(var)
    print(particle_vars)

    df_mc_bin = df_mc_full
    for var in particle_vars:
        if var in df_mc_full.GetColumnNames():
            df_mc_bin = df_mc_bin.Redefine(var, f"{var}[{bin_condition}]")

    df_mc_bin = df_mc_bin.Filter(f"Any({bin_condition})")

    df_data = R.RDataFrame("event", data_rootfile)

    bin_output_dir = os.path.join(batch_output_dir, f"bin_{flat_bin_idx:0{pad_width}d}")
    os.makedirs(bin_output_dir, exist_ok=True)

    comparer = DataMcComparer(splot_rootfile, df_data, df_mc_bin, bin_output_dir)
    comparer.get_splot_map()
    print(f"Loaded {len(comparer.sweight_map)} sWeights")

    print("Generating comparison plots...")
    comparer.start_plotting(var_config, mc_weight_config)

    print(f"Bin {flat_bin_idx} completed. Plots saved to {bin_output_dir}")
    print("=" * 60 + "\n")
