#!/usr/bin/env python3
import numpy as np
import pandas as pd

import ROOT as R

CPP_FUNC_DEFINED = False

def get_pid_eff(rdf:R.RDataFrame, period, is_pion=True, cut_value=6, particle_names = ("pip", "pim"), eff_source="data") -> R.RDataFrame:
    """
    Add PID efficiency branches for the selected efficiency source.
    
    :param rdf: RDataFrame
    :type rdf: R.RDataFrame
    :param period: svd1 or svd2
    :param is_pion: is pionID
    :param cut_value: the cut value for PID , eg 6 for id > 0.6
    :param eff_source: "data" for eff_dt or "mc" for eff_mc
    :param particle_names: tuple of particle name strings, default ("pip", "pim") 
        for branch pip_costheta , pip_p then define branch pip_pid_eff
    """            
    global CPP_FUNC_DEFINED

    if eff_source not in ("data", "mc"):
        raise ValueError("eff_source must be either 'data' or 'mc'")

    eff_column = "eff_dt" if eff_source == "data" else "eff_mc"

    def load_pid_eff(filename):
        df = pd.read_csv(filename, sep=r'\s+', comment='#',
                names=['kind', 'pid', 'map', 'eff_dt', 'statErr_dt', 'systErr_dt', 
                       'eff_mc', 'statErr_mc', 'ratio', 'statErr_ratio', 'systErr_ratio' ,'flag' ]) 
        if is_pion:
            df = df[(df['kind']==2) & (df['pid']== cut_value) & (df['flag']==0)]
        else:
            df = df[(df['kind']==0) & (df['pid']== cut_value) & (df['flag']==0)]
        df['map'] = df['map'].astype(str).str.zfill(4)
        df['pbin'] = df['map'].str[:2].astype(int) - 1
        df['ctbin'] = df['map'].str[2:].astype(int) - 1

        return df

    dir = "/gpfs/home/belle2/wangz/Work/SpinAlignment/offline/draw/belle_kid_eff/"
    if period == "svd1":
        df = load_pid_eff(dir + 'kideff-2006-svd1-all.dat')
    elif period == "svd2":
        df = load_pid_eff(dir + 'kideff-2010.dat')
    else:
        return 0

    if not CPP_FUNC_DEFINED:
        R.gInterpreter.Declare("""
            #include <map>
            #include <tuple>
            #include <vector>
            #include <algorithm>
            #include "ROOT/RVec.hxx"

            std::map<std::pair<int, int>, double> eff_map;

            void fill_eff_map(int pbin, int ctbin, double eff) {
                eff_map[std::make_pair(pbin, ctbin)] = eff;
            }
                               
            void clear_eff_map() {
                eff_map.clear();
            }

            std::vector<double> plab_boundaries = { 0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 4.0, 4.5, 100.0 };

            std::vector<double> costheta_boundaries = { -1, -0.612, -0.511, -0.300, -0.152, 0.017,
                0.209, 0.355, 0.435, 0.542, 0.692, 0.842, 1 };

            double get_eff(double plab, double costheta) {
                int pbin = std::upper_bound(plab_boundaries.begin(), plab_boundaries.end(), plab) - plab_boundaries.begin() - 1;
                int ctbin = std::upper_bound(costheta_boundaries.begin(), costheta_boundaries.end(), costheta) - costheta_boundaries.begin() - 1;
                auto key = std::make_pair(pbin, ctbin);
                if (eff_map.find(key) != eff_map.end()) {
                    return eff_map[key];
                } else {
                    return 1.0;
                }
            }

            ROOT::VecOps::RVec<double> get_eff(const ROOT::VecOps::RVec<double>& plab, const ROOT::VecOps::RVec<double>& costheta) {
                size_t n = std::min(plab.size(), costheta.size());
                ROOT::VecOps::RVec<double> out(n);
                for (size_t i = 0; i < n; ++i) {
                out[i] = get_eff(plab[i], costheta[i]);
                }
                return out;
            }
                                   
        """)
        CPP_FUNC_DEFINED = True

    R.clear_eff_map()
    for _, row in df.iterrows():
        R.fill_eff_map(int(row['pbin']), int(row['ctbin']), float(row[eff_column]))
    
    rdf = rdf.Define(f"{particle_names[0]}_pid_eff", f"get_eff({particle_names[0]}_p, {particle_names[0]}_costheta)")
    rdf = rdf.Define(f"{particle_names[1]}_pid_eff", f"get_eff({particle_names[1]}_p, {particle_names[1]}_costheta)")

    return rdf




