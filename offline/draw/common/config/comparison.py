"""Default variable and binning configuration for data/MC comparisons."""

# Default per-variable binning and axis labels (merged Ks and phi variables)
DEFAULT_VAR_CONFIG = {
    "phi_M": (60, 1.0, 1.06, ";M_{#phi};[MeV]"),
    "phi_p": (100, 0.0, 5.0, ";p(#phi);[GeV/c]"),
    "phi_costheta": (20, -1.0, 1.0, ";cos#theta(#phi);[]"),
    "phi_phi": (20, -3.14, 3.14, ";#varphi(#phi);[rad]"),
    "kp_p": (90, 0.0, 4.5, ";p(K^{+});[MeV]"),
    "km_p": (90, 0.0, 4.5, ";p(K^{-});[MeV]"),
    "kp_costheta": (20, -1.0, 1.0, ";cos#theta(K^{+});[]"),
    "km_costheta": (20, -1.0, 1.0, ";cos#theta(K^{-});[]"),
    "kp_phi": (20, -3.14, 3.14, ";#phi(K^{+});[rad]"),
    "km_phi": (20, -3.14, 3.14, ";#phi(K^{-});[rad]"),
    "Ks_M": (100, 0.47, 0.52, ";M_{K_{s}};[MeV]"),
    "Ks_p": (100, 0, 5, ";p(K_{s});[GeV/c]"),
    "Ks_costheta": (20, -1, 1, ";cos#theta(K_{s});[]"),
    "Ks_phi": (20, -3.14, 3.14, ";#varphi(K_{s});[rad]"),
    "pip_p": (90, 0, 4.5, ";p(#pi^{+});[MeV]"),
    "pim_p": (90, 0, 4.5, ";p(#pi^{-});[MeV]"),
    "pip_costheta": (20, -1, 1, ";cos#theta(#pi^{+});[]"),
    "pim_costheta": (20, -1, 1, ";cos#theta(#pi^{-});[]"),
    "pip_phi": (20, -3.14, 3.14, ";#phi(#pi^{+});[rad]"),
    "pim_phi": (20, -3.14, 3.14, ";#phi(#pi^{-});[rad]"),
}

# Example MC weight config: {"Ks": "Ks_weight", "pip": "pip_pid_weight"}
DEFAULT_KS_WEIGHT_CONFIG = {"Ks": "Ks_weight", "pip": "pip_pid_weight", "pim": "pim_pid_weight"}
DEFAULT_PHI_WEIGHT_CONFIG = {"phi": "phi_weight", "kp": "kp_pid_weight", "km": "km_pid_weight"}

DEFAULT_PLOTS_PER_ROW = 3
DEFAULT_WDITH_PER_PLOT = 0.32
    
