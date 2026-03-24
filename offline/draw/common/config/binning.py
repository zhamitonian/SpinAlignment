"""Shared binning configurations for comparisons and rho00 extraction."""

import numpy as np

# Data/MC comparison binning (vector branches)
DEFAULT_KS_BINNING = {
    "Ks_z": np.linspace(0.0, 1.0, 21).tolist(),
    "Ks_helicity_angle": np.linspace(-1.0, 1.0, 11).tolist(),
}

DEFAULT_PHI_BINNING = {
    "phi_z": np.linspace(0.2, 1.0, 11).tolist(),
    "phi_thrust_pt": [0.0, 0.125, 0.25, 0.375, 0.5, 0.6611, 0.8688, 1.1366, 1.4817, 1.9265, 2.5],
    "phi_helicity_angle": np.linspace(-1.0, 1.0, 11).tolist(),
}

def _normalize_binning(binning_dict):
    """Normalize binning dict to ensure all values are lists."""
    normalized = {}
    for var, bins in binning_dict.items():
        if isinstance(bins, list):
            normalized[var] = bins
        elif isinstance(bins, np.ndarray):
            normalized[var] = bins.tolist()
        else:
            raise ValueError(f"Binning for variable '{var}' must be a list or numpy array.")
    return normalized