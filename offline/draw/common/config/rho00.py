"""Configs for rho00 extraction (binning, defaults, labels)."""

# Default branch/column names (override per call when using other particles)
DEFAULT_Z_BRANCH = "Ks_z"
DEFAULT_HELICITY_BRANCH = "Ks_helicity_angle"
DEFAULT_WEIGHT_BRANCH = "Ks_weight"
DEFAULT_BRANCH_MAP = {"z": DEFAULT_Z_BRANCH, "helicity_angle": DEFAULT_HELICITY_BRANCH}

# Default CSV column names for nsig table
DEFAULT_CSV_COL_MAP = {
    "z": "Ks_z_center",
    "helicity_angle": "Ks_helicity_angle_center",
    "nsig": "nsig",
    "nsig_err": "nsig_err_hi",
}

