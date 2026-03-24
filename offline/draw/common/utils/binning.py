"""Helpers for interpreting binning configs.

Supports both uniform tuple-style (nbins, min, max) and explicit edge lists.
"""

def normalize_bin_cfg(cfg, name):
    """Return (nbins, xmin, xmax) from either (nbins, min, max) or an edge list.

    Args:
        cfg: tuple (nbins, min, max) or list/tuple of edges.
        name: variable name for error messages.
    """
    if isinstance(cfg, tuple) and len(cfg) == 3:
        nbins, xmin, xmax = cfg
        return int(nbins), float(xmin), float(xmax)

    if isinstance(cfg, (list, tuple)):
        boundaries = list(cfg)
        if len(boundaries) < 2:
            raise ValueError(f"Invalid binning for {name}: need at least two boundaries, got {len(boundaries)}")
        return len(boundaries) - 1, float(boundaries[0]), float(boundaries[-1])

    raise ValueError(f"Invalid binning format for {name}: {cfg}")
