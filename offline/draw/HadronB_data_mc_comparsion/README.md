# Hadron B Data-MC Comparison

This directory contains scripts for data-MC comparison after Belle's hadron B selection.

## Overview

Due to the large event size, we pre-process events for each experiment by saving coordinate histograms to ROOT files.

## Main Scripts

### `by_exp.py`
Compares qqbar's MC with data for each experiment individually.

**Purpose:**
- Handle primary vertex Z shift that varies by experiment number
- Use sideband regions to estimate beam wall effects
- Generate ROOT files containing multiple histograms

**Workflow:**
1. Extract the center value of `PrimeVz` by fitting with a Gaussian (`extract_center_primeVz.py`)
2. Correct for the experiment-dependent vertex shift
3. Save processed histograms to ROOT files

### `submit_comparing_by_exp.py`
Submits batch jobs for experiment-by-experiment comparison to LSF.

### `plot_by_combining_hists.py`
Performs the final combined comparison:
- Includes QED MC contributions
- Incorporates beam wall effects
- Combines all MC types across all experiments
- Generates comparison plots

### `qed_mc_distribution_by_type.py`
*(Optional)* Demonstrates QED MC distributions by type for validation purposes.

## Processing Flow

```
Raw Data/MC → by_exp.py (per experiment) → ROOT histograms
                    ↓
          extract_center_primeVz.py (fit PrimeVz)
                    ↓
          run by_exp.py again (with Vertex shift correction)
                    ↓
          plot_by_combining_hists.py (combine all)
                    ↓
          Final comparison plots
```

## Notes

- **PrimeVz shift:** Different experiments have different beam vertex positions, requiring experiment-specific corrections
- **Beam wall background:** Estimated using data sideband regions rather than MC simulation