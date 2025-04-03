# Planar Waveguide Analysis in MATLAB

## Overview

This MATLAB script analyzes the dispersion characteristics of planar optical waveguides for both TE (Transverse Electric) and TM (Transverse Magnetic) modes. It calculates the effective refractive index (`n_eff`) as a function of core thickness (`h`) and determines key parameters such as single-mode propagation ranges, minimum substrate/cladding thicknesses, and propagation constants (`β`).

## Features

- **Waveguide Structures**: Analyzes two predefined waveguides with distinct refractive indices for core, cladding, and substrate for TE and TM modes.
- **Dispersion Curves**: Plots `n_eff` versus core thickness for modes `m = 0` to `m = 5`.
- **Single-Mode Range**: Identifies the thickness range for single-mode propagation.
- **Minimum Thickness Calculation**: Computes minimum substrate and cladding thicknesses based on the upper bound of the single-mode range.
- **Midpoint Analysis**: Calculates `n_eff` and `β` at the midpoint of the single-mode range.
- **Visualization**: Generates plots with labeled axes, legends, and grid for easy interpretation.

## Parameters

- **Wavelength**: `λ = 1550 nm` (fixed).
- **Core Thickness Range**: `h` varies from 0.1 µm to 10 µm with 500 points.
- **Waveguides**:
  - **Waveguide (a)**:
    - TE: `n_cl = 1.5120`, `n_c = 1.5375`, `n_s = 1.4446`
    - TM: `n_cl = 1.5110`, `n_c = 1.5370`, `n_s = 1.4446`
  - **Waveguide (b)**:
    - TE: `n_cl = 1.4446`, `n_c = 2.1381`, `n_s = 1.4446`
    - TM: `n_cl = 1.4446`, `n_c = 2.2213`, `n_s = 1.4446`

## Key Equations

- **TE Mode**:  
  `f_TE(n_eff, h, m, n_c, n_cl, n_s, λ)` solves the transcendental equation for TE modes.
- **TM Mode**:  
  `f_TM(n_eff, h, m, n_c, n_cl, n_s, λ)` solves the transcendental equation for TM modes.
- **Propagation Constants**:  
  - `p = k0 * sqrt(n_eff^2 - n_s^2)`  
  - `q = k0 * sqrt(n_eff^2 - n_cl^2)`  
  where `k0 = 2π/λ`.

## Usage

1. **Requirements**: MATLAB with the Optimization Toolbox (for `fzero`).
2. **Run**: Execute the script in MATLAB. It will:
   - Compute `n_eff` for each waveguide and mode.
   - Display single-mode ranges and thickness parameters in the command window.
   - Generate dispersion curve plots for each waveguide.
3. **Output**:
   - Console output with calculated values (e.g., `n_eff`, `β`, thickness ranges).
   - Figures showing dispersion curves for TE (solid lines) and TM (dashed lines) modes.

## Notes

- The script uses `fzero` to solve the mode equations numerically, with bounds set between the maximum of cladding/substrate indices and just below the core index.
- NaN values indicate modes that do not propagate at a given thickness.
- Units are converted to micrometers (µm) for readability in outputs and plots.

## License

This code is provided for educational purposes. Feel free to modify and adapt it as needed.
