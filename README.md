# MESA School 2026

This repository contains supplementary materials for the 2026 MESA Summer School.

Main school website: <https://mesastar.org/summer-school-2026/>

## Lab Materials

### `lab-1.md`

Supplementary reading notes for Lab 1. Written by Meng Sun Proofread by Dr. Rich Townsend.

### `lab-3.md`

Lab 3: **Beyond the Core: Echoes of Overshoot**

This tutorial was written by Meng Sun for Day 2, Lab 3 of the 2026 MESA School.

The lab focuses on convective-boundary mixing in MESA, including step overshooting, exponential overshooting, and penetration convection, and connects the resulting stellar structures to g-mode pulsations computed with GYRE.

## Solution Files and Working Directory

All example inlists, solution files, and figures for this lab are currently collected here.
These files are:

- `inlist_step_ov_ZAMS_solution`  
  (For step overshooting) Evolves the model from the pre-main sequence to the ZAMS and saves a ZAMS model for later use.

- `inlist_step_ov_MS_solution`  
   (For step overshooting) Starts from the saved ZAMS model and evolves the star to a later main-sequence phase. This run also writes `.GYRE` files for the asteroseismic analysis.

- `inlist_exp_ov_ZAMS_solution`  
  Same as the step-overshooting ZAMS run, but using exponential overshooting.

- `inlist_exp_ov_MS_solution`  
  Same as the step-overshooting main-sequence run, but using exponential overshooting.

- `inlist_penetration_ZAMS_solution`  
  Evolves the penetration-convection model to the ZAMS.

- `inlist_penetration_MS_solution`  
  Evolves the penetration-convection model from the saved ZAMS model to the later main-sequence phase.

- `run_star_extras_solution.f90`  
  The modified `run_star_extras.f90` file used for the convective-penetration prescription. This file should be copied into the `src/` directory of your MESA work directory as

  ```text
  src/run_star_extras.f90
  ```

  Then rebuild the work directory with

  ```bash
  ./clean
  ./mk
  ```
  
## Acknowledgements

Meng Sun from the National Astronomical Observatories, Chinese Academy of Sciences, thanks Daniel Lecoanet for designing this lab. Meng Sun also thanks Lynn Buchele, Caleb Eastlund, Ducheng Lu, Lucas de Sá and Mathijs Vanrespaille for helping test the lab materials, estimate the computational cost and providing useful comments and feedbacks.
