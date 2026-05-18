# Lab 3: Convective Boundary Mixing, Stellar Structure, and g Modes

## Task 0. Goal of This Lab

In this lab, we will study how different convective boundary mixing prescriptions affect stellar evolution, internal structure, and g-mode pulsations. There are three treatments of mixing near the top boundary of a hydrogen burning convective core:

1. step overshooting,
2. exponential overshooting,
3. convective penetration.

The main goal is to understand how these mixing prescriptions modify the near-core chemical-gradient region and the Brunt–Väisälä frequency profile. These structural differences may leave measurable signatures in stellar eigenmodes.

In the first part of the lab, we will build MESA models using different mixing prescriptions. Next, we will inspect their internal structures at an intermediate main-sequence stage. Finally, we will use GYRE to compute g-mode frequencies, compare them with a reference set of modes, and identify the best-fit model.

---

## Task 1. Step and Exponential Overshooting

In MESA, step and exponential overshooting are built-in prescriptions that can be controlled from the inlist. For both prescriptions, we apply overshooting at the top boundary of the convective core:

```fortran
overshoot_zone_type(1) = 'any'
overshoot_zone_loc(1)  = 'core'
overshoot_bdy_loc(1)   = 'top'
```

These lines tell MESA where the overshooting is applied:

- `any`: allow this prescription to be applied to any relevant convective boundary;
- `core`: apply it to a convective core;
- `top`: apply it at the outer edge of the convective core.

---

## Task 2. Step Overshooting

Step overshooting assumes that the material is fully mixed out to a fixed distance beyond the formal convective boundary. Use the same location controls, but change the scheme:

```fortran
overshoot_scheme(1) = 'step'

overshoot_f(1) = 0.2d0
overshoot_f0(1) = 0.005d0
overshoot_D_min = 1d-2
```

The difference between `overshoot_f` and `overshoot_f0` is shown schematically below. `overshoot_D_min` sets the lower cutoff for the overshoot mixing diffusion coefficient.

![MESA overshooting schematic](https://mesa-leuven.4d-star.org/tutorials/monday/overshoot_explanation.png)

*Credit: 2025 MESA School in Leuven Day 1 tutorial material.*

The convective boundary is where the convective diffusion coefficient drops to zero. MESA steps slightly inward from this boundary by a distance `overshoot_f0 * H_p`. The main overshooting length scale is controlled by `overshoot_f * H_p`, where `H_p` is the local pressure scale height.

---

## Task 3. Exponential Overshooting

Exponential overshooting assumes that the mixing coefficient decreases smoothly outside the convective boundary.

A typical setup is

```fortran
overshoot_scheme(1) = 'exponential' ! options: 'exponential', 'step', 'other'

overshoot_f(1) = 0.02d0
overshoot_f0(1) = 0.005d0
overshoot_D_min = 1d-2
```

In the model grid, we will vary `overshoot_f(1)`.

---

## Task 4. Convective Penetration

Convective penetration is different from standard MESA overshooting. Material beyond the convective boundary is chemically mixed, but the thermal structure is usually still treated as radiative. In convective penetration, convective motions penetrate into the formally stable region and can modify both the chemical composition and the thermal stratification. In the implementation used here, the penetration extent is computed inside `run_star_extras.f90`.

For the convective penetration runs, use

```fortran
! Overshooting
overshoot_scheme(1) = 'other'

overshoot_f(1) = 0.00
overshoot_f0(1) = 0.005d0
overshoot_D_min = 1d-2
```

The key line is

```fortran
overshoot_scheme(1) = 'other'
```

This tells MESA to call the user-supplied overshooting routine from `run_star_extras.f90`. You will be given a clean MESA `run_star_extras.f90` file and a modified version that implements convective penetration.

Your task is to identify which parts of `run_star_extras.f90` are needed for the custom penetration scheme. The key pieces are listed below.

---

### Task 4.1 Define Extra Variables

Near the top of the module, after

```fortran
implicit none
```

the modified file defines extra variables that store information about the convective core and the penetration zone. For example, the implementation tracks quantities such as

```fortran
m_core ! the convective core mass
mass_PZ ! the mass of the penetration zone
delta_r_PZ ! the radial width of the penetration zone
alpha_PZ ! the dimensionless penetration extent
r_core ! the radius of the convective core boundary
rho_core_top ! the density at the top of the core
```

---

### Task 4.2 Connect MESA to the Custom Overshooting Routine

Inside `extras_controls`, MESA must be told which custom routine to call when the inlist says

```fortran
overshoot_scheme(1) = 'other'
```

The important line has the form

```fortran
s% other_overshooting_scheme => extended_convective_penetration
```

This is the hook that connects the inlist setting to the custom convective penetration routine.

---

### Task 4.3 Add Extra History Columns

This is done by modifying two routines:

```fortran
how_many_extra_history_columns
data_for_extra_history_columns
```

The modified implementation writes seven extra history columns:

```fortran
m_core
mass_pen_zone
delta_r_pen_zone
alpha_pen_zone
r_core
rho_core_top_pen
r_cb
```

---

### Task 4.4 Add the Custom Overshooting Routine

The main custom overshooting routine is called

```fortran
extended_convective_penetration
```

This routine does:

1. checks that the boundary is the top of a convective core;
2. calls another routine to compute the penetration-zone width;
3. uses the computed `alpha_PZ` as the width of a step like penetration region;
4. optionally attaches an exponential tail controlled by `overshoot_f(1)`;
5. returns the diffusion coefficient profile `D`.

A key line in this routine is

```fortran
call dissipation_balanced_penetration(s, id)
```

This computes the penetration zone extent.

Another important line is

```fortran
alpha_PZ = alpha_PZ + s%overshoot_f0(j)
```

This means that the final step like penetration region includes the computed penetration width plus the small offset set by `overshoot_f0`.

This is why we use

```fortran
overshoot_f(1) = 0.00
overshoot_f0(1) = 0.005
```

---

### Task 4.5 Compute the Penetration Width

The penetration width is computed in the routine

```fortran
dissipation_balanced_penetration
```

This routine estimates how far the convective penetration zone should extend beyond the convective boundary. For this lab, let's focus on identifying how the code computes

```fortran
delta_r_PZ
alpha_PZ
```

The key relation is

```fortran
alpha_PZ = delta_r_PZ / h
```

where `h` is the local pressure scale height near the convective core boundary.

---

### Task 4.6 Optional: Extra Mesh Refinement

The modified implementation also includes an optional mesh refinement routine near the core boundary. This is useful because the Brunt–Väisälä frequency and the composition gradient can vary rapidly near the convective boundary. The relevant hook has the form

```fortran
s% use_other_mesh_delta_coeff_factor = .true.
s% other_mesh_delta_coeff_factor => mesh_delta_coeff_core_boundary
```

---

## Task 5. Model Grid

Run the model grid listed in the shared spreadsheet:

[Lab 3 grid tracker](https://docs.google.com/spreadsheets/d/1v9Dq4AV1ZGssSdy1lQE3uiXW0afyK1mRk9uvBgGOaGI/edit?usp=sharing)

The grid spans

```text
Initial mass: 3.0 to 8.0 Msun, step 0.5 Msun
```

For each model, evolve from ZAMS to TAMS.

For this lab, define TAMS as

```fortran
xa_central_lower_limit_species(1) = 'h1'
xa_central_lower_limit(1) = 0.01
```
## Solution Files and Naming Conventions

Example solution files are provided in the same GitHub directory as this tutorial. The filenames contain placeholders such as `X.X`. Please replace these placeholders with your desired mixing parameters and initial stellar mass before running the models.

For the penetration-convection runs, remember that the main penetration strength parameter is coded in `run_star_extras_solution.f90`. You should change

```fortran
real(dp), parameter :: f = X.Xd0
```
near line 536 to the desired value, for example f = 0.98d0, 0.86d0, or 0.72d0. After changing this value, recompile with:

```bash
./mk
```
In the solution files, we use separate local output directories for the three mixing prescriptions:

```text
LOGS_step_ov
LOGS_exp_ov
LOGS_PC
```

For example, in the exponential overshoot ZAMS run, the saved model may be written as 

```fortran
save_model_filename = './LOGS_exp_ov/exp_ov_zams.model'
```
When you run a different parameter value, you may want to change the output directory or saved model filename to avoid overwriting previous runs.

## Example One Run
For step overshooting, two solution inlists are provided:

inlist_step_ov_ZAMS_solution
inlist_step_ov_MS_solution

Use `inlist_step_ov_ZAMS_solution` for the first-stage run, from the pre-main sequence to ZAMS. This run uses:

```fortran
stop_near_zams = .true.
```

and saves the ZAMS model. Then use `inlist_step_ov_MS_solution` for the second-stage run, from the saved ZAMS model to the late main sequence. This run loads the saved ZAMS file and stops when the central hydrogen abundance reaches 0.1:

```fortran
xa_central_lower_limit_species(1) = 'h1'
xa_central_lower_limit(1) = 0.1
```
The same two stage workflow should be followed for the exponential overshoot and penetration convection cases.

---

## 7. What to Record

For this lab, you need to record the seismic fit quality for each model. The shared Google Sheet already provides the target g-mode frequencies for `n_pg = -20` to `-10`. For each MESA+GYRE model, use the final MESA model, namely the profile with central hydrogen abundance closest to `Xc(H) = 0.1`, extract the corresponding GYRE model frequencies, and compute a single `Chi^2` value. An unweighted Chi^2 can be computed as: 

```Python
Chi2 = np.sum((freq_model - freq_target)**2)
```

Record this `Chi^2` value in the table cell corresponding to the model's initial mass and mixing parameter. After all models are filled in, the cell with the smallest `Chi^2` identifies the best-fit model within this grid.



