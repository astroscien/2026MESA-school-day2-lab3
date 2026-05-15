---
title: Lab 3
nav_order: 3
---

# Lab 3: Convective Boundary Mixing in MESA

In this lab, we compare three ways of treating mixing near the edge of a convective core:

1. step overshooting,
2. exponential overshooting,
3. convective penetration.

The goal is to understand how different convective boundary mixing prescriptions modify the stellar structure near the core boundary, and how those structural changes affect g-mode pulsations.

We will first compare the MESA namelist settings for step and exponential overshooting. Then we will look at how a custom `run_star_extras.f90` file can be used to implement convective penetration.

---

## 1. Step and Exponential Overshooting

In MESA, step and exponential overshooting can be controlled directly from the inlist.

For both prescriptions, we will apply overshooting at the top boundary of the hydrogen-burning convective core:

```fortran
overshoot_zone_type(1) = 'burn_H'
overshoot_zone_loc(1)  = 'core'
overshoot_bdy_loc(1)   = 'top'
```

These three lines tell MESA where the overshooting is applied:

- `burn_H`: apply this prescription during core hydrogen burning;
- `core`: apply it to a convective core;
- `top`: apply it at the outer edge of the convective core.

---

## 2. Exponential Overshooting

Exponential overshooting assumes that the mixing coefficient decreases smoothly outside the formal convective boundary.

A typical setup is

```fortran
overshoot_scheme(1) = 'exponential' ! options: 'exponential', 'step', 'other'
overshoot_zone_type(1) = 'burn_H'   ! options: 'burn_H', 'burn_He', 'burn_Z', 'nonburn', 'any'
overshoot_zone_loc(1) = 'core'      ! options: 'core', 'shell', 'any'
overshoot_bdy_loc(1) = 'top'        ! options: 'bottom', 'top', 'any'

overshoot_f(1) = 1.6d-2
overshoot_f0(1) = 8.0d-3
overshoot_D_min = 1d-2
```

Here, `overshoot_f(1)` controls the scale length of the exponential decay of the mixing coefficient outside the convective boundary.

The parameter `overshoot_f0(1)` tells MESA where to evaluate the reference diffusion coefficient slightly inside the convective region.

In the model grid, we will vary `overshoot_f(1)`.

---

## 3. Step Overshooting

Step overshooting assumes that the material is fully mixed out to a fixed distance beyond the formal convective boundary.

Use the same location controls, but change the scheme:

```fortran
overshoot_scheme(1) = 'step'
overshoot_zone_type(1) = 'burn_H'
overshoot_zone_loc(1) = 'core'
overshoot_bdy_loc(1) = 'top'

overshoot_f(1) = 0.2d0
overshoot_f0(1) = 0.1d0
overshoot_D_min = 1d-2
```

For step overshooting, `overshoot_f(1)` gives the approximate radial extent of the fully mixed overshoot region in units of the local pressure scale height.

In this lab, this parameter is equivalent to the commonly used `alpha_ov`.

---

## 4. Convective Penetration

Convective penetration is different from standard MESA overshooting.

In ordinary overshooting, material beyond the convective boundary is chemically mixed, but the thermal structure is usually still treated as radiative.

In convective penetration, convective motions penetrate into the formally stable region and can modify both the chemical composition and the thermal stratification. In the implementation used here, the penetration extent is not entered directly as a standard MESA namelist parameter. Instead, it is computed inside `run_star_extras.f90`.

For the convective penetration runs, use

```fortran
! Overshooting
overshoot_scheme(1) = 'other'
overshoot_zone_type(1) = 'any'
overshoot_zone_loc(1) = 'core'
overshoot_bdy_loc(1) = 'top'

overshoot_f0(1) = 0.005
overshoot_f(1) = 0.00
! overshoot_D0(1) = 0.005
```

The key line is

```fortran
overshoot_scheme(1) = 'other'
```

This tells MESA not to use one of its built-in overshooting prescriptions. Instead, MESA will call a user-supplied overshooting routine from `run_star_extras.f90`.

Do not treat `alpha_pen` as an input parameter. In this implementation, the penetration extent is computed by the code and written to the history output as

```fortran
alpha_pen_zone
```

---

## 5. What Needs to Be Modified in `run_star_extras.f90`

You will be given a clean MESA `run_star_extras.f90` file and a modified version that implements convective penetration.

Your task is not to blindly copy the full solution. Instead, identify which parts of `run_star_extras.f90` are needed for the custom penetration scheme and understand what each part does.

The key pieces are listed below.

---

### 5.1 Define Extra Variables

Near the top of the module, after

```fortran
implicit none
```

the modified file defines extra variables that store information about the convective core and the penetration zone.

For example, the implementation tracks quantities such as

```fortran
m_core
mass_PZ
delta_r_PZ
alpha_PZ
r_core
rho_core_top
```

These variables are used to store

- the convective core mass,
- the mass of the penetration zone,
- the radial width of the penetration zone,
- the dimensionless penetration extent,
- the radius of the convective core boundary,
- the density at the top of the core.

Some of these should be declared as real variables. Others, such as mesh indices, should be declared as integers. You should check the modified file and decide which type each variable needs.

---

### 5.2 Connect MESA to the Custom Overshooting Routine

Inside `extras_controls`, MESA must be told which custom routine to call when the inlist says

```fortran
overshoot_scheme(1) = 'other'
```

The important line has the form

```fortran
s% other_overshooting_scheme => extended_convective_penetration
```

This is the hook that connects the inlist setting to the custom convective penetration routine.

Without this line, MESA will not know which user-defined overshooting scheme to use.

---

### 5.3 Add Extra History Columns

Because `alpha_pen` is not an input parameter, we need to record the computed value in `history.data`.

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

The most important quantity for this lab is

```fortran
alpha_pen_zone
```

This is the computed penetration extent in units of the local pressure scale height.

In the code, this corresponds to the internal variable

```fortran
alpha_PZ
```

---

### 5.4 Add the Custom Overshooting Routine

The main custom overshooting routine is called

```fortran
extended_convective_penetration
```

This routine does several things:

1. checks that the boundary is the top of a convective core;
2. calls another routine to compute the penetration-zone width;
3. uses the computed `alpha_PZ` as the width of a step-like penetration region;
4. optionally attaches an exponential tail controlled by `overshoot_f(1)`;
5. returns the diffusion coefficient profile `D`.

A key line in this routine is

```fortran
call dissipation_balanced_penetration(s, id)
```

This computes the penetration-zone extent.

Another important line is

```fortran
alpha_PZ = alpha_PZ + s%overshoot_f0(j)
```

This means that the final step-like penetration region includes the computed penetration width plus the small offset set by `overshoot_f0`.

This is why, for the first pass, we use

```fortran
overshoot_f0(1) = 0.005
overshoot_f(1) = 0.00
```

rather than trying to scan `alpha_pen` directly.

---

### 5.5 Compute the Penetration Width

The penetration width is computed in the routine

```fortran
dissipation_balanced_penetration
```

This routine estimates how far the convective penetration zone should extend beyond the formal convective boundary.

The basic idea is that convection produces buoyant work inside the convective core. The penetration zone extends outward until this work is balanced by dissipation and negative buoyant work in the stable region.

For this lab, you do not need to rederive the prescription. Instead, focus on identifying how the code computes

```fortran
delta_r_PZ
alpha_PZ
```

The key relation is

```fortran
alpha_PZ = delta_r_PZ / h
```

where `h` is the local pressure scale height near the convective core boundary.

Thus, `alpha_PZ` is the penetration-zone width measured in units of the local pressure scale height.

---

### 5.6 Optional: Extra Mesh Refinement

The modified implementation also includes an optional mesh refinement routine near the core boundary.

This is useful because the Brunt-Vaisala frequency and the composition gradient can vary rapidly near the convective boundary.

The relevant hook has the form

```fortran
s% use_other_mesh_delta_coeff_factor = .true.
s% other_mesh_delta_coeff_factor => mesh_delta_coeff_core_boundary
```

This part is useful for obtaining cleaner profiles, but it is secondary to the main convective penetration implementation.

---

## 6. Model Grid

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

---

## 7. What to Record

For each run, record the following quantities in the spreadsheet:

- initial mass,
- mixing prescription,
- overshooting or penetration control parameters,
- wall-clock runtime from ZAMS to TAMS,
- final age,
- final luminosity,
- final radius,
- for penetration runs, `alpha_pen_zone`,
- the first 10 to 20 g modes,
- the `g_10` mode as a compact seismic diagnostic.

---

## 8. Discussion Questions

1. How does the choice of convective boundary mixing prescription change the TAMS structure?
2. Does step overshooting produce a sharper or smoother core-boundary structure than exponential overshooting?
3. In the penetration runs, how does the computed `alpha_pen_zone` vary with stellar mass?
4. Do the differences in the Brunt-Vaisala frequency profile show up in the g-mode spectrum?
5. Can a tuned overshooting model mimic the seismic signature of convective penetration?
