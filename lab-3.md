# Lab 3: Convective Boundary Mixing, Stellar Structure, and g Modes

## Overview and Starting Files

In this lab, we will study how different convective boundary mixing prescriptions affect stellar evolution, internal structure, and g-mode pulsations. There are three treatments of mixing near the top boundary of a hydrogen burning convective core:

1. step overshooting,
2. exponential overshooting,
3. convective penetration.

The main goal is to understand how these mixing prescriptions modify the near-core chemical-gradient region and the Brunt–Väisälä frequency profile. These structural differences may leave measurable signatures in stellar eigenmodes.

We use a two-step evolution. Taking the step overshooting case as an example:

1. `inlist_step_ov_ZAMS_solution` evolves the model from the pre-main sequence to the ZAMS and saves the ZAMS model for the next step.
2. `inlist_step_ov_MS_solution` starts from the saved ZAMS model and evolves the star to a later main-sequence phase. During this step, MESA also outputs `.GYRE` files for the asteroseismic analysis.

Before running the models, remember to replace all placeholder values such as

```fortran
X.X
```

with the parameter values assigned in the spreadsheet.

For the current grid, use:

```text
Step overshooting:
0.02, 0.14, 0.28

Exponential overshooting:
0.002, 0.014, 0.028
```

In the first part of the lab, we will build MESA models using different mixing prescriptions. Next, we will inspect their internal structures at an intermediate main-sequence stage. Finally, we will use GYRE to compute g-mode frequencies, compare them with a reference set of modes, and identify the best-fit model.

---

## Background: Standard Overshooting Prescriptions

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

## Step Overshooting

Step overshooting assumes that the material is fully mixed out to a fixed distance beyond the formal convective boundary. Use the same location controls, but change the scheme:

```fortran
overshoot_scheme(1) = 'step'

overshoot_f(1) = X.Xd0
overshoot_f0(1) = 0.005d0
overshoot_D_min = 1d-2
```

The difference between `overshoot_f` and `overshoot_f0` is shown schematically below. `overshoot_D_min` sets the lower cutoff for the overshoot mixing diffusion coefficient.

<img src="https://mesa-leuven.4d-star.org/tutorials/monday/overshoot_explanation.png" width="650">

*Credit: 2025 MESA School in Leuven Day 1 tutorial material (by Dr. Poojan Agrawal).*

The convective boundary is where the convective diffusion coefficient drops to zero. MESA steps slightly inward from this boundary by a distance `overshoot_f0 * H_p`. The main overshooting length scale is controlled by `overshoot_f * H_p`, where `H_p` is the local pressure scale height.

---

## Exponential Overshooting

Exponential overshooting assumes that the mixing coefficient decreases smoothly outside the convective boundary.

A typical setup is

```fortran
overshoot_scheme(1) = 'exponential' ! options: 'exponential', 'step', 'other'

overshoot_f(1) = X.XXd0
overshoot_f0(1) = 0.005d0
overshoot_D_min = 1d-2
```

In the model grid, we will vary `overshoot_f(1)`.

---

## Convective Penetration

This section is a guided code-reading exercise. The goal is to understand how the custom scheme is connected to MESA.

Convective penetration is different from standard MESA overshooting. Material beyond the convective boundary is chemically mixed, but the thermal structure is usually still treated as radiative. In convective penetration, convective motions penetrate into the formally stable region and can modify both the chemical composition and the thermal stratification. In the implementation used here, the penetration extent is computed inside `run_star_extras.f90`. The coding part of this implementation is relatively complicated. For this lab, the task is to identify which parts of `run_star_extras.f90` are needed for the custom penetration scheme, understand what each part does, and then use the supplied solution file as the working implementation.

In the inlists for the convective penetration runs, use

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

Before running the models, find the line in `run_star_extras.f90`

```fortran
real(dp), parameter :: f = X.Xd0
```

and replace `X.Xd0` with the value specified for your run (0.98, 0.86 or 0.72).

---

### Guided Check 1: Define Extra Variables

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


### Guided Check 2: Connect MESA to the Custom Overshooting Routine

Inside `extras_controls`, MESA must be told which custom routine to call when the inlist says

```fortran
overshoot_scheme(1) = 'other'
```

The important line has the form

```fortran
s% other_overshooting_scheme => extended_convective_penetration
```

This is the hook that connects the inlist setting to the custom convective penetration routine. Therefore, the same `run_star_extras.f90` file can be used for the whole lab: the custom routine is only used when `overshoot_scheme(1) = 'other'`, while the step and exponential runs use MESA's built-in schemes.


### Guided Check 3: Add Extra History Columns

This is done by modifying two routines:

```fortran
how_many_extra_history_columns
data_for_extra_history_columns
```

The modified implementation writes seven extra history columns:

```fortran
m_core ! convective core mass
mass_pen_zone ! mass contained in the penetration zone
delta_r_pen_zone ! radial width of the penetration zone
alpha_pen_zone ! penetration-zone width in units of the local pressure scale height
r_core ! radius of the convective core boundary
rho_core_top_pen ! density near the top of the convective penetration zone
r_cb ! radius of the convective boundary
```


### Guided Check 4: Add the Custom Overshooting Routine

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


### Guided Check 5: Compute the Penetration Width

The penetration width is computed in the routine

```fortran
dissipation_balanced_penetration
```

This routine estimates how far the convective penetration zone should extend beyond the convective boundary. For this lab, let's focus on identifying how the code computes

```fortran
delta_r_PZ ! physical radial width of the penetration zone
alpha_PZ ! delta_r_PZ divided by the local pressure scale height
```

The key relation is

```fortran
alpha_PZ = delta_r_PZ / h
```

where `h` is the local pressure scale height near the convective core boundary.

### Guided Check 6: Optional: Extra Mesh Refinement

The modified implementation also includes an optional mesh refinement routine near the core boundary. This is useful because the Brunt–Väisälä frequency and the composition gradient can vary rapidly near the convective boundary. The relevant hook has the form

```fortran
s% use_other_mesh_delta_coeff_factor = .true.
s% other_mesh_delta_coeff_factor => mesh_delta_coeff_core_boundary
```

---

## Task 0. Model Grid

Each student should run only the assigned subset of models and record the result in the shared spreadsheet. A good default is to choose one initial mass, one mixing prescription, and one parameter value. If time allows, students with faster computers can help fill missing entries.

The suggested parameter space to explore is listed in the shared spreadsheet:
[Lab 3 grid tracker](https://docs.google.com/spreadsheets/d/1v9Dq4AV1ZGssSdy1lQE3uiXW0afyK1mRk9uvBgGOaGI/edit?usp=sharing)

For each model, evolve from ZAMS to an evolved main-sequence model with centra hydrogen mass fraction of 0.1.

## Solution Files and Naming Conventions

The solution files are available here: [inlists and run_star_extras solutions](https://github.com/astroscien/2026MESA-school-day2-lab3/tree/main).

Files contain placeholders such as `X.X`. Please replace these placeholders with your desired mixing parameters and initial stellar mass before running the models.

For the penetration-convection runs, remember that the main penetration strength parameter is coded in `run_star_extras_solution.f90`. You should change

```fortran
real(dp), parameter :: f = X.Xd0
```
near line 536 to the desired value, then recompile with:

```bash
./clean
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

## Task 1. What to Record

For this lab, you need to record the seismic fit quality for each model. The shared Google Sheet already provides the target g-mode frequencies for `n_pg = -20` to `-10`. For each MESA+GYRE model, use the final MESA model, namely the profile with central hydrogen abundance closest to `Xc(H) = 0.1`.

### Where to Find the GYRE Eigenmodes

GYRE reports the eigenmodes in two places:

1. in the terminal output while GYRE is running;
2. in the GYRE summary/output files, if summary output is enabled in the GYRE inlist.

For a quick check, the terminal output is often enough. During the root-solving step, GYRE prints a table like this:

```text
Root Solving
   l    m    n_pg    n_p    n_g       Re(omega)       Im(omega)        chi n_iter
   1    0     -20      0     20  0.63087941E+00  0.00000000E+00 0.3165E-12      7
   1    0     -19      0     19  0.66841184E+00  0.00000000E+00 0.3904E-12      8
   1    0     -18      0     18  0.69832663E+00  0.00000000E+00 0.1748E-12      6
   ...
   1    0     -10      0     10  0.12513287E+01  0.00000000E+00 0.1019E-12      7
```

For this lab, use the dipole modes with

```text
l = 1
m = 0
n_pg = -20 to -10
```

Equivalently, for these g modes, you can identify them by

```text
n_g = 20 to 10
```

because `n_pg` is negative for g modes.

### Important: Frequency Units

Be careful with units. By default, GYRE prints the dimensionless angular frequency `omega`.

The terminal column

```text
Re(omega)
```

is the real part of the mode frequency. In this lab, we are running adiabatic oscillation calculations. This means that we ignore heat exchange during one oscillation cycle, so the modes do not grow or decay. Therefore, the eigenfrequencies are real numbers. In a non-adiabatic calculation, the eigenfrequency is generally complex. The real part gives the oscillation frequency, while the imaginary part describes mode growth or damping. A growing mode is unstable and may become observable if it reaches a detectable amplitude, while a damped mode is stable and tends to decay unless it is continuously excited.

### Compute the Fit Quality

Once the model and target frequencies are in the same units, compute an unweighted `Chi^2` value:

```Python
Chi2 = np.sum((freq_model - freq_target)**2)
```

Record this `Chi^2` value in the table cell corresponding to the model's initial mass and mixing parameter. After all models are filled in, the cell with the smallest `Chi^2` identifies the best-fit model within this grid.

---

## Task 2. Compare the Structures at Xc(H) = 0.5

In this task, we compare the internal structures of the three mixing prescriptions at the same evolutionary stage. For each prescription, we will use the model whose central hydrogen abundance is closest to 0.5. The three MESA runs are stored in separate output directories: 

```text
LOGS_step_ov
LOGS_exp_ov
LOGS_PC
```

Run the plotting script from the directory that contains these three folders. The workflow is:

For each LOGS_* directory, find the main sequence history file, such as `step_ov_MS.history`, `exp_ov_MS.history`, or `PC_MS.history`. In that history file, find the `model_number` with `center_h1` is closest to 0.5. Use `profiles.index` to map that model_number to the corresponding profile_number. Read the corresponding `profile*.data` file. Plot the abundance profiles and propagation diagrams for the three mixing prescriptions on the same figure.


### Code Block 1: Basic Setup

This block imports the required packages, defines the target central hydrogen abundance, and lists the three output directories.

<details>
    <summary>Python: imports, constants, and run labels</summary>
    
```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from io import StringIO

DAY = 86400.0
TARGET_XC = 0.5

runs = {
    "step ov": Path("LOGS_step_ov"),
    "exp ov":  Path("LOGS_exp_ov"),
    "PC":      Path("LOGS_PC"),
}

styles = {
    "step ov": dict(color="C0", ls="-"),
    "exp ov":  dict(color="C1", ls="--"),
    "PC":      dict(color="C2", ls=":"),
}
```

</details>

### Code Block 2: Read a MESA History or Profile File

MESA history and profile files have a header section, then a blank line, then the main data table. This helper function reads the main table into a pandas DataFrame.

<details>
    <summary>Python: helper function for MESA tables</summary>
    
```python
def read_mesa_table(path):

    path = Path(path)
    lines = path.read_text(errors="replace").splitlines()

    blank = next(i for i, line in enumerate(lines) if not line.strip())

    cols = lines[blank + 2].split()
    data = "\n".join(lines[blank + 3:])

    return pd.read_csv(StringIO(data), sep=r"\s+", names=cols)
```

</details>

### Code Block 3: Read profiles.index

The history file tells us which model_number is closest to `Xc(H)=0.5`. The `profiles.index` file tells us which profile file corresponds to that model.

<details>
    <summary>Python: read profiles.index</summary>
    
```python
def read_profiles_index(logdir):
    """
    Read profiles.index or profile.index in a LOGS directory.

    The columns are:
        model_number, priority, profile_number
    """
    for name in ["profiles.index", "profile.index"]:
        path = logdir / name
        if path.exists():
            return pd.read_csv(
                path,
                sep=r"\s+",
                skiprows=1,
                names=["model_number", "priority", "profile_number"],
                comment="#",
            )

    raise FileNotFoundError(
        f"Cannot find profiles.index or profile.index in {logdir}"
    )
```

</details>

### Code Block 4: Find the History and Profile File

<details>
    <summary>Python: locate history and profile files</summary>
    
```python
def find_history_file(logdir):

    candidates = sorted(logdir.glob("*MS.history"))

    if not candidates:
        candidates = sorted(logdir.glob("*.history"))

    if not candidates:
        raise FileNotFoundError(f"No history file found in {logdir}")

    return candidates[0]

def find_profile_file(logdir, profile_number):

    profile_number = int(profile_number)

    candidates = [
        logdir / f"profile{profile_number}.data",
        logdir / f"profile{profile_number:05d}.data",
    ]

    for path in candidates:
        if path.exists():
            return path

    matches = sorted(logdir.glob(f"profile*{profile_number}*.data"))

    if matches:
        return matches[0]

    raise FileNotFoundError(
        f"Cannot find profile file for profile_number={profile_number} in {logdir}"
    )
```

</details>

### Code Block 5: Find the Profile Closest to Xc(H) = 0.5

<details>
    <summary>Python: select the profile closest to Xc(H)=0.5</summary>
    
```python
def find_profile_at_xc(logdir, target=0.5):

    hist_file = find_history_file(logdir)
    hist = read_mesa_table(hist_file)
    pidx = read_profiles_index(logdir)

    if "center_h1" not in hist.columns:
        raise KeyError(f"{hist_file} does not contain center_h1")

    if "model_number" not in hist.columns:
        raise KeyError(f"{hist_file} does not contain model_number")

    idx = np.argmin(np.abs(hist["center_h1"].to_numpy() - target))
    row_h = hist.iloc[idx]
    target_model_number = int(row_h["model_number"])

    matched = pidx[pidx["model_number"] == target_model_number]

    if len(matched) == 0:
        j = np.argmin(
            np.abs(pidx["model_number"].to_numpy() - target_model_number)
        )
        matched = pidx.iloc[[j]]

    row_p = matched.iloc[0]
    saved_model_number = int(row_p["model_number"])
    profile_number = int(row_p["profile_number"])

    prof_file = find_profile_file(logdir, profile_number)
    prof = read_mesa_table(prof_file)

    print(
        f"{logdir.name:12s} | "
        f"history = {hist_file.name:20s} | "
        f"target model = {target_model_number:6d} | "
        f"saved model = {saved_model_number:6d} | "
        f"profile = {prof_file.name:16s} | "
        f"center_h1 = {row_h['center_h1']:.6f}"
    )

    return prof
```

</details>

### Code Block 6: Make the Comparison Plot

This final block loops over the three mixing prescriptions, reads the selected profile, and plots:

`h1` and `he4` in the abundance diagram;
`N/2pi` and `S_l/2pi` in the propagation diagram.

The unit conversion for the Brunt-Vaisala frequency is:

```text
brunt_N2:        rad^2/s^2
sqrt(brunt_N2):  rad/s
N/(2pi):         cycles/s
N/(2pi)*86400:   cycles/day
```

For `lamb_Sl1`, MESA stores the Lamb frequency in microHz, so:

```text
lamb_Sl1 * 1e-6 * 86400 = cycles/day
```

<details>
    <summary>Python: make the abundance and propagation plots</summary>
    
```python
fig, (ax_abun, ax_prop) = plt.subplots(
    2,
    1,
    figsize=(7.2, 7.6),
    sharex=True,
    gridspec_kw={"height_ratios": [1.0, 1.15], "hspace": 0.05},
)

for label, logdir in runs.items():
    prof = find_profile_at_xc(logdir, TARGET_XC)

    mfrac = prof["mass"].to_numpy() / prof["mass"].max()

    # Abundance diagram
    ax_abun.plot(
        mfrac,
        prof["h1"],
        lw=1.0,
        alpha=0.6,
        label=rf"{label}: $X_\mathrm{{H}}$",
        **styles[label],
    )

    ax_abun.plot(
        mfrac,
        prof["he4"],
        lw=1.0,
        alpha=0.6,
        label=rf"{label}: $Y_\mathrm{{He}}$",
        **styles[label],
    )

    # Propagation diagram
    small = 1e-30

    N = np.sqrt(np.maximum(prof["brunt_N2"].to_numpy(), small))
    N_cpd = N / (2.0 * np.pi) * DAY

    S1_cpd = np.maximum(prof["lamb_Sl1"].to_numpy(), small) * 1e-6 * DAY

    ax_prop.plot(
        mfrac,
        np.log10(N_cpd),
        lw=1.0,
        alpha=0.6,
        label=rf"{label}: $N/2\pi$",
        **styles[label],
    )

    ax_prop.plot(
        mfrac,
        np.log10(S1_cpd),
        lw=1.0,
        alpha=0.6,
        label=rf"{label}: $S_{{\ell=1}}/2\pi$",
        **styles[label],
    )

ax_abun.set_ylabel("Mass fraction")
ax_abun.set_ylim(-0.03, 1.03)
ax_abun.legend(frameon=False, fontsize=10, ncol=2)

ax_prop.set_xlabel(r"$m/M_\star$")
ax_prop.set_ylabel(r"$\log_{10}(\mathrm{frequency}/\mathrm{day}^{-1})$")
ax_prop.set_ylim(-0.5, 2.2)
ax_prop.legend(frameon=False, fontsize=10, ncol=2)
fig.savefig("compare_XcH050_structure.png", dpi=300, bbox_inches="tight")
```

</details>

### Example Output

The figure below shows an example comparison at approximately `Xc(H) = 0.5`. The upper panel compares the hydrogen and helium abundance profiles, while the lower panel shows the corresponding propagation diagram for the three mixing prescriptions. The complete solution script for Task 8 is provided as `diff_mixing_profiles_in_mass.py`.

<img src="https://github.com/astroscien/2026MESA-school-day2-lab3/blob/main/compare_XcH050_structure_mass_fraction_f0p86.png?raw=true" width="750">

---

## References

- [Johnston et al. (2024)](https://ui.adsabs.harvard.edu/abs/2024ApJ...964..170J/abstract), *Modelling Time-dependent Convective Penetration in 1D Stellar Evolution*, ApJ, 964, 170.
