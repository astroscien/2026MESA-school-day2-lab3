import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from io import StringIO

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 16,
    "font.size": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
})

DAY = 86400.0 # day to sec
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


def read_mesa_table(path):
    """
    Read a MESA history/profile file with standard MESA header format.
    """
    path = Path(path)
    lines = path.read_text(errors="replace").splitlines()

    blank = next(i for i, line in enumerate(lines) if not line.strip())

    cols = lines[blank + 2].split()
    data = "\n".join(lines[blank + 3:])

    return pd.read_csv(StringIO(data), sep=r"\s+", names=cols)


def read_profiles_index(logdir):
    """
    Read profiles.index or profile.index in a LOGS directory.
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


def find_history_file(logdir):
    """
    Find the main-sequence history file in a LOGS directory.
    Prefer *MS.history if available.
    """
    candidates = sorted(logdir.glob("*MS.history"))

    if not candidates:
        candidates = sorted(logdir.glob("*.history"))

    if not candidates:
        raise FileNotFoundError(f"No history file found in {logdir}")

    return candidates[0]


def find_profile_file(logdir, profile_number):
    """
    Find the profile file corresponding to a profile number.
    """
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


def find_profile_at_xc(logdir, target=0.5):
    """
    Find the profile closest to target center_h1.

    Logic:
    1. Read the MS history file.
    2. Find model_number where center_h1 is closest to target.
    3. Use profiles.index to map model_number to profile_number.
    4. Read the corresponding profile file.
    """
    hist_file = find_history_file(logdir)
    hist = read_mesa_table(hist_file)
    pidx = read_profiles_index(logdir)

    if "center_h1" not in hist.columns:
        raise KeyError(f"{hist_file} does not contain center_h1")

    if "model_number" not in hist.columns:
        raise KeyError(f"{hist_file} does not contain model_number")

    # Find the target model in history.
    idx = np.argmin(np.abs(hist["center_h1"].to_numpy() - target))
    row_h = hist.iloc[idx]
    target_model_number = int(row_h["model_number"])

    # Find the corresponding saved profile.
    matched = pidx[pidx["model_number"] == target_model_number]

    if len(matched) == 0:
        # If the exact model was not saved as a profile, use the nearest saved model.
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


def main():
    fig, (ax_abun, ax_prop) = plt.subplots(
        2,
        1,
        figsize=(7.2, 7.6),
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.15], "hspace": 0.05},
    )

    for label, logdir in runs.items():
        prof = find_profile_at_xc(logdir, TARGET_XC)

        r = prof["radius"].to_numpy()

        # Abundance diagram
        ax_abun.plot(
            r,
            prof["h1"],
            lw=2.0,
            label=rf"{label}: $X_\mathrm{{H}}$",
            **styles[label],
        )

        ax_abun.plot(
            r,
            prof["he4"],
            lw=1.5,
            alpha=0.7,
            label=rf"{label}: $Y_\mathrm{{He}}$",
            **styles[label],
        )

        # Propagation diagram
        small = 1e-30

        # brunt_N2 is in rad^2/s^2, so sqrt gives rad/s.
        N = np.sqrt(np.maximum(prof["brunt_N2"].to_numpy(), small))
        N_cpd = N / (2.0 * np.pi) * DAY

        # lamb_Sl1 is usually in microHz, i.e. 1e-6 cycles/s.
        S1_cpd = np.maximum(prof["lamb_Sl1"].to_numpy(), small) * 1e-6 * DAY

        ax_prop.plot(
            r,
            np.log10(N_cpd),
            lw=2.0,
            label=rf"{label}: $N/2\pi$",
            **styles[label],
        )

        ax_prop.plot(
            r,
            np.log10(S1_cpd),
            lw=1.5,
            alpha=0.7,
            label=rf"{label}: $S_{{\ell=1}}$",
            **styles[label],
        )

    ax_abun.set_ylabel("Mass fraction")
    ax_abun.set_ylim(-0.03, 1.03)
    ax_abun.legend(frameon=False, fontsize=10, ncol=2)

    ax_prop.set_xlabel(r"$r/R_\odot$")
    ax_prop.set_ylabel(r"$\log_{10}(\mathrm{frequency}/\mathrm{day}^{-1})$")
    ax_prop.set_ylim(-0.5, 2.2)
    ax_prop.legend(frameon=False, fontsize=10, ncol=2)

    fig.savefig("compare_XcH050_structure.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
