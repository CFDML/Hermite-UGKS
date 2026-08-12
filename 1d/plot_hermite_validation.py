"""Plot one Hermite/Gauss-Hermite result for each 1-D paper case."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent

# Final output from the three paper-matched runs.  Change only these paths when
# a case is rerun with a different final iteration number.
HYDRO_FILE = ROOT / "results/hydrostatic_kn1_gh36_release/macro944.dat"
SOD_FILE = ROOT / "results/sod_kn001_gh34_release/macro2309.dat"
FOURIER_FILE = ROOT / "results/fourier_force2_kn01_gh22_release2/macro213439.dat"


def read_tecplot_block(path):
    """Read the simple BLOCK-format files written by src/output.f90."""
    with path.open(encoding="utf-8") as stream:
        names = [
            name.strip()
            for name in stream.readline().split("=", 1)[1].split(",")
        ]
        nx = int(stream.readline().split("I =", 1)[1].split(",", 1)[0])
        values = np.fromstring(stream.read(), sep=" ")

    expected = len(names) * nx
    if values.size != expected:
        raise ValueError(f"{path}: found {values.size} values, expected {expected}")
    blocks = values.reshape(len(names), nx)
    return {name.upper(): blocks[i] for i, name in enumerate(names)}


hydro = read_tecplot_block(HYDRO_FILE)
sod = read_tecplot_block(SOD_FILE)
fourier = read_tecplot_block(FOURIER_FILE)

fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.35))

axes[0].plot(hydro["X"], hydro["PERTURBATION"], color="black", lw=1.5)
axes[0].set(xlabel=r"$x$", ylabel=r"$p-p_{\rm ref}$")
axes[0].text(0.05, 0.92, r"$N_m=6$, GH36", transform=axes[0].transAxes)

axes[1].plot(sod["X"], sod["RHO"], color="black", lw=1.5)
axes[1].set(xlabel=r"$x$", ylabel=r"$\rho$")
axes[1].text(0.05, 0.92, r"$N_m=6$, GH34", transform=axes[1].transAxes)

axes[2].plot(fourier["X"], fourier["RHO"], color="black", lw=1.5)
axes[2].set(xlabel=r"$x$", ylabel=r"$\rho$")
axes[2].text(
    0.05, 0.92, r"$F=2$, $N_m=10$, GH22", transform=axes[2].transAxes
)

for label, ax in zip(("(a)", "(b)", "(c)"), axes):
    ax.set_xlim(0.0, 1.0)
    ax.tick_params(direction="in", top=True, right=True)
    ax.text(0.5, -0.25, label, ha="center", transform=ax.transAxes)

fig.tight_layout()
out_dir = ROOT / "figures"
out_dir.mkdir(exist_ok=True)
fig.savefig(out_dir / "hermite_1d_validation.png", dpi=600, bbox_inches="tight")
fig.savefig(out_dir / "hermite_1d_validation.pdf", bbox_inches="tight")
plt.close(fig)

print(out_dir / "hermite_1d_validation.png")
