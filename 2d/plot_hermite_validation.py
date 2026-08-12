from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"


def read_tecplot(path):
    """Read the block-format Tecplot files written by output2d.f90."""
    lines = Path(path).read_text().splitlines()
    zone_index = next(i for i, line in enumerate(lines) if line.lstrip().startswith("ZONE"))
    match = re.search(r"I\s*=\s*(\d+).*J\s*=\s*(\d+)", lines[zone_index])
    ni, nj = map(int, match.groups())
    values = np.fromstring(" ".join(lines[zone_index + 1 :]), sep=" ")

    node_count = ni * nj
    nx, ny = ni - 1, nj - 1
    cell_count = nx * ny
    expected = 2 * node_count + 11 * cell_count
    if values.size != expected:
        raise ValueError(f"{path}: expected {expected} values, found {values.size}")

    xnode = values[:node_count].reshape((ni, nj), order="F")
    ynode = values[node_count : 2 * node_count].reshape((ni, nj), order="F")
    raw = values[2 * node_count :].reshape((11, cell_count))
    fields = [raw[k].reshape((nx, ny), order="F") for k in range(11)]
    x = 0.25 * (xnode[:-1, :-1] + xnode[1:, :-1] + xnode[:-1, 1:] + xnode[1:, 1:])
    y = 0.25 * (ynode[:-1, :-1] + ynode[1:, :-1] + ynode[:-1, 1:] + ynode[1:, 1:])
    names = ("rho", "u", "v", "T", "p", "qx", "qy", "radius", "radius0", "rho0", "error")
    return x, y, dict(zip(names, fields))


def numbered_files(directory):
    return sorted(directory.glob("macro*.dat"), key=lambda p: int(re.search(r"\d+", p.stem).group()))


def save_all(fig, stem):
    FIGURES.mkdir(exist_ok=True)
    for suffix in ("png", "pdf", "eps"):
        fig.savefig(FIGURES / f"{stem}.{suffix}", dpi=600, bbox_inches="tight")


def plot_poiseuille():
    steady = RESULTS / "poiseuille_steady"
    directory = steady if numbered_files(steady) else RESULTS / "poiseuille_fast"
    path = numbered_files(directory)[-1]
    _, y, field = read_tecplot(path)
    yc = y.mean(axis=0)

    fig, axes = plt.subplots(1, 3, figsize=(8.2, 3.1))
    for ax, key, label in zip(axes, ("rho", "u", "T"), (r"$\rho$", r"$U$", r"$T$")):
        value = field[key].mean(axis=0)
        ax.plot(yc, value, "o--", ms=2.8, mfc="white", mec="black", color="black", lw=0.8)
        ax.set_xlabel(r"$Y$")
        ax.set_ylabel(label)
        ax.grid(True, color="0.8", ls="--", lw=0.5)
        ax.set_xlim(0.0, 1.0)
    fig.tight_layout()
    save_all(fig, "poiseuille_hermite_kn002")
    plt.close(fig)

    u = field["u"].mean(axis=0)
    print(f"Poiseuille file: {path}")
    for key in ("rho", "u", "T"):
        print(f"  {key} range = [{np.min(field[key]):.6e}, {np.max(field[key]):.6e}]")
    print(f"  max |V|                = {np.max(np.abs(field['v'])):.6e}")
    print(f"  max x variation of U   = {np.max(np.ptp(field['u'], axis=0)):.6e}")
    print(f"  midplane symmetry of U = {np.max(np.abs(u-u[::-1])):.6e}")


def plot_rt():
    production = RESULTS / "rt_120"
    directory = production if len(numbered_files(production)) >= 4 else RESULTS / "rt_60_fast"
    files = numbered_files(directory)
    if len(files) >= 5:
        files = [files[0], files[1], files[2], files[-1]]
    times = (0.0, 0.04, 0.08, 0.14)

    records = [read_tecplot(path) for path in files]
    vmin, vmax = 0.0, 0.45
    fig, axes = plt.subplots(1, len(records), figsize=(3.0 * len(records), 2.8), constrained_layout=True)
    axes = np.atleast_1d(axes)
    contour = None
    for ax, (x, y, field), time in zip(axes, records, times):
        contour = ax.pcolormesh(x, y, field["rho"], shading="nearest", cmap="jet", vmin=vmin, vmax=vmax)
        ax.set_aspect("equal")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax.text(0.04, 0.94, rf"$t={time:g}$", transform=ax.transAxes, color="white", va="top")
        ax.set_xlabel(r"$x$")
    axes[0].set_ylabel(r"$y$")
    fig.colorbar(contour, ax=axes, shrink=0.82)
    save_all(fig, f"{directory.name}_density_hermite")
    plt.close(fig)

    final = records[-1][2]
    print(f"RT directory: {directory}")
    print(f"  density range = [{np.min(final['rho']):.6e}, {np.max(final['rho']):.6e}]")
    print(f"  temperature range = [{np.min(final['T']):.6e}, {np.max(final['T']):.6e}]")
    print(f"  all fields finite = {all(np.all(np.isfinite(v)) for v in final.values())}")


if __name__ == "__main__":
    plot_poiseuille()
    plot_rt()
    print(f"Figures written to {FIGURES}")
