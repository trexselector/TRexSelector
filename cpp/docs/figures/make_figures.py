#!/usr/bin/env python3
"""Draw the figures of the C++ user guide from the data the example programs write into
examples/output, with the matplotlib style of the Python guide (its guide/matplotlibrc). The PNG files
go to pages/figures and are committed, so the pages build without matplotlib as well.

    python3 cpp/docs/figures/make_figures.py          # all figures
"""
import pathlib
import sys

import matplotlib
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
DOCS = HERE.parent
OUT = DOCS / "examples" / "output"
FIG = DOCS / "pages" / "figures"
matplotlib.use("Agg")
matplotlib.rc_file(DOCS.parents[1] / "Python" / "trexselector" / "docs" / "guide" / "matplotlibrc")
import matplotlib.pyplot as plt  # noqa: E402  (after the style is loaded)


def read_meta(path):
    """Key-value lines, the value kept as a string."""
    out = {}
    for line in path.read_text().splitlines():
        if line.strip():
            key, _, value = line.partition(" ")
            out[key] = value.strip()
    return out


def tsolvers_path():
    """Page 9: two snapshots of one warm-started T-LARS path."""
    meta = read_meta(OUT / "tsolvers.figure.meta.txt")
    p = int(meta["p"])
    true_support = [int(v) for v in meta["true_support"].split()]

    def load(tag):
        # the file holds one row per step, the all-zero start first; the figure wants variables in rows
        path = np.loadtxt(OUT / f"tsolvers.figure_{tag}.csv", delimiter=",", ndmin=2).T
        codes = (OUT / f"tsolvers.figure_{tag}.actions.txt").read_text().split()
        entering = [abs(int(c)) - 1 for c in codes]     # lars codes +-(j + 1) to the 0-based j
        return path, entering

    def draw(ax, path, entering, title):
        steps = np.arange(path.shape[1])
        for j in np.flatnonzero(np.any(path != 0, axis=1)):
            if j >= p:
                ax.plot(steps, path[j], color="C1", linestyle="--", linewidth=1.2)
            elif j in true_support:
                ax.plot(steps, path[j], color="C0", linewidth=1.8)
            else:
                ax.plot(steps, path[j], color="0.45", linewidth=1.0)
        for step, j in enumerate(entering, start=1):
            if j >= p:
                ax.axvline(step, color="C1", linestyle=":", linewidth=0.9)
        top = ax.secondary_xaxis("top")
        top.set_xticks(range(1, len(entering) + 1))
        top.set_xticklabels([f"D{j - p}" if j >= p else str(j) for j in entering], fontsize=7)
        top.tick_params(length=0)
        top.spines["top"].set_visible(False)
        ax.set_title(title)
        ax.set_xlabel("Step")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    draw(ax1, *load("T1"), "T_stop = 1")
    draw(ax2, *load("T3"), "T_stop = 3, continued")
    ax1.set_ylabel("Coefficient")
    handles = [plt.Line2D([], [], color="C0", linewidth=1.8, label="true predictor"),
               plt.Line2D([], [], color="0.45", linewidth=1.0, label="other predictor"),
               plt.Line2D([], [], color="C1", linestyle="--", label="dummy"),
               plt.Line2D([], [], color="C1", linestyle=":", linewidth=0.9, label="dummy enters")]
    fig.legend(handles=handles, loc="outside center right")
    fig.savefig(FIG / "tsolvers_path.png")
    plt.close(fig)
    return "tsolvers_path.png"


def trex_selector_phi():
    """Page 3: the relative occurrences after the last stopping time with the voting threshold."""
    rows = np.loadtxt(OUT / "trex_selector.figure_phi.csv", delimiter=",", ndmin=2)
    meta = read_meta(OUT / "trex_selector.figure_phi.meta.txt")
    variable, phi, is_true = rows[:, 0], rows[:, 1], rows[:, 2] > 0
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.vlines(variable[is_true], 0, phi[is_true], colors="C0", linewidth=2, label="true support")
    ax.vlines(variable[~is_true], 0, phi[~is_true], colors="0.6", linewidth=2, label="other")
    ax.axhline(float(meta["v_thresh"]), linestyle="--", color="C1", label="voting threshold")
    ax.set_xlabel("Variable")
    ax.set_ylabel("Relative occurrence")
    ax.legend()
    fig.savefig(FIG / "trex_selector_phi.png")
    plt.close(fig)
    return "trex_selector_phi.png"


def trex_selector_mc():
    """Page 3: the estimated FDR against the target and the TPR of the small simulation."""
    rows = np.loadtxt(OUT / "trex_selector.figure_mc.csv", delimiter=",", ndmin=2)
    target, fdr, tpr = rows[:, 0], rows[:, 1], rows[:, 2]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.5))
    ax1.plot(target, fdr, "o-", label="mean over 20 trials")
    ax1.plot([0, 0.35], [0, 0.35], "--", color="C1", label="FDR equal to the target")
    ax1.set_xlim(0.05, 0.35)
    ax1.set_ylim(0, 0.5)
    ax1.set_xlabel("Target FDR")
    ax1.set_ylabel("FDR")
    ax2.plot(target, tpr, "o-")
    ax2.set_xlim(0.05, 0.35)
    ax2.set_ylim(0, 1.05)
    ax2.set_xlabel("Target FDR")
    ax2.set_ylabel("TPR")
    fig.legend(*ax1.get_legend_handles_labels(), loc="outside center right")
    fig.savefig(FIG / "trex_selector_mc.png")
    plt.close(fig)
    return "trex_selector_mc.png"


def trex_spca_loadings():
    """Page 7: the two sparse loading vectors."""
    rows = np.loadtxt(OUT / "trex_spca.figure_loadings.csv", delimiter=",", ndmin=2)
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))
    variable = rows[:, 0]
    supports = [range(0, 10), range(30, 40)]
    for m, ax in enumerate(axes):
        on = np.isin(variable, list(supports[m]))
        ax.vlines(variable[on], 0, rows[on, m + 1], color="C0", linewidth=2, label="true support")
        ax.vlines(variable[~on], 0, rows[~on, m + 1], color="C1", linewidth=2, label="other")
        ax.axhline(0, color="0.3", linewidth=0.5)
        ax.set_xlabel("Variable")
        ax.set_ylabel(f"Loading {m + 1}")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside center right")
    fig.savefig(FIG / "trex_spca_loadings.png")
    plt.close(fig)
    return "trex_spca_loadings.png"


def ml_methods_ridge():
    """Page 10: the cross-validation error of the ridge penalty with the two choices."""
    rows = np.loadtxt(OUT / "ml_methods.figure_ridge.csv", delimiter=",", ndmin=2)
    meta = read_meta(OUT / "ml_methods.figure_ridge.meta.txt")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(rows[:, 0], rows[:, 1], label="CV error, 5 folds")
    ax.axvline(float(meta["cv_min"]), color="C1", linestyle="-", label="minimum")
    ax.axvline(float(meta["cv_1se"]), color="C1", linestyle="--", label="one standard error")
    ax.set_xscale("log")
    ax.set_xlabel("Lambda")
    ax.set_ylabel("CV error")
    ax.legend()
    fig.savefig(FIG / "ml_methods_ridge.png")
    plt.close(fig)
    return "ml_methods_ridge.png"


FIGURES = [trex_selector_phi, trex_selector_mc, trex_spca_loadings, tsolvers_path, ml_methods_ridge]

if __name__ == "__main__":
    FIG.mkdir(parents=True, exist_ok=True)
    for make in FIGURES:
        print("wrote", make())
    sys.exit(0)
