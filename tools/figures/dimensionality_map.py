"""Draw the regime map of the dimensionality page: the (p, n) plane up to p = 1e8 with the boundaries
p = n and p = 1000 n, the three regimes shaded in one hue from light to dark, and typical settings as
markers.

usage: python dimensionality_map.py <output folder>
"""
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TEAL_DARK, TEAL, TEAL_LIGHT, TINT = "#336C68", "#6BB2AD", "#BFDDDA", "#E6F1F0"
INK, MUTED = "#24292f", "#57606a"

out = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
out.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})
fig, ax = plt.subplots(figsize=(7.4, 6.0))
LO, P_HI, N_HI = 1e1, 1e8, 1e6
RULE = 1e3  # the dashed boundary p = RULE * n
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(LO, P_HI)
ax.set_ylim(LO, N_HI)

# The three regimes as polygons in data coordinates (the axes clip them).
ax.fill([LO, N_HI, LO], [LO, N_HI, N_HI], color=TINT, lw=0, zorder=0)                         # n >= p
ax.fill([LO, P_HI, P_HI, LO * RULE], [LO, P_HI, P_HI / RULE, LO], color=TEAL_LIGHT, lw=0, zorder=0)  # n < p <= 1000 n
ax.fill([LO * RULE, P_HI, P_HI], [LO, P_HI / RULE, LO], color=TEAL, lw=0, zorder=0)           # p > 1000 n

# The boundaries.
ax.plot([LO, N_HI], [LO, N_HI], color=TEAL_DARK, lw=1.8, zorder=2)
ax.plot([LO * RULE, P_HI], [LO, P_HI / RULE], color=TEAL_DARK, lw=1.4, ls=(0, (5, 4)), zorder=2)


def boundary_label(x, y, text, va):
    """A label along a boundary of slope one in the log-log plane, rotated by the on-screen angle of
    that slope, which depends on the axes' aspect (set after the layout, see below)."""
    return ax.text(x, y, text, color=TEAL_DARK, rotation_mode="anchor", ha="left", va=va, fontsize=10.5)


boundary_labels = [boundary_label(1.5e5, 1.5e5 * 1.25, "p = n", "bottom"),
                   boundary_label(1.5e5, 1.5e2 / 1.3, "p = 1000 n", "top")]  # below the dashed line

# Regime labels: the title, and the keywords a factor `gap` below it (log axis).
def regime(x, y, title, sub, sub_color=MUTED, gap=3.0):
    ax.text(x, y, title, color=INK, ha="center", va="center", fontsize=12.5, fontweight="bold", zorder=3)
    ax.text(x, y / gap, sub, color=sub_color, ha="center", va="top", fontsize=9.5, zorder=3, linespacing=1.4)

def shift(percent):
    """A horizontal shift by a share of the plot width, seven decades wide, as a factor in p."""
    return 10 ** (percent / 100 * 7)


def lift(percent):
    """A vertical shift by a share of the plot height, five decades tall, as a factor in n."""
    return 10 ** (percent / 100 * 5)


regime(1.2e2 * shift(15), 1.5e5, "Low dimensionality\n$n \\geq p$",
       "least squares, p-values,\nclassical multiple testing")
regime(4.0e4 * shift(10), 1.0e4, "High dimensionality\n$p > n$",
       "sparsity, path solvers,\ndummies for FDR control", gap=2.5)
regime(5.0e6, 1.6e2, "Highest\ndimensionality\n$p \\gg n$",
       "screening, memory mapping,\nthread budgets", sub_color=INK, gap=5.5)

# Typical settings: the label beside the marker, or above it where a boundary runs close by.
def point(x, y, label, dx=1.5, ha="left", above=False):
    ax.plot([x], [y], marker="o", ms=7, mfc="white", mec=INK, mew=1.4, zorder=4)
    if above:
        ax.text(x, y * 1.7, label, color=INK, ha="center", va="bottom", fontsize=9.5, zorder=4)
    else:
        ax.text(x * dx, y, label, color=INK, ha=ha, va="center", fontsize=9.5, zorder=4)

point(4e1, 1e4, "classical regression")
point(3e2, 2e3, "financial factors", above=True)
point(8e7, 8e4, "genomics biobank", dx=0.65, ha="right")  # millions of variants, about 1e5 participants

ax.set_xlabel("Number of variables $p$", color=INK, fontsize=11.5)
ax.set_ylabel("Number of observations $n$", color=INK, fontsize=11.5)
ax.tick_params(colors=MUTED, labelsize=9.5)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
for side in ("left", "bottom"):
    ax.spines[side].set_color(MUTED)
ax.grid(False)
fig.tight_layout()
fig.canvas.draw()  # fixes the axes position, so the on-screen slope of the boundaries is known
(x0, y0), (x1, y1) = ax.transData.transform([[1e2, 1e2], [1e5, 1e5]])
slope = math.degrees(math.atan2(y1 - y0, x1 - x0))
for label in boundary_labels:
    label.set_rotation(slope)
fig.savefig(out / "dimensionality_map.svg")
fig.savefig(out / "dimensionality_map.png", dpi=200)
print("written", out / "dimensionality_map.svg", out / "dimensionality_map.png")
