"""Draw the diagram of the memory-mapping page: a design in memory and a design in a file reach the
selectors, the T-solvers and the ML methods through the same mapped view, so one implementation
serves both and the results agree.

usage: python memory_agnostic.py <output folder>
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

TEAL_DARK, TEAL, TEAL_LIGHT, TINT = "#336C68", "#6BB2AD", "#BFDDDA", "#E6F1F0"
INK, MUTED = "#24292f", "#57606a"

out = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
out.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 10.5, "font.family": "DejaVu Sans"})
fig, ax = plt.subplots(figsize=(9.0, 3.6))
ax.set_xlim(0, 100)
ax.set_ylim(0, 40)
ax.axis("off")


def box(x, y, w, h, title, sub, fill, edge=TEAL_DARK, title_color=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.6",
                                fc=fill, ec=edge, lw=1.4, zorder=2))
    ax.text(x + w / 2, y + h * 0.66, title, ha="center", va="center", fontsize=11,
            fontweight="bold", color=title_color, zorder=3)
    ax.text(x + w / 2, y + h * 0.30, sub, ha="center", va="center", fontsize=9, color=MUTED,
            zorder=3, linespacing=1.35)


def arrow(x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=1,
                arrowprops=dict(arrowstyle="-|>", color=TEAL_DARK, lw=1.6, shrinkA=0, shrinkB=0,
                                mutation_scale=14))


# The two sources on the left.
box(2, 23, 30, 13, "Design in memory", "an R matrix, a NumPy array,\nan Eigen matrix", TINT)
box(2, 4, 30, 13, "Design in a file", "a memory-mapped file\n(Boost.Interprocess),\npages read on demand",
    TINT)

# The one view and the one implementation.
box(42, 13.5, 24, 13, "One mapped view", "Eigen::Map of the columns,\nno copy of the design", TEAL_LIGHT)
box(74, 13.5, 24, 13, "One implementation", "selectors, T-solvers,\nML methods", TEAL,
    title_color="white")

arrow(32, 29.5, 42, 22.5)
arrow(32, 10.5, 42, 17.5)
arrow(66, 20, 74, 20)

ax.text(54, 36.5, "the same code path, the same results", ha="center", va="center", fontsize=10,
        color=TEAL_DARK, fontweight="bold")
ax.text(86, 8.5, "where the matrix lives\nmakes no difference", ha="center", va="center",
        fontsize=9, color=MUTED, linespacing=1.35)

fig.tight_layout(pad=0.4)
fig.savefig(out / "memory_agnostic.svg")
fig.savefig(out / "memory_agnostic.png", dpi=200)
print("written", out / "memory_agnostic.svg", out / "memory_agnostic.png")
