"""Shared figure style for the supplementary figure scripts.

Keeps the supplementary panels visually consistent with the main figures:
same family, same sizes, same frame treatment, same neutral grey for
annotation. Self-contained, so the archive does not depend on anything outside
it.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

META_GREY = "#666666"
OKABE_ITO = ["#009E73", "#0072B2", "#D55E00", "#CC79A7",
             "#E69F00", "#56B4E9", "#F0E442", "#000000"]


def apply_figure_style(sizes=(8, 8.5, 7.5)):
    """sizes = (base font, axis-label font, tick font)."""
    base, label, tick = sizes
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": base,
        "axes.labelsize": label,
        "axes.titlesize": label + 0.5,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "legend.fontsize": tick,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "savefig.dpi": 300,
        "figure.dpi": 110,
        "axes.prop_cycle": mpl.cycler(color=OKABE_ITO),
    })


def set_frame(ax, left=True, bottom=True):
    """Keep only the axes a reader needs to decode the plot."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    ax.tick_params(length=3, pad=2)


def panel_letter(ax, letter, x=-0.12, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="left")


def panel_crops(*args, **kwargs):
    """No-op kept for call compatibility with the figure scripts."""
    return None
