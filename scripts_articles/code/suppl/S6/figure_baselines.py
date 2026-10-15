"""x60_figure.py -- two-panel figure: Pi vs reference rankings.

Panel a: simulated benchmark, top-1 rate pooled over the 10 positive scenarios.
Panel b: rank assigned to the true reservoir in each of the three real outbreaks.
"""
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from figstyle import apply_figure_style, panel_letter, set_frame, META_GREY
apply_figure_style(sizes=(8, 7, 6))

tab = pd.read_csv("../../../data/suppl_S6_baselines/pi_vs_baselines_table.csv")
per = pd.read_csv("../../../data/suppl_S6_baselines/bench_per_scenario.csv")
N_REP = int(tab["n_rep_per_scenario"].iloc[0])

SHORT = {
    "pi": "Composite $\\Pi$",
    "pcr": "Animal PCR prevalence",
    "sero": "Animal seroprevalence",
    "C7": "C7 catalytic FOI",
    "C14": "C14 spatial independence",
    "C5": "C5 exposure fraction",
    "or": "Exposure odds ratio",
    "C12": "C12 inter-epidemic",
    "taxo": "Taxonomic order (proxy)",
}
FOCAL = "#1a4f8b"     # composite Pi
OTHER = "#b0b7bf"     # simple reference rankings
PROXY = "#d98d3a"     # the caricature proxy
ALARM = "#c0392b"

def bar_color(b):
    return FOCAL if b == "pi" else (PROXY if b == "taxo" else OTHER)

fig = plt.figure(figsize=(7.2, 4.3))
gs = fig.add_gridspec(1, 2, width_ratios=[1.30, 1.0], wspace=0.55,
                      left=0.215, right=0.965, top=0.795, bottom=0.235)

# ---------------- Panel a: simulated benchmark -------------------------------
axa = fig.add_subplot(gs[0, 0])
pa = tab.sort_values("sim_top1_pct", ascending=True, na_position="first").copy()
ypos = np.arange(len(pa))
vals = pa["sim_top1_pct"].to_numpy(dtype=float)
plotted = np.where(np.isnan(vals), 0.0, vals)

for y, b, v in zip(ypos, pa["baseline"], plotted):
    axa.barh(y, v, height=0.68, color=bar_color(b),
             edgecolor="white", linewidth=0.6, zorder=2)
# zero-valued bar gets a visible stub (figure-style 6.1)
for y, v in zip(ypos, plotted):
    if v <= 0:
        axa.plot([0], [y], marker="|", color="#4a4a4a", ms=7, mew=1.3, zorder=3)

# taxonomic proxy is not evaluable in silico -> hatched ghost + n.d.
na_idx = np.where(pa["sim_top1_pct"].isna().to_numpy())[0]
for y in na_idx:
    axa.barh(y, 8, height=0.68, facecolor="none", edgecolor=PROXY,
             linewidth=0.8, hatch="////", zorder=2, alpha=0.85)
    axa.text(10.5, y + 0.06, "n.d. — simulated hosts carry no taxonomy",
             va="bottom", ha="left", fontsize=6, color=PROXY, zorder=4)

# top-2 as an open marker on the same row
axa.plot(pa["sim_top2_pct"], ypos, ls="none", marker="o", mfc="white",
         mec="#33383d", mew=0.9, ms=4.2, zorder=4)
# direct-label the top-2 glyph on the row where it is most separated from its bar
_i = int(np.nanargmax(pa["sim_top2_pct"].to_numpy(dtype=float) - plotted))
axa.annotate("top-2 rate", xy=(pa["sim_top2_pct"].iloc[_i], ypos[_i]),
             xytext=(pa["sim_top2_pct"].iloc[_i] - 1.0, ypos[_i] + 0.62),
             fontsize=6, color="#33383d", ha="center", va="center",
             arrowprops=dict(arrowstyle="-", lw=0.6, color="#33383d",
                             shrinkA=1, shrinkB=3), zorder=5)

axa.set_yticks(ypos)
axa.set_yticklabels([SHORT[b] for b in pa["baseline"]])
for tick, b in zip(axa.get_yticklabels(), pa["baseline"]):
    if b == "pi":
        tick.set_color(FOCAL); tick.set_fontweight("bold")
    if b == "taxo":
        tick.set_color(PROXY)
axa.set_xlim(0, 108)
axa.set_xticks([0, 25, 50, 75, 100])
axa.set_xlabel("Reps ranking the true reservoir first (%)")
_n_beat = int((pa["sim_top1_pct"] > float(tab.loc[tab.baseline == "pi",
                                                  "sim_top1_pct"].iloc[0])).sum())
axa.set_title(f"a  In silico, {_n_beat} of 8 reference rankings beat $\\Pi$",
              loc="left", fontweight="normal", pad=15)
axa.margins(y=0.045)
axa.annotate("higher = better", xy=(1.0, -0.185), xycoords="axes fraction",
             fontsize=6, color="#6b7280", ha="right", va="top")

# ---------------- Panel b: three real outbreaks ------------------------------
axb = fig.add_subplot(gs[0, 1])
CASES = [("rank_sn", "Sin Nombre ($P.\\ maniculatus$)"),
         ("rank_nv", "Nipah ($Pteropus$)"),
         ("rank_rv", "Rift Valley fever (wild ungulates)")]
ob = ["pi", "C12", "taxo", "C14", "sero", "pcr", "C7", "C5", "or"]
pb = tab.set_index("baseline").loc[ob]
yb = np.arange(len(pb))[::-1]

axb.axvspan(0.55, 1.45, color="#1a4f8b", alpha=0.055, zorder=0, lw=0)
MARK = ["o", "s", "^"]
for j, (col, _) in enumerate(CASES):
    r = pb[col].to_numpy(dtype=float)
    hit = r == 1
    axb.plot(r[hit], yb[hit], ls="none", marker=MARK[j], ms=5.0,
             mfc=FOCAL, mec=FOCAL, zorder=3)
    axb.plot(r[~hit], yb[~hit], ls="none", marker=MARK[j], ms=5.0,
             mfc="white", mec=ALARM, mew=1.0, zorder=3)

for y, b in zip(yb, pb.index):
    n = int(pb.loc[b, "real_cases_top1_of3"])
    axb.text(3.78, y, f"{n}/3", va="center", ha="right", fontsize=6,
             color=FOCAL if n == 3 else "#4a4a4a",
             fontweight="bold" if n == 3 else "normal")

axb.set_yticks(yb)
axb.set_yticklabels([SHORT[b] for b in pb.index])
for tick, b in zip(axb.get_yticklabels(), pb.index):
    if b == "pi":
        tick.set_color(FOCAL); tick.set_fontweight("bold")
    if b == "taxo":
        tick.set_color(PROXY)
axb.set_xlim(0.4, 3.9)
axb.set_xticks([1, 2, 3])
axb.set_xlabel("Rank given to the true reservoir", labelpad=7)
axb.set_title("b  In vivo, only $\\Pi$ and C12 hold", loc="left",
              fontweight="normal", pad=15)
axb.set_ylim(-0.55, len(pb) - 0.38)
axb.annotate("correct", xy=(1.0, 1.008), xycoords=("data", "axes fraction"),
             fontsize=6, color=FOCAL, ha="center", va="bottom")
axb.annotate("cases solved", xy=(3.78, 1.008),
             xycoords=("data", "axes fraction"),
             fontsize=6, color="#4a4a4a", ha="right", va="bottom")
axb.annotate("lower = better", xy=(1.0, -0.185), xycoords="axes fraction",
             fontsize=6, color="#6b7280", ha="right", va="top")

hnd = [plt.Line2D([], [], ls="none", marker=MARK[j], ms=4.6, mfc=FOCAL,
                  mec=FOCAL, label=lab)
       for j, (_, lab) in enumerate(CASES)]
hnd.append(plt.Line2D([], [], ls="none", marker="o", ms=4.6, mfc="white",
                      mec=ALARM, mew=1.0,
                      label="open marker = reservoir not ranked first"))
fig.legend(handles=hnd, loc="lower left", bbox_to_anchor=(0.215, 0.008),
           frameon=False, ncol=2, handletextpad=0.35, columnspacing=1.6,
           labelspacing=0.30, fontsize=6)

fig.suptitle(
    "Simple reference rankings outperform $\\Pi$ in simulation but fail the "
    "real outbreaks $\\Pi$ was built for",
    x=0.012, y=0.965, ha="left", fontsize=8.6)
fig.text(0.012, 0.895,
         f"Simulated: {N_REP} replicates × 10 positive scenarios "
         f"(n = {N_REP*10} per ranking); filled = true reservoir ranked first.",
         ha="left", fontsize=6.2, color="#5b6167")

fig.savefig("../../../figures/suppl_figure_10_baselines.png", dpi=300, bbox_inches="tight")
fig.savefig("../../../figures/suppl_figure_10_baselines.pdf", bbox_inches="tight")

# ---- §9.1 geometric check ---------------------------------------------------
r = fig.canvas.get_renderer()
texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
         if t.get_text().strip() and t.get_visible()]
spines = [(s, s.get_window_extent(r)) for ax in fig.axes
          for s in ax.spines.values() if s.get_visible()]
tl = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
      for ax in fig.axes}
ov = [(a.get_text(), b.get_text()) for i, (a, ba) in enumerate(texts)
      for b, bb in texts[i+1:] if ba.overlaps(bb)]
ov += [(t.get_text(), "spine") for t, bt in texts for s, bs in spines
       if bt.overlaps(bs) and t not in tl[s.axes]]
oob = [t.get_text() for t, bt in texts if not fig.bbox.containsx(bt.x0)
       or not fig.bbox.containsx(bt.x1) or not fig.bbox.containsy(bt.y0)
       or not fig.bbox.containsy(bt.y1)]
print("OVERLAPS:", ov)
print("OUT_OF_BOUNDS:", oob)
