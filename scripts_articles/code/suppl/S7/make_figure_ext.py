"""
make_figure_ext.py -- Figure for the taxonomic-resolution and
multi-reservoir scenarios.

Four panels:
  a  species-level vs genus-level top-1 accuracy, per scenario
  b  S12: P(exact species correct) as a function of the congeners'
     ecological similarity rho, with and without a resolving assay
  c  S13/S13b: does the score identify the reservoir SET, or arbitrate
     between two co-equal reservoirs?
  d  S14 specificity: distribution of max(Pi) when the true reservoir
     is not sampled, against the null and the single-reservoir arm

Reads ../../../data/suppl_S7_taxonomic/ext_per_replicate.csv and ext_benchmark_results.json.
Writes ../../../figures/suppl_figure_11_taxonomic_resolution.png
"""
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

OUT = "../../../data/suppl_S7_taxonomic"
FIG = "../../../figures"
ev = pd.read_csv(f"{OUT}/ext_per_replicate.csv")
with open(f"{OUT}/ext_benchmark_results.json") as fh:
    res = json.load(fh)

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from figstyle import apply_figure_style, panel_letter, set_frame, META_GREY
apply_figure_style(sizes=(8, 7, 6))

# --- colour threading: one colour per entity, reused in every panel ---
C_SPECIES = "#1f4e79"   # species-level resolution (the strict claim)
C_GENUS   = "#8fb8de"   # genus-level resolution (the weaker claim)
C_XR      = "#c1420f"   # cross-reactive assay  (S12)
C_RES     = "#2a7f62"   # species-resolving assay (S12b)
C_NULL    = "#8a8a8a"   # no-spillover null
C_HIDDEN  = "#c1420f"   # incomplete inventory
C_SINGLE  = "#1f4e79"   # single true reservoir

fig = plt.figure(figsize=(7.2, 7.0))
gs = fig.add_gridspec(2, 2, hspace=0.60, wspace=0.32,
                      left=0.115, right=0.975, top=0.855, bottom=0.075)
axa = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[0, 1])
axc = fig.add_subplot(gs[1, 0])
axd = fig.add_subplot(gs[1, 1])

# ============================ panel a ============================
# Species-exact vs genus-correct top-1, for the arms where a species
# level truth exists.
order_a = ["S1c_ref_single", "S12b_congeners_resolved", "S12_congeners"]
lab_a = ["Single reservoir\n(no congener)",
         "Congeners,\nassay resolves them",
         "Congeners,\ncross-reactive assay"]
sp, gn, nrep_a = [], [], []
for s in order_a:
    d = ev[ev.scenario == s]
    sp.append(d.top1_species.mean() * 100)
    gn.append(d.top1_genus.mean() * 100)
    nrep_a.append(len(d))

y = np.arange(len(order_a))[::-1]
h = 0.34
axa.barh(y + h/2, gn, height=h, color=C_GENUS, edgecolor="none",
         label="correct genus", zorder=2)
axa.barh(y - h/2, sp, height=h, color=C_SPECIES, edgecolor="none",
         label="exact species", zorder=2)
for yi, v in zip(y + h/2, gn):
    axa.text(v + 1.5, yi, f"{v:.0f}", va="center", ha="left", fontsize=6,
             color=C_GENUS)
for yi, v in zip(y - h/2, sp):
    axa.text(v + 1.5, yi, f"{v:.0f}", va="center", ha="left", fontsize=6,
             color=C_SPECIES, fontweight="bold")
axa.set_yticks(y); axa.set_yticklabels(lab_a, fontsize=6)
axa.set_xlim(0, 118); axa.set_xticks([0, 25, 50, 75, 100])
axa.set_xlabel("top-1 accuracy (% of replicates)")
axa.set_title("A congener costs species-level accuracy;\ngenus-level accuracy is unaffected",
              loc="left", fontsize=8)
# direct labels on the topmost pair instead of a legend box
axa.text(gn[0] - 3, y[0] + h/2, "correct genus", va="center", ha="right",
         fontsize=5.8, color="white", fontweight="bold")
axa.text(sp[0] - 3, y[0] - h/2, "exact species", va="center", ha="right",
         fontsize=5.8, color="white", fontweight="bold")
axa.grid(axis="x", lw=0.4, alpha=0.35, zorder=0)
axa.set_axisbelow(True)
panel_letter(axa, "a")

# ============================ panel b ============================
# Which of the 14 criteria actually distinguishes the true congener from
# its sibling? Percentage of replicates in which the criterion scores the
# TRUE species above the congener. 50% = no species-level information;
# below 50% = the criterion systematically prefers the wrong congener.
cs = pd.read_csv(f"{OUT}/ext_criteria_contrast.csv")
cs = cs[cs.scenario == "S12_congeners"].copy()
# effective weights actually applied in these runs (C11 returns NA, its
# weight is redistributed) -- taken from the run metadata, not re-derived
w_eff = res["meta"]["per_scenario"]["S12_congeners"]["weights_used"]
w_map = {f"C{i+1}": w_eff[i] for i in range(14)}
cs["w"] = cs.criterion.map(w_map)
cs = cs.sort_values("pct_favours_true", ascending=True)

yb = np.arange(len(cs))
cols_b = [C_SPECIES if v >= 50 else C_XR for v in cs.pct_favours_true]
axb.hlines(yb, 50, cs.pct_favours_true, color=cols_b, lw=1.1, zorder=2)
axb.scatter(cs.pct_favours_true, yb, s=8 + 240 * cs.w, c=cols_b, zorder=3,
            edgecolors="white", linewidths=0.4)
axb.axvline(50, color="#444444", lw=0.8, ls=(0, (4, 2)), zorder=1)
axb.set_yticks(yb)
axb.set_yticklabels([f"{c}  ({w:.02f})".replace("0.", ".")
                     for c, w in zip(cs.criterion, cs.w)], fontsize=5.8)
axb.set_xlabel("replicates scoring the TRUE species\nabove its congener (%)")
axb.set_xlim(0, 100)
axb.set_xticks([0, 25, 50, 75, 100])
axb.set_ylim(-0.8, len(cs) - 0.2)
axb.annotate("no species-level\ninformation", xy=(50, 8.4), xytext=(24, 9.7),
             fontsize=5.5, color="#444444", va="center", ha="center",
             arrowprops=dict(arrowstyle="-", lw=0.5, color="#888888",
                             shrinkA=1, shrinkB=1))
axb.set_title("Only the lighter criteria resolve species;\nthe heaviest, C12, prefers the congener",
              loc="left", fontsize=8)
axb.grid(axis="x", lw=0.4, alpha=0.35); axb.set_axisbelow(True)
axb.legend(handles=[
    Line2D([], [], marker="o", ls="", color=C_SPECIES, ms=3.6,
           label="favours true species"),
    Line2D([], [], marker="o", ls="", color=C_XR, ms=3.6,
           label="favours congener")],
    frameon=False, fontsize=5.6, loc="lower right", handletextpad=0.3,
    borderpad=0.15)
panel_letter(axb, "b")

# ============================ panel c ============================
# Two co-maintaining reservoirs: is the SET recovered, and how large is
# the Pi gap the score opens between two co-equal reservoirs?
d13 = ev[ev.scenario == "S13_dual_symmetric"]
d13b = ev[ev.scenario == "S13b_dual_asymmetric"]

# left half: recall of the true set; right half: within-truth Pi gap
xpos = [0, 1]
both = [d13.both_top2.mean() * 100, d13b.both_top2.mean() * 100]
anyt = [d13.any_top1.mean() * 100, d13b.any_top1.mean() * 100]
axc.bar([x - 0.19 for x in xpos], anyt, width=0.36, color=C_GENUS,
        edgecolor="none", zorder=2)
axc.bar([x + 0.19 for x in xpos], both, width=0.36, color=C_SPECIES,
        edgecolor="none", zorder=2)
for x, v in zip([p - 0.19 for p in xpos], anyt):
    axc.text(x, v + 1.6, f"{v:.0f}", ha="center", va="bottom", fontsize=6,
             color=C_GENUS)
for x, v in zip([p + 0.19 for p in xpos], both):
    axc.text(x, v + 1.6, f"{v:.0f}", ha="center", va="bottom", fontsize=6,
             color=C_SPECIES, fontweight="bold")
axc.set_xticks(xpos)
axc.set_xticklabels(["symmetric\n(50/50 spillover)", "asymmetric\n(60–90% / 40–10%)"],
                    fontsize=6.5)
axc.set_ylabel("% of replicates")
axc.set_ylim(0, 148); axc.set_yticks([0, 25, 50, 75, 100])
d13b_dom = (d13b.pi_h1 > d13b.pi_h2).mean() * 100
axc.text(1.42, 106, f"under asymmetric spillover the top-ranked host\n"
                    f"is the dominant reservoir in only {d13b_dom:.0f}% of replicates",
         fontsize=5.6, color="#444444", ha="right", va="bottom")
axc.set_title("With two co-maintaining reservoirs the score\nfinds one reliably, the pair less often",
              loc="left", fontsize=8)
axc.text(-0.42, 133, "top-1 is a true reservoir", fontsize=5.8, color=C_GENUS,
         va="bottom", ha="left", fontweight="bold")
axc.text(-0.42, 121, "both true reservoirs in top-2", fontsize=5.8,
         color=C_SPECIES, va="bottom", ha="left", fontweight="bold")
axc.grid(axis="y", lw=0.4, alpha=0.35); axc.set_axisbelow(True)
panel_letter(axc, "c")

# ============================ panel d ============================
# Specificity when the reservoir is unsampled. Distribution of max(Pi).
arms = [("S11c_ref_null", C_NULL, "no spillover\n(negative control)"),
        ("S14_hidden_reservoir", C_HIDDEN, "reservoir NOT\nsampled"),
        ("S1c_ref_single", C_SINGLE, "reservoir sampled\n(reference)")]
parts_x = np.arange(len(arms))
for i, (scen, col, nm) in enumerate(arms):
    v = ev[ev.scenario == scen].max_pi.dropna().values
    rng = np.random.default_rng(11 + i)
    xj = i + rng.uniform(-0.16, 0.16, size=len(v))
    axd.scatter(xj, v, s=2.2, color=col, alpha=0.28, lw=0, zorder=2)
    med = np.median(v)
    axd.plot([i - 0.30, i + 0.30], [med, med], color=col, lw=2.0, zorder=4,
             solid_capstyle="butt")
    q1, q3 = np.percentile(v, [25, 75])
    axd.plot([i, i], [q1, q3], color=col, lw=0.9, zorder=3)
    axd.text(i + 0.34, med, f"{med:.3f}", va="center", ha="left", fontsize=6,
             color=col, fontweight="bold")
axd.axhline(0.55, color="#444444", lw=0.8, ls=(0, (4, 2)), zorder=1)
axd.text(-0.42, 0.556, "call threshold  \u03a0 = 0.55", fontsize=5.8,
         color="#444444", va="bottom", ha="left")
axd.set_xticks(parts_x)
axd.set_xticklabels([a[2] for a in arms], fontsize=6.2)
axd.set_ylabel(r"max $\Pi$ across candidate hosts")
axd.set_xlim(-0.5, len(arms) - 0.22)
axd.set_title("An unsampled reservoir does not behave\nlike the negative control",
              loc="left", fontsize=8)
axd.grid(axis="y", lw=0.4, alpha=0.35); axd.set_axisbelow(True)
panel_letter(axd, "d")

n_rep = int(ev.groupby("scenario").size().min())
fig.text(0.012, 0.985,
         "The composite score resolves genus, not species, and does not "
         "detect an unsampled reservoir",
         fontsize=9.2, ha="left", va="top")
fig.text(0.012, 0.958,
         f"{n_rep} stochastic replicates per scenario; K = 5 candidate hosts "
         "(4 in the reference and unsampled-reservoir arms)",
         fontsize=6.6, ha="left", va="top", color="#555555")

fig.savefig(f"{FIG}/suppl_figure_11_taxonomic_resolution.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{FIG}/suppl_figure_11_taxonomic_resolution.pdf", bbox_inches="tight")
print("panel a  species/genus:", [f"{a:.1f}/{b:.1f}" for a, b in zip(sp, gn)])
print("panel c  any_top1:", [f"{v:.1f}" for v in anyt], " both_top2:", [f"{v:.1f}" for v in both])
print("panel d  medians:", {a[0]: round(float(np.median(ev[ev.scenario == a[0]].max_pi)), 3) for a in arms})
print("n_rep min:", n_rep)
