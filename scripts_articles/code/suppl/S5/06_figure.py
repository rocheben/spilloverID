"""
Supplementary Figure 9, two panels:
  (a) reliability curve -- observed accuracy against the observable margin;
  (b) the no-spillover null under the two availability regimes that matter,
      the simulated benchmark (13 active criteria) and the three outbreak
      reconstructions (11 active: C6, C8 and C11 unavailable, C9 present but
      returning zero), with the nominal and the recalibrated thresholds.
"""
import numpy as np, pandas as pd, matplotlib as mpl
from pathlib import Path
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from figstyle import (apply_figure_style, panel_letter, set_frame,
                      panel_crops, META_GREY)

OUT = Path("../../../data/suppl_S5_uncertainty")
apply_figure_style()

rec = pd.read_csv(OUT / "benchmark_recomputed_margins.csv")
s11 = pd.read_csv(OUT / "regime_per_rep.csv")
thr = pd.read_csv(OUT / "threshold_recalibration.csv")
ce  = pd.read_csv(OUT / "confident_error_rates.csv")

pos = rec[rec.k_star.notna()].copy(); pos["top1"] = pos["top1"].astype(int)

C_OBS = "#1f4e79"   # observable margin
C_B13 = "#4c72b0"   # 13 discriminating criteria (calibration regime)
C_R11 = "#c44e52"   # 11 active criteria (reconstruction regime)
REGIME_APP = "reconstruction11"   # the availability pattern of the three cases

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.3))

# ---------------- panel a: reliability curve ---------------------------------
ax = axes[0]
edges = np.unique(np.quantile(pos.margin_obs, np.linspace(0, 1, 11)))
b = pd.cut(pos.margin_obs, edges, include_lowest=True)
g = pos.groupby(b, observed=True).agg(n=("top1", "size"), acc=("top1", "mean"),
                                      x=("margin_obs", "median")).reset_index(drop=True)
z = 1.959964; n, p = g.n, g.acc
den = 1 + z**2/n; ctr = (p + z**2/(2*n))/den
hw = z*np.sqrt(p*(1-p)/n + z**2/(4*n**2))/den
g["lo"], g["hi"] = (ctr-hw).clip(0, 1), (ctr+hw).clip(0, 1)

ax.axvline(0.15, color=META_GREY, lw=0.9, ls="--", zorder=1)
ax.errorbar(g.x, g.acc, yerr=[g.acc-g.lo, g.hi-g.acc], fmt="o-", color=C_OBS,
            ms=4.5, lw=1.6, capsize=2, zorder=3)
ax.text(0.1445, 0.345, "0.15 threshold\n(Discussion)", fontsize=6,
        color=META_GREY, ha="right", va="bottom")

row = ce[(ce.margin_def.str.startswith("observable")) & (ce.threshold == 0.15)].iloc[0]
ax.text(0.163, 0.545,
        f"margin $\\geq$ 0.15:  {row.accuracy_given_above*100:.1f}% correct\n"
        f"    ({int(row.n_confident_errors)} confident errors / {int(row.n_above_threshold)} calls)",
        fontsize=6.5, color=C_OBS, ha="left", va="center")
ax.text(0.163, 0.435,
        f"margin < 0.15:  {row.accuracy_given_below*100:.1f}% correct",
        fontsize=6.5, color=META_GREY, ha="left", va="center")

ax.set_xlabel("observable margin  $\\Pi_{(1)}-\\Pi_{(2)}$", labelpad=6)
ax.set_ylabel("P(true reservoir ranked first)", labelpad=4)
ax.set_title("Wrong calls come with small margins", loc="left")
ax.set_ylim(0.33, 1.055); ax.set_xlim(-0.014, 0.368)
set_frame(ax)
panel_letter(ax, "a")

# ---------------- panel b: negative-control null ------------------------------
ax = axes[1]
bins = np.linspace(0.25, 0.66, 42)
v13 = s11.loc[(s11.scenario == "S11_no_spillover") & (s11.regime == "bench13"), "max_pi"]
v11 = s11.loc[(s11.scenario == "S11_no_spillover") & (s11.regime == REGIME_APP), "max_pi"]
ax.hist(v13, bins=bins, color=C_B13, alpha=0.40, zorder=2)
ax.hist(v13, bins=bins, histtype="step", color=C_B13, lw=1.3, zorder=3)
ax.hist(v11, bins=bins, color=C_R11, alpha=0.40, zorder=2)
ax.hist(v11, bins=bins, histtype="step", color=C_R11, lw=1.3, zorder=3)

t_rec = float(thr.loc[(thr.regime == REGIME_APP) &
                      (thr.published_threshold == 0.55), "equivalent_threshold"].iloc[0])
ax.axvline(0.55,  color=C_B13, lw=1.5, zorder=5)
ax.axvline(t_rec, color=C_R11, lw=1.5, ls="--", zorder=5)
ymax = ax.get_ylim()[1]
ax.set_ylim(0, ymax * 1.22)
ymax = ax.get_ylim()[1]
ax.text(0.55 - 0.006, ymax*0.88, "0.55\nnominal", color=C_B13, fontsize=6.5,
        ha="right", va="top")
ax.text(t_rec + 0.006, ymax*0.88, f"{t_rec:.3f}\nrecalibrated", color=C_R11,
        fontsize=6.5, ha="left", va="top")
ax.text(0.268, ymax*0.96, "13 active criteria\n(simulated benchmark)",
        color=C_B13, fontsize=6.5, ha="left", va="top")
ax.text(0.268, ymax*0.76, "11 active criteria\n(reconstructions; C6, C8, C11 NA, C9 = 0)",
        color=C_R11, fontsize=6.5, ha="left", va="top")
ax.set_xlabel("max $\\Pi$ under the no-spillover control (S11)", labelpad=6)
ax.set_ylabel("replicates (of 300)", labelpad=4)
ax.set_title("The null shifts with the criteria available", loc="left")
ax.set_xlim(0.25, 0.66)
set_frame(ax)
panel_letter(ax, "b")

fig.tight_layout(w_pad=2.0)
fig.savefig("../../../figures/suppl_figure_9_uncertainty_calibration.png", dpi=300, bbox_inches="tight")
fig.savefig("../../../figures/suppl_figure_9_uncertainty_calibration.pdf", bbox_inches="tight")

# ---- 9.1 geometric check ----
r = fig.canvas.get_renderer()
texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
         if t.get_text().strip() and t.get_visible()]
spines = [(s, s.get_window_extent(r)) for a in fig.axes
          for s in a.spines.values() if s.get_visible()]
tl = {a: set(a.get_xticklabels(which="both") + a.get_yticklabels(which="both"))
      for a in fig.axes}
ov = [(a.get_text()[:28], b.get_text()[:28]) for i, (a, ba) in enumerate(texts)
      for b, bb in texts[i+1:] if ba.overlaps(bb)]
ov += [(t.get_text()[:28], "spine") for t, bt in texts for s, bs in spines
       if bt.overlaps(bs) and t not in tl[s.axes]]
print("OVERLAPS:", ov)
print("null medians  13 active: %.4f   11 active: %.4f" % (v13.median(), v11.median()))
print("recalibrated threshold: %.4f" % t_rec)
print(g.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
