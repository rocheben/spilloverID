"""Tipping-point figure: rank margin of the designated host vs discount factor q
applied to the maintenance signatures, under published and reweighted schemes."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from figstyle import apply_figure_style, panel_letter, set_frame, META_GREY
apply_figure_style(sizes=(8, 7, 6))

SW = pd.read_csv("../../../data/suppl_S4_reliability_weighting/sweep_margin_vs_q.csv")
TIP = pd.read_csv("../../../data/suppl_S4_reliability_weighting/table_tipping_points.csv")

CASE_COL = {"sn": "#1b6ca8", "nv": "#c1440e", "rv": "#2e7d32"}
CASE_LAB = {"sn": "Sin Nombre\n(P. maniculatus)", "nv": "Nipah\n(Pteropus)",
            "rv": "Rift Valley fever\n(Wild ungulates)"}
LS = {"published": "-", "intermediate": (0, (5, 2)), "conservative": (0, (1.6, 1.6))}
PANEL = [("triad", "Maintenance triad C12–C14\ndown-weighted"),
         ("C12", "Inter-epidemic persistence C12\nalone down-weighted"),
         ("sero", "All serology-derived evidence\ndown-weighted")]

fig, axes = plt.subplots(1, 3, figsize=(11.4, 4.3), sharey=True)

for ax, (tgt, title) in zip(axes, PANEL):
    ax.axhline(0, color="0.25", lw=1.0, zorder=1)
    for case in ["sn", "nv", "rv"]:
        for wt in ["published", "intermediate", "conservative"]:
            s = SW[(SW.target == tgt) & (SW.case == case) & (SW.weighting == wt)].sort_values("q")
            ax.plot(s.q, s.margin, color=CASE_COL[case], ls=LS[wt],
                    lw=2.0 if wt == "published" else 1.3,
                    alpha=1.0 if wt == "published" else 0.85, zorder=3)
    # Mark every crossing; label only the published-weighting ones (the headline
    # numbers). Reweighted crossings are smaller unlabeled markers; their exact
    # q* is in out/table_tipping_points.csv. Published labels are staggered above
    # and below the zero line so that two nearby crossings cannot overprint.
    t = (TIP[(TIP.target == tgt) & TIP.q_star.notna() & (TIP.q_star < 1.0)]
         .sort_values("q_star").reset_index(drop=True))
    for _, r in t[t.weighting != "published"].iterrows():
        ax.plot([r.q_star], [0], marker="o", ms=3.6, mfc="white",
                mec=CASE_COL[r.case], mew=1.1, zorder=5, clip_on=False)
    pub = t[t.weighting == "published"].reset_index(drop=True)
    for i, r in pub.iterrows():
        ax.plot([r.q_star], [0], marker="o", ms=6.0, mfc="white",
                mec=CASE_COL[r.case], mew=1.7, zorder=5, clip_on=False)
        above = (i % 2 == 0)
        ha = "left" if r.q_star < 0.10 else ("right" if r.q_star > 0.92 else "center")
        dx = {"left": -4, "right": 4, "center": 0}[ha]
        ax.annotate(f"$q^*$={r.q_star:.2f}", (r.q_star, 0),
                    textcoords="offset points", xytext=(dx, 9 if above else -9),
                    ha=ha, va="bottom" if above else "top",
                    fontsize=6.5, color=CASE_COL[r.case], zorder=6,
                    bbox=dict(fc="white", ec="none", pad=0.9, alpha=0.9))
    ax.set_title(title, loc="left", fontsize=8)
    ax.set_xlabel("retention factor $q$")
    ax.set_xlim(-0.03, 1.03)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.margins(y=0.10)

axes[0].set_ylabel("margin of designated host\nover next-ranked host  ($\\Delta\\Pi$)")
axes[0].text(0.30, 0.145, "below 0 = designated host mis-ranked", transform=axes[0].transAxes,
             fontsize=6, color="0.3", va="bottom", ha="left")

# direct labels at q = 1 (published weighting), then decluttered vertically so
# that series whose margins are close cannot overprint each other.
series_ann = []
pts = []
for case in ["sn", "nv", "rv"]:
    s_ = SW[(SW.target == "sero") & (SW.case == case) & (SW.weighting == "published")]
    pts.append((case, float(s_[s_.q == 1.0].margin.iloc[0])))
pts.sort(key=lambda kv: kv[1], reverse=True)
for case, y in pts:
    series_ann.append(axes[2].annotate(
        CASE_LAB[case], (1.0, y), xytext=(6, 0), textcoords="offset points",
        color=CASE_COL[case], fontsize=6.5, va="center", ha="left"))


def declutter_vertical(fig, anns, pad=2.5, max_iter=400):
    """Spread offset-point annotations apart until their boxes no longer touch."""
    px_to_pt = 72.0 / fig.dpi
    for _ in range(max_iter):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        order = sorted(anns, key=lambda a: a.get_window_extent(r).y0, reverse=True)
        moved = False
        for upper, lower in zip(order, order[1:]):
            bu = upper.get_window_extent(r)
            bl = lower.get_window_extent(r)
            gap = bu.y0 - bl.y1
            if gap < pad:
                shift = (pad - gap) * px_to_pt / 2.0
                dxu, dyu = upper.xyann
                dxl, dyl = lower.xyann
                upper.xyann = (dxu, dyu + shift)
                lower.xyann = (dxl, dyl - shift)
                moved = True
        if not moved:
            return
    raise RuntimeError("declutter_vertical did not converge")

h = [mpl.lines.Line2D([], [], color="0.25", ls=LS[k], lw=2.0 if k == "published" else 1.3)
     for k in ["published", "intermediate", "conservative"]]
fig.legend(h, ["weights as reported", "reliability-weighted (intermediate)",
               "reliability-weighted (conservative)"],
           loc="lower center", bbox_to_anchor=(0.47, 0.005), ncol=3,
           frameon=False, fontsize=7, handlelength=2.8, columnspacing=2.4)

fig.text(0.47, 0.135, "$q$ = 1 keeps the evidence at full trust; $q$ = 0 discards it entirely. "
                      "The discount applied at a tipping point is $1-q^*$.",
         ha="center", fontsize=6.5, color="0.3")

fig.subplots_adjust(left=0.09, right=0.855, wspace=0.08, bottom=0.30, top=0.87)

# Labels are placed once the panel geometry is final, so the offsets reflect
# the geometry that is actually rendered.
declutter_vertical(fig, series_ann)

# --- 9.1 geometric check ---
fig.canvas.draw()
r = fig.canvas.get_renderer()
texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
         if t.get_text().strip() and t.get_visible()]
spines = [(s, s.get_window_extent(r)) for ax in fig.axes
          for s in ax.spines.values() if s.get_visible()]
tl = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
      for ax in fig.axes}
ov = [(a.get_text(), b.get_text()) for i, (a, ba) in enumerate(texts)
      for b, bb in texts[i + 1:] if ba.overlaps(bb)]
ov += [(t.get_text(), "spine") for t, bt in texts for s, bs in spines
       if bt.overlaps(bs) and t not in tl[s.axes]]
print("overlaps:", ov)
if ov:
    raise SystemExit("FAIL: %d overlapping label pair(s) remain: %s" % (len(ov), ov))

fig.savefig("../../../figures/suppl_figure_B_reliability_tipping.png", dpi=300, bbox_inches="tight")
fig.savefig("../../../figures/suppl_figure_B_reliability_tipping.pdf", bbox_inches="tight")

print("tipping points shown:", TIP[TIP.q_star.notna() & (TIP.q_star < 1)]
      [["case", "weighting", "target", "q_star"]].to_string(index=False))
