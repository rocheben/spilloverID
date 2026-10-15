"""
Confidence/accuracy calibration of the composite score on the 11-scenario
benchmark.

When the framework ranks the wrong host first, is the margin
that accompanies the call small (honest uncertainty) or large (confident error)?

Two margin definitions are contrasted:
  margin_oracle = Pi(true host) - max Pi(other hosts)   [deposited `margin`;
                  requires knowing the truth, so unavailable in application]
  margin_obs    = Pi_(1) - Pi_(2)                        [observable gap]
Only margin_obs can serve as a confidence measure in a real analysis.
"""
import json, numpy as np, pandas as pd
from pathlib import Path

OUT = Path("../../../data/suppl_S5_uncertainty")
THR_DISCUSSION = 0.15          # margin threshold recommended in the Discussion

# ---- per-replicate benchmark table -----------------------------------------
dep = pd.DataFrame(json.load(open("../../../data/benchmark_results.json"))["per_rep"])
dep = dep.rename(columns={"margin": "margin_oracle"})
pos_dep = dep[dep["k_star"].notna()].copy()
pos_dep["top1"] = pos_dep["top1"].astype(int)

# ---- recomputed table with the observable margin ---------------------------
rec_path = OUT / "benchmark_recomputed_margins.csv"
rec = pd.read_csv(rec_path) if rec_path.exists() else None

def reliability(df, margin_col, bins):
    """Observed accuracy within bins of the confidence measure."""
    b = pd.cut(df[margin_col], bins=bins, include_lowest=True)
    g = df.groupby(b, observed=True).agg(
        n=("top1", "size"), n_correct=("top1", "sum"),
        margin_mean=(margin_col, "mean"), margin_med=(margin_col, "median"))
    g["accuracy"] = g["n_correct"] / g["n"]
    # Wilson 95% interval on the accuracy
    z, n, p = 1.959964, g["n"], g["accuracy"]
    den = 1 + z**2 / n
    ctr = (p + z**2 / (2 * n)) / den
    hw = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / den
    g["acc_lo"], g["acc_hi"] = (ctr - hw).clip(0, 1), (ctr + hw).clip(0, 1)
    return g.reset_index().rename(columns={margin_col: "bin"})

rows_rel, rows_conf = [], []

def analyse(df, label, margin_col):
    d = df.dropna(subset=[margin_col, "top1"]).copy()
    d["top1"] = d["top1"].astype(int)
    lo, hi = d[margin_col].min(), d[margin_col].max()
    edges = np.unique(np.round(np.quantile(d[margin_col], np.linspace(0, 1, 11)), 6))
    rel = reliability(d, margin_col, edges)
    rel.insert(0, "margin_def", label); rel["bin"] = rel["bin"].astype(str)
    rows_rel.append(rel)

    # confident-error rate: wrong top-1 while the margin exceeds the threshold
    for thr in (0.10, THR_DISCUSSION, 0.20):
        above = d[margin_col] >= thr
        n_above = int(above.sum())
        n_conf_err = int(((d["top1"] == 0) & above).sum())
        rows_conf.append(dict(
            margin_def=label, threshold=thr, n_total=len(d),
            n_above_threshold=n_above,
            frac_above_threshold=n_above / len(d),
            n_confident_errors=n_conf_err,
            confident_error_rate_overall=n_conf_err / len(d),
            error_rate_given_above=n_conf_err / n_above if n_above else np.nan,
            accuracy_given_above=1 - n_conf_err / n_above if n_above else np.nan,
            accuracy_given_below=(d.loc[~above, "top1"].mean()
                                  if (~above).sum() else np.nan),
            n_below=int((~above).sum())))
    return d

d_or = analyse(pos_dep, "oracle (Pi_true - max Pi_other)", "margin_oracle")
if rec is not None:
    pos_rec = rec[rec["k_star"].notna()].copy()
    d_ob = analyse(pos_rec, "observable (Pi_1 - Pi_2)", "margin_obs")
    d_or2 = analyse(pos_rec, "oracle (recomputed)", "margin_oracle")

rel_all = pd.concat(rows_rel, ignore_index=True)
conf_all = pd.DataFrame(rows_conf)
rel_all.to_csv(OUT / "reliability_bins.csv", index=False)
conf_all.to_csv(OUT / "confident_error_rates.csv", index=False)

# ---- per-scenario breakdown -------------------------------------------------
def per_scen(d, margin_col, label):
    g = d.groupby("scenario").apply(lambda x: pd.Series({
        "n": len(x), "top1_rate": x["top1"].mean(),
        "margin_mean": x[margin_col].mean(),
        "margin_mean_correct": x.loc[x.top1 == 1, margin_col].mean(),
        "margin_mean_error": x.loc[x.top1 == 0, margin_col].mean(),
        "frac_errors_above_0.15": ((x.top1 == 0) & (x[margin_col] >= THR_DISCUSSION)).sum() /
                                  max((x.top1 == 0).sum(), 1),
        "confident_error_rate": ((x.top1 == 0) & (x[margin_col] >= THR_DISCUSSION)).mean(),
    }), include_groups=False).reset_index()
    g.insert(0, "margin_def", label)
    return g

scen_rows = [per_scen(d_or, "margin_oracle", "oracle (Pi_true - max Pi_other)")]
if rec is not None:
    scen_rows.append(per_scen(d_ob, "margin_obs", "observable"))
scen = pd.concat(scen_rows, ignore_index=True)
scen.to_csv(OUT / "calibration_by_scenario.csv", index=False)

# ---- discrimination: AUC of the margin for predicting a correct call --------
def auc(y, s):
    y = np.asarray(y); s = np.asarray(s)
    r = pd.Series(s).rank().values
    n1, n0 = (y == 1).sum(), (y == 0).sum()
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

auc_rows = [dict(margin_def="oracle (Pi_true - max Pi_other)", auc=auc(d_or.top1, d_or.margin_oracle),
                 n=len(d_or))]
if rec is not None:
    auc_rows.append(dict(margin_def="observable", auc=auc(d_ob.top1, d_ob.margin_obs),
                         n=len(d_ob)))
aucs = pd.DataFrame(auc_rows)
aucs.to_csv(OUT / "margin_auc.csv", index=False)

print("== confident-error rates ==")
print(conf_all.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
print("\n== per-scenario (observable margin where available) ==")
print(scen.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print("\n== AUC ==")
print(aucs.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
