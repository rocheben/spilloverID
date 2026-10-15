"""Generate Supplementary Figure 1 (weight perturbation sensitivity) for the spilloverID Nature Microbiology submission.

Requires:
  - suppl_sensitivity.json (produced by run_suppl.R)
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.patches as mpatches

rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 8
rcParams['axes.labelsize'] = 8.5
rcParams['axes.titlesize'] = 9
rcParams['xtick.labelsize'] = 7.5
rcParams['ytick.labelsize'] = 7.5
rcParams['savefig.dpi'] = 300
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

with open('../data/suppl_sensitivity.json') as f:
    sens = json.load(f)

cases_short = {'sn': 'Sin Nombre', 'nv': 'Nipah', 'rv': 'RVF'}
case_color = {'sn': '#009E73', 'nv': '#0072B2', 'rv': '#D55E00'}

# =====================================================================
# SUPPL FIG 1: Weight perturbation sensitivity
# =====================================================================
fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.2),
                         gridspec_kw={'width_ratios': [1, 1.2], 'wspace': 0.4})

# Panel a: rank-1 stability
ax = axes[0]
case_keys = ['sn', 'nv', 'rv']
case_labels = ['Sin Nombre\n(canonical)', 'Nipah\n(bridge host)', 'RVF\n(vector-borne)']
pct = [sens[k]['pct_rank1'] for k in case_keys]
colors = [case_color[k] for k in case_keys]
bars = ax.bar(range(3), pct, color=colors, edgecolor='black', linewidth=0.8, width=0.6)
for i, v in enumerate(pct):
    ax.text(i, v + 1.5, f'{v:.1f}%', ha='center', fontsize=9, fontweight='bold')
ax.axhline(25, color='#000000', linestyle=':', linewidth=1.0, alpha=0.55)
ax.text(2.4, 26, 'chance (1/4)', fontsize=7, color='#000000', ha='right')
ax.set_xticks(range(3))
ax.set_xticklabels(case_labels, fontsize=8)
ax.set_ylabel('Rank-1 stability under\nweight perturbation (%)')
ax.set_ylim(0, 115)
ax.set_title('a   Rank stability under random weight perturbation\n(n = 200 Dirichlet draws per case)',
             loc='left', fontsize=9, fontweight='bold')

# Panel b: Pi distributions
ax = axes[1]
# Per-case distribution from quantiles (approximate boxplot)
for i, k in enumerate(case_keys):
    s = sens[k]
    baseline = s['baseline_pi']
    med = s['pi_median']; lo = s['pi_2.5']; hi = s['pi_97.5']
    ax.plot([i, i], [lo, hi], color=case_color[k], linewidth=2.5)
    ax.plot([i-0.18, i+0.18], [lo, lo], color=case_color[k], linewidth=2)
    ax.plot([i-0.18, i+0.18], [hi, hi], color=case_color[k], linewidth=2)
    ax.plot([i-0.22, i+0.22], [med, med], color=case_color[k], linewidth=3)
    ax.scatter([i], [baseline], color='#E69F00', s=70, zorder=5, marker='D',
               edgecolor='black', linewidth=0.8,
               label='baseline (default weights)' if i == 0 else None)

# 0.55 threshold
ax.axhline(0.55, color='#888888', linestyle='--', linewidth=1.0, alpha=0.7)
ax.text(2.4, 0.555, 'threshold \u03A0 \u2248 0.55', fontsize=7, color='#666666', ha='right')

ax.set_xticks(range(3))
ax.set_xticklabels(case_labels, fontsize=8)
ax.set_ylabel(r'$\Pi_{\mathrm{true\ host}}$ under weight perturbation')
ax.set_ylim(0.30, 0.95)
ax.set_title(r'b   $\Pi$ distribution of true host under weight perturbation',
             loc='left', fontsize=9, fontweight='bold')
ax.legend(loc='lower right', fontsize=7, frameon=False)

plt.savefig('../figures/suppl_figure_1_sensitivity.png',
            dpi=300, bbox_inches='tight')
plt.savefig('../figures/suppl_figure_1_sensitivity.pdf', bbox_inches='tight')
plt.close()
print("Supplementary Figure 1 saved.")


