"""Generate Supplementary Figure 3 (per-replicate benchmark scatter) for the spilloverID Nature Microbiology submission.

Requires:
  - benchmark_results.json (produced by run_benchmark.R)
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

with open('../data/benchmark_results.json') as f:
    bench = json.load(f)

case_keys = ['sn', 'nv', 'rv']

# =====================================================================
# SUPPL FIG 3: Per-replicate scatter plot benchmark
# Show every replicate's true-host Pi vs max-other-Pi, colored by scenario
# =====================================================================
per_rep = bench['per_rep']
fig, ax = plt.subplots(figsize=(6.5, 5.5))

scen_order = ['S1_canonical', 'S2_bridge', 'S3_cryptic', 'S4_confounded',
              'S5_codistributed', 'S6_strong_seasonal', 'S7_low_power',
              'S8_measurement_error', 'S9_weak_signal', 'S10_short_immunity']

cmap = plt.cm.tab10
for i, scen in enumerate(scen_order):
    rows = [r for r in per_rep if r['scenario'] == scen]
    pi_true = [r['pi_true'] for r in rows]
    pi_max_other = [r['pi_max_other'] for r in rows]
    label = scen.replace('_', ' ').replace('S', 'S').lstrip()
    ax.scatter(pi_max_other, pi_true, color=cmap(i), s=35, alpha=0.85,
               edgecolor='black', linewidth=0.4, label=label)

# Add S11 negative control as black X
s11 = [r for r in per_rep if r['scenario'] == 'S11_no_spillover']
s11_max_pi = [r['max_pi'] for r in s11]
ax.scatter(s11_max_pi, s11_max_pi, color='black', s=50, marker='x',
           label='S11 no spillover (max Π on identity line)', linewidth=2)

# Diagonal line: above = truth wins
ax.plot([0, 1], [0, 1], color='gray', linestyle='--', linewidth=1, alpha=0.6)

# 0.55 threshold lines
ax.axhline(0.55, color='#D55E00', linestyle=':', alpha=0.7, linewidth=1)
ax.axvline(0.55, color='#D55E00', linestyle=':', alpha=0.7, linewidth=1)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel(r'$\max_{k\neq \mathrm{truth}} \Pi_k$')
ax.set_ylabel(r'$\Pi_{\mathrm{truth}}$')
ax.set_title('Per-replicate benchmark scatter: positive scenarios (10 colors)\nand S11 negative control (black ×)',
             fontsize=9, fontweight='bold', loc='left')
ax.legend(loc='upper left', fontsize=6.5, frameon=False, ncol=2,
          bbox_to_anchor=(1.02, 1.0))
ax.text(0.95, 0.55, '\u03A0 = 0.55 threshold', fontsize=7, color='#D55E00',
        ha='right', va='bottom')

plt.savefig('../figures/suppl_figure_3_per_rep_scatter.png',
            dpi=300, bbox_inches='tight')
plt.savefig('../figures/suppl_figure_3_per_rep_scatter.pdf', bbox_inches='tight')
plt.close()
print("Supplementary Figure 3 saved.")


