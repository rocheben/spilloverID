"""Generate Supplementary Figure 2 (leave-one-criterion-out heatmap) for the spilloverID Nature Microbiology submission.

Requires:
  - suppl_loco.json (produced by run_suppl.R)
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

with open('../data/suppl_loco.json') as f:
    loco = json.load(f)

cases_short = {'sn': 'Sin Nombre', 'nv': 'Nipah', 'rv': 'RVF'}
case_color = {'sn': '#009E73', 'nv': '#0072B2', 'rv': '#D55E00'}

case_keys = ['sn', 'nv', 'rv']

# =====================================================================
# SUPPL FIG 2: Leave-one-criterion-out heatmap
# =====================================================================
crit_names = [f'C{i}' for i in range(1, 15)]
fig, ax = plt.subplots(figsize=(7.5, 3.3))

# Build delta-Pi matrix: rows = cases, cols = criteria
delta_matrix = np.zeros((3, 14))
rank_matrix = np.zeros((3, 14), dtype=int)
for ci, ck in enumerate(case_keys):
    rows = loco[ck]['rows']
    for c in rows:
        ji = crit_names.index(c['criterion'])
        delta_matrix[ci, ji] = c['delta_pi']
        rank_matrix[ci, ji] = c['rank_true']

# Colormap: blue (negative dPi = removing this hurts truth = good criterion) → red (positive)
cmap = plt.cm.RdBu_r
vmax = max(abs(delta_matrix.min()), abs(delta_matrix.max()))
im = ax.imshow(delta_matrix, aspect='auto', cmap=cmap, vmin=-vmax, vmax=vmax)

# Annotate cells: dPi + flag if rank dropped
for ci in range(3):
    for ji in range(14):
        v = delta_matrix[ci, ji]
        r = rank_matrix[ci, ji]
        txt = f'{v:+.02f}'
        if r > 1:
            txt += f'\n\u2192 #{r}'
            txt_color = 'white' if abs(v) > 0.05 else 'black'
            ax.text(ji, ci, txt, ha='center', va='center', fontsize=6.5,
                    color=txt_color, fontweight='bold')
        else:
            txt_color = 'white' if abs(v) > 0.06 else 'black'
            ax.text(ji, ci, txt, ha='center', va='center', fontsize=6.5,
                    color=txt_color)

ax.set_xticks(range(14))
ax.set_xticklabels(crit_names, fontsize=8)
ax.set_yticks(range(3))
ax.set_yticklabels(['Sin Nombre\n(canonical)',
                    'Nipah\n(bridge host)',
                    'RVF\n(vector-borne)'], fontsize=8)
ax.set_title('Leave-one-criterion-out: change in $\\Pi$ of true reservoir, with rank flag if not #1',
             fontsize=9, fontweight='bold', loc='left')
cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
cbar.set_label(r'$\Delta\Pi$ of true host when criterion is removed', fontsize=8)

plt.savefig('../figures/suppl_figure_2_loco.png',
            dpi=300, bbox_inches='tight')
plt.savefig('../figures/suppl_figure_2_loco.pdf', bbox_inches='tight')
plt.close()
print("Supplementary Figure 2 saved.")


