"""Supplementary Figure 7 -- bootstrap convergence.

Width of the 95% site-level bootstrap percentile interval on the true host's
composite score as the number of bootstrap replicates B increases, for the three
reconstructed outbreaks.

Requires: ../data/boot_converge.json (produced by boot_converge.R)
"""
import json
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 8
rcParams['axes.labelsize'] = 8.5
rcParams['xtick.labelsize'] = 7.5
rcParams['ytick.labelsize'] = 7.5
rcParams['savefig.dpi'] = 300
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

with open('../data/boot_converge.json') as f:
    bc = json.load(f)

COL = {'sn': '#009E73', 'nv': '#0072B2', 'rv': '#D55E00'}
LAB = {'sn': 'Sin Nombre (canonical)', 'nv': 'Nipah (bridge host)',
       'rv': 'Rift Valley fever (vector-borne)'}

fig, ax = plt.subplots(figsize=(4.6, 3.2))
for k in ('sn', 'nv', 'rv'):
    z = bc[k]
    ax.plot(z['B'], z['width'], 'o-', color=COL[k], markersize=4,
            linewidth=1.4, label=LAB[k])
    ax.axhline(z['width_at_Bmax'], color=COL[k], linestyle=':', linewidth=0.8,
               alpha=0.6)

ax.axvline(100, color='0.4', linestyle='--', linewidth=1)
ax.text(105, ax.get_ylim()[1] * 1.10, 'B = 100 used in the main analyses',
        fontsize=7, va='top', color='0.3')
ax.set_xlabel('bootstrap replicates B')
ax.set_ylabel("width of the 95% interval on $\\Pi$\n(true host)", fontsize=9)
ax.set_ylim(bottom=0, top=ax.get_ylim()[1] * 1.12)
ax.legend(frameon=False, fontsize=7, loc='lower right')

fig.savefig('../figures/suppl_figure_7_bootconv.png', dpi=300, bbox_inches='tight')
fig.savefig('../figures/suppl_figure_7_bootconv.pdf', bbox_inches='tight')
print('Wrote suppl_figure_7_bootconv')
