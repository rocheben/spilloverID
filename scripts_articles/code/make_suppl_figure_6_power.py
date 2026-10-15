"""Supplementary Figure 6 -- prospective power analysis.

Top-1 identification rate as a function of the number of sampled sites and the
number of animals sampled per host per site, for the canonical (a) and
bridge-host (b) scenarios.

Requires: ../data/power_grid.json (produced by power_grid.R)
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 8
rcParams['axes.labelsize'] = 8.5
rcParams['axes.titlesize'] = 9
rcParams['xtick.labelsize'] = 7.5
rcParams['ytick.labelsize'] = 7.5
rcParams['savefig.dpi'] = 300

with open('../data/power_grid.json') as f:
    pg = json.load(f)

g = pg['grid']
target = pg['target_pct']
scen_lab = {'S1_canonical': 'a  S1 canonical',
            'S2_bridge': 'b  S2 bridge host'}
scenarios = ['S1_canonical', 'S2_bridge']
sites = sorted({r['n_sites'] for r in g})
animals = sorted({r['n_animals'] for r in g})

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1),
                         gridspec_kw={'wspace': 0.30})
VERMILLION = '#D55E00'

for ax, sc in zip(axes, scenarios):
    M = np.full((len(animals), len(sites)), np.nan)
    for r in g:
        if r['scenario'] != sc:
            continue
        M[animals.index(r['n_animals']), sites.index(r['n_sites'])] = r['top1_pct']
    im = ax.imshow(M, origin='lower', cmap='YlGnBu', vmin=0, vmax=100,
                   aspect='auto')
    ax.set_xticks(range(len(sites)), [str(s) for s in sites])
    ax.set_yticks(range(len(animals)), [str(a) for a in animals])
    ax.set_xlabel('sampled sites')
    ax.set_ylabel('animals per host per site')
    ax.set_title(scen_lab[sc], loc='left', fontweight='bold')
    for i in range(len(animals)):
        for j in range(len(sites)):
            v = M[i, j]
            if np.isnan(v):
                continue
            ax.text(j, i, f'{v:.0f}', ha='center', va='center', fontsize=7,
                    color='white' if v > 55 else 'black')
            if v >= target:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                           edgecolor=VERMILLION, linewidth=2))

cb = fig.colorbar(im, ax=axes, fraction=0.030, pad=0.02)
cb.set_label('top-1 identification rate (%)', fontsize=8)
fig.suptitle(
    f"{pg['n_rep_per_cell']} replicates per cell, "
    f"{pg['n_humans_per_site']} humans per site; "
    f"cells outlined in vermillion reach the {target:.0f}% target",
    fontsize=7.5, y=1.03)

fig.savefig('../figures/suppl_figure_6_power.png', dpi=300, bbox_inches='tight')
fig.savefig('../figures/suppl_figure_6_power.pdf', bbox_inches='tight')
print('Wrote suppl_figure_6_power')
