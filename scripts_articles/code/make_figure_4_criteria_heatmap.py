"""Generate Figure 4 (per-criterion contribution heatmap) for the spilloverID Nature Microbiology submission.

Requires:
  - case_results.json (produced by run_cases.R)

Run from a directory containing case_results.json.
"""
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib import rcParams

rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 8
rcParams['axes.labelsize'] = 8.5
rcParams['axes.titlesize'] = 9
rcParams['xtick.labelsize'] = 7.5
rcParams['ytick.labelsize'] = 7.5
rcParams['legend.fontsize'] = 7.5
rcParams['savefig.dpi'] = 300
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

with open('../data/case_results.json', 'r') as f:
    data = json.load(f)

COLOR_TRUE = '#2ca02c'
COLOR_AMP = '#d62728'
COLOR_OTHER = '#a8a8a8'
from host_labels import labels as host_labels

AMPLIFIERS = {'nv': ['Pigs'], 'rv': ['Cattle'], 'sn': []}

case_meta = {
    'sn': ('a', 'Sin Nombre hantavirus',    'USA, 1993'),
    'nv': ('b', 'Nipah virus',              'Malaysia, 1998-99'),
    'rv': ('c', 'Rift Valley fever',        'Kenya, 2006-07')
}

# ===================== FIGURE 3: Per-criterion heatmap =====================
crit_names = ['C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14']

fig, axes = plt.subplots(3, 1, figsize=(8.5, 7.5),
                          gridspec_kw={'hspace': 0.65})
cmap = plt.cm.RdYlGn

for ax, key in zip(axes, ['sn', 'nv', 'rv']):
    c = data[key]
    hosts = c['hosts']
    crit_df = c['criteria']
    
    mat = np.full((len(hosts), len(crit_names)), np.nan)
    for i, row in enumerate(crit_df):
        for j, cn in enumerate(crit_names):
            v = row.get(cn)
            mat[i, j] = float(v) if v is not None else np.nan
    
    masked = np.ma.masked_invalid(mat)
    cmap_with_bad = plt.cm.RdYlGn.copy()
    cmap_with_bad.set_bad('#dddddd')
    
    im = ax.imshow(masked, aspect='auto', cmap=cmap_with_bad, vmin=0, vmax=1)
    
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            if np.isnan(v):
                ax.text(j, i, 'NA', ha='center', va='center', fontsize=7, color='#666666')
            else:
                txt_color = 'white' if (v > 0.70 or v < 0.20) else 'black'
                ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=7.2,
                        color=txt_color)
    
    ax.set_xticks(range(len(crit_names)))
    ax.set_xticklabels(crit_names, fontsize=8)
    ax.set_yticks(range(len(hosts)))
    ax.set_yticklabels(host_labels(hosts), fontsize=8.5)
    
    true_host = c['true_host']
    for tick, lab in zip(ax.get_yticklabels(), hosts):
        if lab == true_host:
            tick.set_fontweight('bold')
            tick.set_color(COLOR_TRUE)
    
    panel, name, when = case_meta[key]
    ax.set_title(f'{panel}   {name} ({when})', fontsize=9, fontweight='bold', loc='left')

# Bottom colorbar
cbar_ax = fig.add_axes([0.30, 0.02, 0.40, 0.018])
sm = plt.cm.ScalarMappable(cmap=plt.cm.RdYlGn, norm=plt.Normalize(0, 1))
sm.set_array([])
cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal')
cbar.set_label('Criterion score $s_{k,c}$ (normalized to [0, 1])', fontsize=8.5)

plt.savefig('../figures/figure_4_criteria_heatmap.png', dpi=300, bbox_inches='tight')
plt.savefig('../figures/figure_4_criteria_heatmap.pdf', bbox_inches='tight')
plt.close()
print("Figure 4 (criteria heatmap) saved.")


