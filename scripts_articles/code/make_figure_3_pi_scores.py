"""Generate Figure 3 (composite plausibility scores across three real outbreaks)."""
import json
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

# Try several likely locations for case_results.json and give a helpful error
candidates = ['../data/case_results.json']
data = None
for p in candidates:
    if os.path.exists(p):
        with open(p, 'r') as f:
            data = json.load(f)
        break
if data is None:
    raise FileNotFoundError(
        "case_results.json not found. Run `run_cases.R` from the scripts directory "
        "or place the file in the current directory (or data/). Searched: " + ", ".join(candidates)
    )

from host_labels import labels as host_labels

COLOR_TRUE = '#009E73'
COLOR_AMP = '#D55E00'
COLOR_OTHER = '#999999'
AMPLIFIERS = {'nv': ['Pigs'], 'rv': ['Cattle'], 'sn': []}

case_meta = {
    'sn': ('a', 'Sin Nombre hantavirus',    'USA, 1993'),
    'nv': ('b', 'Nipah virus',              'Malaysia, 1998-99'),
    'rv': ('c', 'Rift Valley fever',        'Kenya, 2006-07')
}

# ===================== FIGURE 2: Π scores =====================
fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.4), gridspec_kw={'wspace': 1.05})

for pos_i, (ax, key) in enumerate(zip(axes, ['sn', 'nv', 'rv'])):
    c = data[key]
    hosts = c['hosts']
    pi_vals = c['pi']
    true_host = c['true_host']
    ci_low = c.get('ci_low') or [None]*len(hosts)
    ci_high = c.get('ci_high') or [None]*len(hosts)
    
    order = list(np.argsort(pi_vals)[::-1])
    hosts_ord = [hosts[i] for i in order]
    pi_ord = [pi_vals[i] for i in order]
    cilo_ord = [ci_low[i] for i in order]
    cihi_ord = [ci_high[i] for i in order]
    
    colors = []
    for h_ in hosts_ord:
        if h_ == true_host: colors.append(COLOR_TRUE)
        elif h_ in AMPLIFIERS[key]: colors.append(COLOR_AMP)
        else: colors.append(COLOR_OTHER)
    
    y_pos = np.arange(len(hosts_ord))
    ax.barh(y_pos, pi_ord, color=colors, edgecolor='black', linewidth=0.6, height=0.7)
    
    for i, (lo, hi) in enumerate(zip(cilo_ord, cihi_ord)):
        if lo is not None and hi is not None:
            ax.plot([lo, hi], [i, i], 'k-', linewidth=1.0)
            ax.plot([lo, lo], [i-0.10, i+0.10], 'k-', linewidth=1.0)
            ax.plot([hi, hi], [i-0.10, i+0.10], 'k-', linewidth=1.0)
    
    for i, v in enumerate(pi_ord):
        offset = max(0.025, (cihi_ord[i] or v) - v + 0.025)
        ax.text(v + offset, i, f'{v:.2f}', va='center', fontsize=7.5)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(host_labels(hosts_ord), fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.set_xlabel(r'Composite plausibility score $\Pi_k$')
    
    panel = ['a','b','c'][pos_i]
    _, name, when = case_meta[key]
    ax.set_title(f'{name}\n{when}', fontsize=8.8, fontweight='bold')
    ax.text(-0.50, 1.13, panel, transform=ax.transAxes, fontsize=11, fontweight='bold')
    
    # Star marker (inside the reservoir bar, so it never collides with the host labels)
    truth_idx_ord = hosts_ord.index(true_host)
    ax.annotate('★', xy=(0.045, truth_idx_ord), xytext=(0.045, truth_idx_ord),
                xycoords='data', textcoords='data',
                fontsize=11, ha='center', va='center', color='white', zorder=5)

legend_handles = [
    mpatches.Patch(facecolor=COLOR_TRUE, edgecolor='black',
                   label='Maintenance reservoir (★, molecular evidence)'),
    mpatches.Patch(facecolor=COLOR_AMP, edgecolor='black', label='Amplifier / bridge host'),
    mpatches.Patch(facecolor=COLOR_OTHER, edgecolor='black', label='Other candidate host')
]
fig.legend(handles=legend_handles, loc='lower center', ncol=3, frameon=False,
           bbox_to_anchor=(0.5, -0.06), fontsize=8)
plt.savefig('../figures/figure_3_pi_scores.png', dpi=300, bbox_inches='tight')
plt.savefig('../figures/figure_3_pi_scores.pdf', bbox_inches='tight')
plt.close()
print("Figure 3 (pi scores) saved.")


