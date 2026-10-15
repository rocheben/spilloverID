"""Generate Figure 5 (discrimination margin and C12 signal) for the spilloverID Nature Microbiology submission.

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

COLOR_TRUE = '#009E73'
COLOR_AMP = '#D55E00'
COLOR_OTHER = '#999999'
AMPLIFIERS = {'nv': ['Pigs'], 'rv': ['Cattle'], 'sn': []}

case_meta = {
    'sn': ('a', 'Sin Nombre hantavirus',    'USA, 1993'),
    'nv': ('b', 'Nipah virus',              'Malaysia, 1998-99'),
    'rv': ('c', 'Rift Valley fever',        'Kenya, 2006-07')
}

# ===================== FIGURE 4: Margins + C12 signal =====================
fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.9), gridspec_kw={'wspace': 0.40})
plt.subplots_adjust(bottom=0.26, top=0.86)

case_keys = ['sn', 'nv', 'rv']
case_short = ['Sin Nombre\n(canonical)', 'Nipah\n(bridge host)', 'RVF\n(vector-borne)']

# --- Panel a: margin ---
ax = axes[0]
margins = []
ci_low_margin = []
ci_high_margin = []
for key in case_keys:
    c = data[key]
    pi_vals = np.array(c['pi'])
    hosts = c['hosts']
    truth_idx = hosts.index(c['true_host'])
    others = [i for i in range(len(hosts)) if i != truth_idx]
    truth_pi = pi_vals[truth_idx]
    best_other = max(pi_vals[others])
    best_other_idx = others[int(np.argmax(pi_vals[others]))]
    margins.append(c.get('margin', truth_pi - best_other))
    # Percentile interval of the DIFFERENCE across bootstrap replicates, written
    # by run_cases.R. Differencing the two marginal intervals instead would
    # ignore the correlation between hosts within a replicate and overstate the
    # width by a wide margin.
    if c.get('margin_ci_low') is not None:
        ci_low_margin.append(c['margin_ci_low'])
        ci_high_margin.append(c['margin_ci_high'])
    else:
        raise SystemExit('case_results.json lacks margin_ci_low / margin_ci_high; '
                         're-run run_cases.R')

x = np.arange(len(case_keys))
ax.bar(x, margins, color=COLOR_TRUE, edgecolor='black', linewidth=0.8, width=0.55)
ax.errorbar(x, margins,
            yerr=[[margins[i]-ci_low_margin[i] for i in range(3)],
                  [ci_high_margin[i]-margins[i] for i in range(3)]],
            fmt='none', ecolor='black', capsize=4, linewidth=1.0)
for i, v in enumerate(margins):
    ax.text(i, max(v + 0.04, 0.05), f'{v:+.2f}', ha='center', fontsize=8.5, fontweight='bold')

ax.axhline(0, color='#D55E00', linestyle='--', linewidth=1.0, alpha=0.6)
ax.set_xticks(x)
ax.set_xticklabels(case_short, fontsize=8.5)
ax.set_ylabel(r'$\Pi_{\mathrm{truth}} - \max \Pi_{\mathrm{others}}$')
ax.set_ylim(-0.10, 0.65)
ax.set_title('Discrimination margin', fontsize=9, fontweight='bold')
ax.text(-0.18, 1.10, 'a', transform=ax.transAxes, fontsize=11, fontweight='bold')

# --- Panel b: C12 raw signal ---
ax = axes[1]
width = 0.32
for i, key in enumerate(case_keys):
    c = data[key]
    hosts = c['hosts']
    crit_rows = c['criteria']
    c12 = [r.get('C12') for r in crit_rows]
    truth_idx = hosts.index(c['true_host'])
    amps = AMPLIFIERS[key]
    if amps:
        amp_idx = hosts.index(amps[0])
        amp_color = COLOR_AMP
    else:
        order = np.argsort(c['pi'])[::-1]
        amp_idx = order[1] if order[0] == truth_idx else order[0]
        amp_color = COLOR_OTHER
    
    truth_v = c12[truth_idx] if c12[truth_idx] is not None else 0
    amp_v = c12[amp_idx] if c12[amp_idx] is not None else 0
    
    ax.bar(i - width/2, truth_v, width, color=COLOR_TRUE, edgecolor='black', linewidth=0.7,
           label='Maintenance reservoir' if i == 0 else None)
    ax.bar(i + width/2, amp_v, width, color=amp_color, edgecolor='black', linewidth=0.7,
           label='Amplifier / 2nd-ranked' if i == 0 else None)
    ax.text(i - width/2, truth_v + 0.02, f'{truth_v:.2f}', ha='center', fontsize=7.5)
    ax.text(i + width/2, amp_v + 0.02, f'{amp_v:.2f}', ha='center', fontsize=7.5)

ax.set_xticks(range(len(case_keys)))
ax.set_xticklabels(case_short, fontsize=8.5)
ax.set_ylim(0, 1.18)
ax.set_ylabel('Inter-epidemic persistence\n(criterion C12)')
ax.set_title('Maintenance signature separates\nreservoir from amplifier', fontsize=9, fontweight='bold')
ax.text(-0.18, 1.10, 'b', transform=ax.transAxes, fontsize=11, fontweight='bold')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=7.5, frameon=False)

plt.savefig('../figures/figure_5_margins_and_C12.png', dpi=300, bbox_inches='tight')
plt.savefig('../figures/figure_5_margins_and_C12.pdf', bbox_inches='tight')
plt.close()
print("Figure 5 (margins + C12) saved.")


