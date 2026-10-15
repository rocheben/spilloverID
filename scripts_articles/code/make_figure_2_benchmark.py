"""Generate Figure 2 (mechanistic benchmark across 11 simulated scenarios)."""
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

THRESHOLD = 0.55

with open('../data/benchmark_results.json', 'r') as f:
    bench = json.load(f)

summary = bench['summary']
per_rep = bench['per_rep']

scen_order = ['S1_canonical', 'S2_bridge', 'S3_cryptic', 'S4_confounded',
              'S5_codistributed', 'S6_strong_seasonal', 'S7_low_power',
              'S8_measurement_error', 'S9_weak_signal', 'S10_short_immunity',
              'S11_no_spillover']

short_labels = {
    'S1_canonical': 'S1\ncanonical',
    'S2_bridge': 'S2\nbridge\nhost',
    'S3_cryptic': 'S3\ncryptic',
    'S4_confounded': 'S4\nconfounded',
    'S5_codistributed': 'S5\nco-distrib.',
    'S6_strong_seasonal': 'S6\nseasonal',
    'S7_low_power': 'S7\nlow\npower',
    'S8_measurement_error': 'S8\nmeas.\nerror',
    'S9_weak_signal': 'S9\nweak\nsignal',
    'S10_short_immunity': 'S10\nshort\nimmunity',
    'S11_no_spillover': 'S11\nNEG\ncontrol'
}

summ_by_key = {s['scenario']: s for s in summary}
pos_scens = [s for s in scen_order if s != 'S11_no_spillover']
top1 = [summ_by_key[k]['top1_rate'] for k in pos_scens]
top2 = [summ_by_key[k]['top2_rate'] for k in pos_scens]
margins_mean = [summ_by_key[k]['mean_margin'] for k in pos_scens]

margins_by_scen = {}
max_pi_by_scen = {}
for row in per_rep:
    s = row['scenario']
    if row.get('margin') is not None:
        margins_by_scen.setdefault(s, []).append(row['margin'])
    if row.get('max_pi') is not None:
        max_pi_by_scen.setdefault(s, []).append(row['max_pi'])

# Negative-control false-positive rate (computed from data, not invented)
s11_maxpi = max_pi_by_scen['S11_no_spillover']
fp_rate = 100.0 * sum(1 for m in s11_maxpi if m >= THRESHOLD) / len(s11_maxpi)

COL_TOP1 = '#009E73'
COL_TOP2 = '#56B4E9'
COL_NEG = '#888888'

fig = plt.figure(figsize=(9.0, 10.2))
gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 0.95, 0.95], hspace=1.15,
                      top=0.885, bottom=0.06, left=0.10, right=0.97)

x_pos = np.arange(len(pos_scens))          # positions for S1..S10
x_neg = len(pos_scens)                     # position for S11
all_ticks = list(range(len(scen_order)))   # 0..10

# ----- Panel a: identification rate (+ S11 false-positive rate) -----
ax = fig.add_subplot(gs[0])
width = 0.38
bars1 = ax.bar(x_pos - width/2, top1, width, color=COL_TOP1, edgecolor='black',
               linewidth=0.7, label='Top-1 (true reservoir at rank 1)')
bars2 = ax.bar(x_pos + width/2, top2, width, color=COL_TOP2, edgecolor='black',
               linewidth=0.7, label='Top-2 (true reservoir in top 2)')
bar_neg = ax.bar(x_neg, fp_rate, width, color=COL_NEG, edgecolor='black',
                 linewidth=0.7, hatch='//',
                 label='S11: false-positive rate (max \u03A0 \u2265 0.55)')

ax.axhline(25, color='#D55E00', linestyle=':', linewidth=1.0, alpha=0.7)
ax.axhline(50, color='#D55E00', linestyle=':', linewidth=1.0, alpha=0.5)
ax.text(len(scen_order) - 1.55, 27.0, 'top-1 chance (1/4)', fontsize=6.5, color='#D55E00', ha='right')
ax.text(len(scen_order) - 1.55, 52.0, 'top-2 chance (2/4)', fontsize=6.5, color='#D55E00', ha='right')

for b, v in zip(bars1, top1):
    ax.text(b.get_x() + b.get_width()/2, v + 1.8, f'{v:.0f}', ha='center',
            fontsize=6.8, fontweight='bold')
for b, v in zip(bars2, top2):
    ax.text(b.get_x() + b.get_width()/2, v + 1.8, f'{v:.0f}', ha='center',
            fontsize=6.8, color='dimgray')
ax.text(x_neg, fp_rate + 1.8, f'{fp_rate:.0f}', ha='center', fontsize=6.8,
        fontweight='bold', color='#444444')

ax.set_xticks(all_ticks)
ax.set_xticklabels([short_labels[k] for k in scen_order], fontsize=7)
ax.set_ylabel('Identification rate (%)\n(S11: false-positive rate)')
ax.set_ylim(0, 118)
ax.set_xlim(-0.7, len(scen_order) - 0.3)
ax.set_title('a   Reservoir identification rate per scenario (300 replicates each)',
             loc='left', fontsize=9, fontweight='bold')
ax.legend(loc='lower center', ncol=3, frameon=False, fontsize=7.2,
          bbox_to_anchor=(0.5, -0.62))

# ----- Panel b: discrimination margin (S11 flagged n/a) -----
ax = fig.add_subplot(gs[1])
data_for_violins = [margins_by_scen[k] for k in pos_scens]
parts = ax.violinplot(data_for_violins, positions=list(x_pos), showmeans=True,
                      showmedians=False, widths=0.7)
for i, pc in enumerate(parts['bodies']):
    mean_m = margins_mean[i]
    if mean_m > 0.15: col = '#009E73'
    elif mean_m > 0.05: col = '#E69F00'
    else: col = '#D55E00'
    pc.set_facecolor(col); pc.set_edgecolor('black'); pc.set_alpha(0.8)
parts['cmeans'].set_color('black'); parts['cmeans'].set_linewidth(1.5)
parts['cmins'].set_color('gray'); parts['cmaxes'].set_color('gray')
parts['cbars'].set_color('gray')

ax.axhline(0, color='#D55E00', linestyle='--', linewidth=1.0, alpha=0.6)
# S11: margin undefined (no true reservoir) -> explicit marker, no fabricated value
ax.text(x_neg, 0.0, 'n/a\n(no true\nreservoir)', ha='center', va='center',
        fontsize=6.8, color='#555555', style='italic',
        bbox=dict(boxstyle='round,pad=0.25', facecolor='#eeeeee',
                  edgecolor='#aaaaaa', linewidth=0.6))

ax.set_xticks(all_ticks)
ax.set_xticklabels([short_labels[k] for k in scen_order], fontsize=7)
ax.set_ylabel(r'Margin $\Pi_{\mathrm{truth}} - \max \Pi_{\mathrm{others}}$')
ax.set_title('b   Discrimination margin distribution per scenario',
             loc='left', fontsize=9, fontweight='bold')
ax.set_ylim(-0.35, 0.65)
ax.set_xlim(-0.7, len(scen_order) - 0.3)
mean_patches = [
    mpatches.Patch(facecolor='#009E73', edgecolor='black', label='mean margin > 0.15'),
    mpatches.Patch(facecolor='#E69F00', edgecolor='black', label='0.05\u20130.15'),
    mpatches.Patch(facecolor='#D55E00', edgecolor='black', label='< 0.05'),
]
ax.legend(handles=mean_patches, loc='lower center', ncol=3, frameon=False,
          fontsize=7.2, bbox_to_anchor=(0.5, -0.62))

# ----- Panel c: max Pi per replicate (all 11 incl. S11) -----
ax = fig.add_subplot(gs[2])
data_max = [max_pi_by_scen[k] for k in scen_order]
parts = ax.violinplot(data_max, positions=all_ticks,
                      showmeans=True, showmedians=False, widths=0.7)
for i, pc in enumerate(parts['bodies']):
    if scen_order[i] == 'S11_no_spillover':
        pc.set_facecolor(COL_NEG)
    else:
        pc.set_facecolor('#009E73')
    pc.set_edgecolor('black'); pc.set_alpha(0.85)
parts['cmeans'].set_color('black'); parts['cmeans'].set_linewidth(1.5)
parts['cmins'].set_color('gray'); parts['cmaxes'].set_color('gray')
parts['cbars'].set_color('gray')

s11 = max_pi_by_scen['S11_no_spillover']
s11_lo, s11_hi = np.quantile(s11, [0.025, 0.975])
ax.axhspan(s11_lo, s11_hi, color='gray', alpha=0.16, zorder=0)

ax.axhline(THRESHOLD, color='#D55E00', linestyle='--', linewidth=1.0, alpha=0.8)
ax.text(len(scen_order) - 0.4, THRESHOLD + 0.01,
        f'suggested threshold \u03A0 \u2248 {THRESHOLD:.2f}',
        fontsize=6.8, color='#D55E00', va='bottom', ha='right')

ax.set_xticks(all_ticks)
ax.set_xticklabels([short_labels[k] for k in scen_order], fontsize=7)
ax.set_ylabel(r'$\max_k \Pi_k$ per replicate')
ax.set_title('c   Maximum \u03A0 per replicate: positive scenarios vs. negative control (S11)',
             loc='left', fontsize=9, fontweight='bold')
ax.set_ylim(0.15, 0.95)
ax.set_xlim(-0.7, len(scen_order) - 0.3)

neg_patch = mpatches.Patch(facecolor=COL_NEG, edgecolor='black',
                           label='Negative control S11 (no spillover)')
pos_patch = mpatches.Patch(facecolor='#009E73', edgecolor='black',
                           label='Positive scenarios (true reservoir exists)')
band_patch = mpatches.Patch(facecolor='gray', alpha=0.16,
                            label=f'S11 95% range [{s11_lo:.2f}, {s11_hi:.2f}]')
ax.legend(handles=[pos_patch, neg_patch, band_patch], loc='lower left',
          fontsize=7.0, frameon=False)

# Suptitle removed: keep figure focused on panels only

plt.savefig('../figures/figure_2_benchmark.png', dpi=300, bbox_inches='tight')
plt.savefig('../figures/figure_2_benchmark.pdf', bbox_inches='tight')
plt.close()
print("Figure 2 saved.")
