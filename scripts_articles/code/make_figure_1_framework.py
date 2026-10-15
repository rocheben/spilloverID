"""Generate Figure 1 (14-criterion framework diagram) for the spilloverID Nature Microbiology submission."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
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

# ===================== FIGURE 1: Framework diagram =====================
fig, ax = plt.subplots(figsize=(9.0, 5.0))

categories = [
    ('I. Spatial co-distribution',         ['C1', 'C2'],         '#a6cee3'),
    ('II. Individual exposure',            ['C3', 'C4', 'C5'],   '#1f78b4'),
    ('III. Multi-level coupling',          ['C6'],               '#b2df8a'),
    ('IV. Demographic force of infection', ['C7'],               '#33a02c'),
    ('V. Temporal precedence',             ['C8'],               '#fb9a99'),
    ('VI. Mechanistic causal models',      ['C9', 'C10'],        '#e31a1c'),
    ('VII. Specificity controls',          ['C11'],              '#fdbf6f'),
    ('VIII. Maintenance signatures',       ['C12','C13','C14'],  '#ff7f00'),
]

x0 = 0.4
y_top = 8.0
h = 0.78
gap = 0.05
w_cat = 6.5

for i, (name, crits, col) in enumerate(categories):
    y = y_top - i * (h + gap)
    rect = FancyBboxPatch((x0, y), w_cat, h, boxstyle='round,pad=0.02',
                          facecolor=col, edgecolor='black', linewidth=0.8, alpha=0.92)
    ax.add_patch(rect)
    ax.text(x0 + 0.15, y + h/2, name, fontsize=8.5, va='center', fontweight='bold')
    ax.text(x0 + w_cat - 0.15, y + h/2, ', '.join(crits), fontsize=8, va='center', ha='right',
            family='monospace')

arrow_y = y_top - 3.5 * (h + gap) + h/2
ax.annotate('', xy=(9.4, arrow_y), xytext=(x0 + w_cat + 0.25, arrow_y),
            arrowprops=dict(arrowstyle='->', lw=2, color='black'))

sc_x, sc_y, sc_w, sc_h = 9.55, arrow_y - 2.10, 3.0, 4.2
sc_box = FancyBboxPatch((sc_x, sc_y), sc_w, sc_h, boxstyle='round,pad=0.05',
                        facecolor='#fff5b8', edgecolor='black', linewidth=1.3)
ax.add_patch(sc_box)
ax.text(sc_x + sc_w/2, sc_y + sc_h - 0.50,
        r'$\Pi_k = \sum_{c=1}^{14} w_c \cdot s_{k,c}$',
        fontsize=12, ha='center', va='center', fontweight='bold')
ax.text(sc_x + sc_w/2, sc_y + sc_h - 1.85, 'Composite\nplausibility\nscore',
        fontsize=9, ha='center', va='center')
ax.text(sc_x + sc_w/2, sc_y + 0.50, 'for each\ncandidate host\n$k = 1,\\dots,K$',
        fontsize=8, ha='center', va='center', style='italic')

ax.text(x0 + w_cat/2, y_top + h + 0.4, '14 epidemiological criteria across 8 categories',
        fontsize=10, ha='center', fontweight='bold')

ax.set_xlim(0, 13.0)
ax.set_ylim(min(0.6, sc_y - 0.25), 9.5)
ax.axis('off')
plt.savefig('../figures/figure_1_framework.png', dpi=300, bbox_inches='tight')
plt.savefig('../figures/figure_1_framework.pdf', bbox_inches='tight')
plt.close()
print("Figure 1 (framework) saved.")


