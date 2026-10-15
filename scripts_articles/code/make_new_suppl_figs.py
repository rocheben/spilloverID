"""New supplementary figures for the two added sensitivity analyses.
S-Fig A: composite vs single-criterion identification across regimes (simulation vs real cases).
S-Fig B: robustness of the 0.55 threshold under varying background seroprevalence x specificity.
Wong (2011) colour-blind-safe palette; sans-serif; no gridlines; RGB; PNG (300 dpi) + vector PDF.
"""
import json, csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'DejaVu Sans'  # sans-serif (matches existing figure set)
rcParams['font.size'] = 7
rcParams['axes.labelsize'] = 7.5
rcParams['axes.titlesize'] = 8.5
rcParams['xtick.labelsize'] = 6.5
rcParams['ytick.labelsize'] = 6.5
rcParams['legend.fontsize'] = 6.5
rcParams['savefig.dpi'] = 300
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False
rcParams['axes.grid'] = False

# Wong palette
W = dict(green='#009E73', sky='#56B4E9', verm='#D55E00', orange='#E69F00',
         blue='#0072B2', pink='#CC79A7', yellow='#F0E442', black='#000000', grey='#999999')

DATA='../data'
OUT='../figures'

# ---------- read sim top-1 per item ----------
sim={}
with open(f'{DATA}/fig_compindiv_sim.csv') as f:
    for row in csv.DictReader(f):
        v=row['sim_top1']
        sim[row['item']]=float(v) if v not in ('NA','') else np.nan

# ---------- real-case identification (/3), recomputed from case_results.json ----------
cd=json.load(open(f'{DATA}/case_results.json'))
crit=['C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14']
def real_id(c):
    s=0; seen=0
    for k in ('sn','nv','rv'):
        cc=cd[k]; hosts=cc['hosts']; ti=hosts.index(cc['true_host'])
        vals=[cc['criteria'][i].get(c) for i in range(len(hosts))]
        if all(v is None for v in vals): return None
        vv=[-1e9 if v is None else v for v in vals]; mx=max(vv)
        s+= int(vv[ti]==mx and sum(1 for x in vv if x==mx)==1); seen+=1
    return s
realid={c:real_id(c) for c in crit}
# composite: the number of real cases in which it ranks the true host strictly first
def composite_real():
    s=0
    for k in ('sn','nv','rv'):
        cc=cd[k]; ti=cc['hosts'].index(cc['true_host'])
        s+=int(cc['rank'][ti]==1)
    return s
comp_real=composite_real()

# =================== FIGURE A: cross-regime ===================
figA, ax = plt.subplots(figsize=(5.0, 4.2))
# criteria points (only those evaluable on real cases)
xs=[]; ys=[]; labs=[]
rng=np.random.default_rng(1)
for c in crit:
    if realid[c] is None or np.isnan(sim.get(c,np.nan)): continue
    xs.append(sim[c]); ys.append(realid[c]); labs.append(c)
xs=np.array(xs); ys=np.array(ys, float)
yj=ys + rng.uniform(-0.10,0.10,size=len(ys))
ax.scatter(xs, yj, s=42, facecolor=W['sky'], edgecolor='black', linewidth=0.5, zorder=3)
off={'C12':(7,-9),'C13':(7,4),'C1':(-14,2),'C6':(6,-9),'C9':(6,3),'C2':(6,3),'C14':(6,3),'C7':(6,3),'C5':(6,-10),'C3':(6,3),'C4':(6,3),'C10':(6,3)}
for x,y,l in zip(xs,yj,labs):
    dx,dy=off.get(l,(5,3))
    ax.annotate(l, (x,y), textcoords='offset points', xytext=(dx,dy), fontsize=6.3, color='#333333')
# composite
ax.scatter([sim['composite']],[comp_real], s=200, marker='*', facecolor=W['green'],
           edgecolor='black', linewidth=0.8, zorder=5)
ax.annotate('composite '+r'$\Pi_k$', (sim['composite'],comp_real),
            textcoords='offset points', xytext=(-6,-16), fontsize=7.5,
            fontweight='bold', color=W['green'], ha='right')
# guide region (good in both)
ax.axhspan(2.5,3.4, xmin=0.0, xmax=1.0, color=W['green'], alpha=0.06, zorder=0)
ax.text(40, 3.28, 'strong on real cases (top band)', fontsize=6.2, color=W['green'], va='center')
ax.set_xlabel('Top-1 identification on simulated benchmark (%, 10 positive scenarios)')
ax.set_ylabel('Real cases identified at rank 1 (out of 3)')
ax.set_xlim(-3,104); ax.set_ylim(-0.4,3.5)
ax.set_yticks([0,1,2,3])
ax.set_title('Composite vs single criteria: only the composite is strong in both regimes',
             loc='left', fontweight='bold', fontsize=8)
ax.text(0.0,-0.20,
        'Single criteria that excel on simulated data (e.g. C7, C5) rely on idealized signals and fail on real cases;\n'
        'criteria that identify all real reservoirs (C12, C13) are weak alone on simulated data. Only $\\Pi_k$ is high on both axes.\n'
        'Criteria that cannot be evaluated on the reconstructions (C6, C8, C11) are omitted.',
        transform=ax.transAxes, fontsize=5.6, color='#555555', va='top')
figA.subplots_adjust(left=0.11, right=0.97, top=0.92, bottom=0.20)
figA.savefig(f'{OUT}/suppl_figure_4_composite_vs_individual.png', dpi=300, bbox_inches='tight')
figA.savefig(f'{OUT}/suppl_figure_4_composite_vs_individual.pdf', bbox_inches='tight')
plt.close(figA)

# =================== FIGURE B: threshold robustness ===================
rows=[]
with open(f'{DATA}/fig_threshold_grid.csv') as f:
    for r in csv.DictReader(f):
        rows.append({k:(float(v) if k!='n' else int(float(v))) for k,v in r.items()})
alphas=sorted(set(r['alpha'] for r in rows))
specs=sorted(set(r['spec'] for r in rows), reverse=True)  # 0.99 at top
M=np.full((len(specs),len(alphas)), np.nan)   # colour = 97.5th pct max Pi
BG=np.full((len(specs),len(alphas)), np.nan)
FP=np.full((len(specs),len(alphas)), np.nan)
MED=np.full((len(specs),len(alphas)), np.nan)
for r in rows:
    i=specs.index(r['spec']); j=alphas.index(r['alpha'])
    M[i,j]=r['maxpi_975']; BG[i,j]=r['bg_sero']; FP[i,j]=r['fp_055']; MED[i,j]=r['maxpi_median']

figB, ax = plt.subplots(figsize=(5.4,3.6))
im=ax.imshow(M, aspect='auto', cmap='cividis', vmin=0.35, vmax=0.58)  # CVD-safe sequential
for i in range(len(specs)):
    for j in range(len(alphas)):
        ax.text(j,i, f'med {MED[i,j]:.2f}\n97.5% {M[i,j]:.2f}\nFP {FP[i,j]:.0f}%\nbg {BG[i,j]:.0f}%',
                ha='center', va='center', fontsize=5.6,
                color='white' if M[i,j]<0.50 else 'black')
ax.set_xticks(range(len(alphas)))
ax.set_xticklabels([f'{a:g}' for a in alphas])
ax.set_yticks(range(len(specs)))
ax.set_yticklabels([f'{s:.2f}' for s in specs])
ax.set_xlabel('Cryptic background FOI (alpha_cryptic_bg)')
ax.set_ylabel('Assay specificity')
ax.set_title('Negative-control max '+r'$\Pi$'+' stays below the 0.55 threshold across background and specificity',
             loc='left', fontweight='bold', fontsize=7.6)
cb=figB.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
cb.set_label('97.5th percentile of max '+r'$\Pi$'+' (negative control)', fontsize=6.5)
cb.ax.axhline(0.55, color=W['verm'], linewidth=1.6)
cb.set_ticks([0.40,0.45,0.50,0.55])
cb.ax.text(0.5, 0.55, ' 0.55 threshold', color=W['verm'], fontsize=6, va='bottom', ha='left', transform=cb.ax.get_yaxis_transform())
figB.text(0.01,0.01,
          'Each cell: 20 negative-control replicates. med = median max '+r'$\Pi$'+'; 97.5% = upper percentile; '
          'FP = false-positive rate at '+r'$\Pi\geq0.55$'+'; bg = realized human background seroprevalence.',
          fontsize=5.4, color='#555555')
figB.subplots_adjust(left=0.10, right=1.0, top=0.90, bottom=0.16)
figB.savefig(f'{OUT}/suppl_figure_5_threshold_robustness.png', dpi=300, bbox_inches='tight')
figB.savefig(f'{OUT}/suppl_figure_5_threshold_robustness.pdf', bbox_inches='tight')
plt.close(figB)
print("Saved both supplementary figures (PNG + PDF).")
print("composite sim_top1=%.1f%% real=3/3 ; best single sim=C7 %.0f%% real=%d/3 ; C12 real=3/3 sim=%.0f%%"
      % (sim['composite'], sim['C7'], realid['C7'], sim['C12']))
