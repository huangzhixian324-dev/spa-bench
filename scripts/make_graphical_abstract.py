"""Graphical Abstract for SPATBench (Cell Systems submission).
Spec: 5:2 aspect ratio, 300 DPI, 3750x1500 px (= 317.5 x 127 mm). Layout:
CDS workflow (left) + three findings (right) + headline claim (bottom).
Style matches make_figures_v33.py.
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

FIGD = Path(__file__).resolve().parents[1] / 'results' / 'figures' / 'v33'
FIGD.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    'font.size': 13, 'font.family': 'DejaVu Sans',
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': None})

FIG_W, FIG_H = 12.5, 5.0  # inches -> 5:2 at 300dpi = 3750x1500
C_MAIN = '#4C72B0'
C_WARN = '#C44E52'
C_OK = '#55A868'
C_GREY = '#8A8A8A'

fig = plt.figure(figsize=(FIG_W, FIG_H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 40)
ax.axis('off')

# ---------- Title banner ----------
ax.add_patch(FancyBboxPatch((1, 35.2), 98, 4.2, boxstyle='round,pad=0.3',
                            fc='#1F3B60', ec='none'))
ax.text(50, 37.3, 'SPATBench: Methodological Choices Can Outweigh '
        'Algorithmic Choices in ICI Transcriptomic Prediction',
        ha='center', va='center', color='white', fontsize=15, weight='bold')

# ---------- Left panel: CDS workflow ----------
ax.add_patch(FancyBboxPatch((1.5, 4.5), 26, 28.5, boxstyle='round,pad=0.4',
                            fc='#F4F6FA', ec='#1F3B60', lw=1.2))
ax.text(14.5, 31.2, 'Circularity Detection Score (CDS)', ha='center',
        fontsize=12.5, weight='bold', color='#1F3B60')

steps = [
    ('Input: expression matrix\n+ response labels\n+ DECLARE endpoint genes', '#FFFFFF'),
    ('C1  MI(genes, label)   30%\nC2  |Spearman (genes, label)|   30%\nC3  top-500 / random MI   40%', '#FFFFFF'),
    ('CDS = 0.3·C1n + 0.3·C2 + 0.4·C3', '#DCE6F2'),
]
y = 27.5
for txt, fc in steps:
    ax.add_patch(FancyBboxPatch((3, y - 4.2), 23, 5.2, boxstyle='round,pad=0.25',
                                fc=fc, ec='#9AA7B8', lw=0.9))
    ax.text(14.5, y - 1.6, txt, ha='center', va='center', fontsize=9.5)
    if y > 13:
        ax.add_patch(FancyArrowPatch((14.5, y - 4.4), (14.5, y - 5.9),
                                     arrowstyle='-|>', mutation_scale=14,
                                     color='#1F3B60', lw=1.4))
    y -= 6.2
ax.add_patch(FancyBboxPatch((3, y - 4.2), 23, 5.0, boxstyle='round,pad=0.25',
                            fc='#1F3B60', ec='none'))
ax.text(14.5, y - 1.7, 'null-referenced:\nvs cohort-matched label shuffles',
        ha='center', va='center', fontsize=9.5, color='white', weight='bold')

# ---------- Right: three findings ----------
cards = [
    dict(x0=30.5, title='1  Endpoint switch collapses',
         lines=[('Cytolytic surrogate vs RECIST v1.1:', 'k', 8.6),
                ('EN(Var) 0.985 → 0.465', '#C44E52', 10),
                ('GEP  0.968 → 0.503', '#C44E52', 10),
                ('IMPRES (no overlap): unchanged', '#55A868', 8.6),
                ('Effect = 3.5x method gap', 'k', 8.6)]),
    dict(x0=53.0, title='2  Feature selection rewrites',
         lines=[('Var- vs MI-based selection:', 'k', 8.6),
                ('27% gene overlap', 'k', 8.6),
                ('yet indistinguishable AUROC', '#1F3B60', 9.6),
                ('→ the narrative changes,', 'k', 8.6),
                ('not the ranking', 'k', 8.6)]),
    dict(x0=75.5, title='3  CDS on real cohorts',
         lines=[('Circular cytolytic:', '#55A868', 8.6),
                ('CDS 0.954 >> null 0.842', '#55A868', 9.6),
                ('Clinical: inside nulls', '#C44E52', 8.6),
                ('Labels non-specific', 'k', 8.6),
                ('(60/60 pseudo-endpoints HIGH)', 'k', 8.0)]),
]
for c in cards:
    ax.add_patch(FancyBboxPatch((c['x0'], 16.5), 21, 16.5,
                                boxstyle='round,pad=0.4', fc='#F7F9FB',
                                ec='#1F3B60', lw=1.1))
    ax.text(c['x0'] + 10.5, 31.2, c['title'], ha='center', fontsize=9.8,
            weight='bold', color='#1F3B60')
    yy = 28.8
    for txt, col, fs in c['lines']:
        ax.text(c['x0'] + 10.5, yy, txt, ha='center', va='center',
                fontsize=fs, color=col)
        yy -= 2.55

# ---------- Findings 1 & 2 arrow to 3 ----------
ax.add_patch(FancyArrowPatch((42.5, 24.5), (44.8, 24.5), arrowstyle='-|>',
                             mutation_scale=16, color=C_GREY, lw=1.6))
ax.add_patch(FancyArrowPatch((65.0, 24.5), (67.3, 24.5), arrowstyle='-|>',
                             mutation_scale=16, color=C_GREY, lw=1.6))

# ---------- Bottom: Gide highlight + headline ----------
ax.add_patch(FancyBboxPatch((30.5, 4.5), 68.5, 9.0, boxstyle='round,pad=0.4',
                            fc='#EAF3EA', ec='#55A868', lw=1.4))
ax.text(64.75, 11.7, 'High-resolution validation: ElasticNet (MI) on Gide '
        '2019 (n = 73) becomes FDR-significant',
        ha='center', fontsize=10.0, weight='bold', color='#2F5D3A')
ax.text(64.75, 9.3, 'at 5,000 full-pipeline shuffles',
        ha='center', fontsize=10.0, weight='bold', color='#2F5D3A')
ax.text(64.75, 6.6, 'p = 0.0232, within-cohort BH q = 0.0278      |      '
        'no trainable method survives FDR on any clinical endpoint at n ≤ 43',
        ha='center', fontsize=9.2, color='#2F5D3A')
ax.add_patch(FancyArrowPatch((14.5, 4.2), (14.5, 4.2), arrowstyle='-',
                             color='none'))

ax.add_patch(FancyBboxPatch((1.5, 4.5), 26, 0.01, boxstyle='round,pad=0.1',
                            fc='none', ec='none'))

fig.savefig(FIGD / 'graphical_abstract.png', dpi=300, facecolor='white')
fig.savefig(FIGD / 'graphical_abstract.pdf', facecolor='white')
print('graphical abstract written to', FIGD)
print('graphical abstract written (PNG + PDF)')
