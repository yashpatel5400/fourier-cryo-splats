#!/usr/bin/env python3
"""Vector schematic of the conditional procedure; contains no empirical data."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['STIXGeneral'],
                     'mathtext.fontset': 'stix', 'pdf.fonttype': 42})
fig, ax = plt.subplots(figsize=(7.2, 1.60))
ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
titles = ['1  Specify inputs', '2  Choose weights', '3  Audit bounds', '4  Apply to particles']
text = ['Pilot, target, geometry\nDensity and pose classes\nNoise model',
        'Continuous weights $w$\nConservative objective\nPrimal–dual gap',
        'Continuous bias $b$\nIntegration / spectral bounds\nNoise SD (or upper bound) $s$',
        'Inference particles $y$\n$\\widehat{t} \\pm q(s,b,\\alpha)$\nNo-data fallback when wider']
colors = ['#f3f3f3', '#eef4fa', '#eef4fa', '#edf6f2']
for i, (title, body, color) in enumerate(zip(titles, text, colors)):
    x = .012+i*.251; width = .222
    ax.add_patch(FancyBboxPatch((x, .25), width, .68, boxstyle='round,pad=.008,rounding_size=.014',
        linewidth=.7, edgecolor='#52616b', facecolor=color))
    ax.text(x+width/2, .79, title, ha='center', va='center', fontsize=10, fontweight='bold')
    ax.text(x+width/2, .505, body, ha='center', va='center', fontsize=9.2, linespacing=1.5)
    if i < 3:
        ax.add_patch(FancyArrowPatch((x+width+.009, .585), (x+.251-.009, .585),
            arrowstyle='-|>', mutation_scale=9, linewidth=.7, color='#52616b'))
ax.text(.5, .095, 'Required: independent design, valid noise assumptions and membership in the stated classes.',
        ha='center', va='center', fontsize=9.5)
fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
out = ROOT/'paper/figures/conditional-workflow'
fig.savefig(out.with_suffix('.pdf'))
fig.savefig(out.with_suffix('.png'), dpi=180)
plt.close(fig)
