"""Figure of component micro-benchmarks (measured 2026-11-03, 6 CPU threads, laptop).

Numbers come from benchmark/README.md; rerun the commands listed there to refresh them.
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OLD, NEW = '#eb6834', '#2a78d6'
INK, INK2, GRID, SURFACE = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
PANELS = [  # title, unit, old, new
    ('InfoNCE loss step\n(batch 1024 × 6 views)', 's', 2.26, 0.50),
    ('Loss step peak memory', 'GB', 2.72, 0.92),
    ('Coverage features\n(4,000 contigs, 1 BAM)', 's', 12.9, 7.4),
    ('Batch loading per epoch\n(68,260 contigs)', 's', 0.394, 0.030),
]
plt.rcParams.update({'font.size': 9.5, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.edgecolor': GRID, 'xtick.color': INK2, 'ytick.color': INK2,
                     'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE})
fig, axes = plt.subplots(1, len(PANELS), figsize=(10, 2.6))
for ax, (title, unit, old, new) in zip(axes, PANELS):
    bars = ax.barh([0, 1], [old, new], color=[OLD, NEW], height=0.62)
    for bar, value in zip(bars, [old, new]):
        ax.text(bar.get_width() * 1.03, bar.get_y() + bar.get_height() / 2,
                f'{value:g} {unit}', va='center', fontsize=8.5, color=INK2)
    ax.set_yticks([0, 1], ['v1.1.0', 'this branch'])
    ax.invert_yaxis()
    ax.set_xlim(0, old * 1.45)
    ax.set_xticks([])
    ax.spines['bottom'].set_visible(False)
    ax.set_title(f'{title}\n{old / new:.1f}× less', fontsize=9.5, color=INK, loc='left')
for ax in axes[1:]:
    ax.set_yticklabels([])
fig.tight_layout()
out = Path(__file__).resolve().parent / 'figures' / 'components.png'
out.parent.mkdir(exist_ok=True)
fig.savefig(out, dpi=150)
print(out)
