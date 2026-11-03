"""Compare two COMEBin benchmark runs (old vs new) and draw README figures.

usage: python analyze.py <bench_dir> <figure_dir>
Expects <bench_dir>/{old,new}/{run.log,mem.tsv,out/,checkm.tsv}.
"""
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import adjusted_rand_score

BENCH, FIG = Path(sys.argv[1]), Path(sys.argv[2])
FIG.mkdir(parents=True, exist_ok=True)
VERSIONS = {'old': 'v1.1.0 (comebin-optimizations-v11)', 'new': 'perf-robustness'}
COLORS = {'old': '#eb6834', 'new': '#2a78d6'}
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
STAGES = ['data augmentation', 'training', 'clustering', 'final-result selection']
TS = re.compile(r'^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d),\d+ - (.*)$')

plt.rcParams.update({
    'font.size': 10, 'axes.edgecolor': GRID, 'axes.labelcolor': INK2, 'xtick.color': INK2,
    'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False,
    'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb', 'savefig.facecolor': '#fcfcfb'})


def first_time(path, pattern):
    for line in path.read_text(errors='replace').splitlines():
        match = TS.match(line)
        if match and re.search(pattern, match.group(2)):
            return datetime.strptime(match.group(1), '%Y-%m-%d %H:%M:%S').timestamp()
    raise ValueError(f'{pattern!r} not found in {path}')


def stage_times(run):
    out = BENCH / run / 'out'
    end = int(re.search(r'end=(\d+)', (BENCH / run / 'run.log').read_text()).group(1))
    start = int(re.search(r'start=(\d+)', (BENCH / run / 'run.log').read_text()).group(1))
    log = out / 'comebin_res' / 'comebin.log'
    marks = [start, first_time(log, r'^train$'), first_time(log, r'^bin$'),
             first_time(log, r'^get_result$'), end]
    return {stage: marks[i + 1] - marks[i] for i, stage in enumerate(STAGES)}


def memory(run):
    data = np.loadtxt(BENCH / run / 'mem.tsv', ndmin=2)
    return data[:, 0] - data[0, 0], data[:, 1] / 1024 ** 2   # seconds, GiB


def checkm_quality(run):
    rows = [line.split('\t') for line in (BENCH / run / 'checkm.tsv').read_text().splitlines()[1:] if line]
    comp = np.array([float(r[-3]) for r in rows]) if rows else np.zeros(0)
    cont = np.array([float(r[-2]) for r in rows]) if rows else np.zeros(0)
    return comp, cont


def assignments(run):
    lines = (BENCH / run / 'out' / 'comebin_res' / 'comebin_res.tsv').read_text().split()
    return dict(zip(lines[0::2], lines[1::2]))


results = {}
for run in VERSIONS:
    comp, cont = checkm_quality(run)
    _, rss = memory(run)
    results[run] = {
        'stages_s': stage_times(run), 'peak_rss_gib': float(rss.max()),
        'bins': int(len(comp)),
        'hq': int(((comp > 90) & (cont < 5)).sum()),
        'mq': int(((comp > 50) & (cont < 10)).sum() - ((comp > 90) & (cont < 5)).sum()),
        'mean_completeness': float(comp.mean()) if len(comp) else 0.0,
        'mean_contamination': float(cont.mean()) if len(cont) else 0.0,
    }
    results[run]['total_s'] = sum(results[run]['stages_s'].values())
old_a, new_a = assignments('old'), assignments('new')
shared = sorted(set(old_a) | set(new_a))
results['ari_old_vs_new'] = float(adjusted_rand_score(
    [old_a.get(c, 'unbinned') for c in shared], [new_a.get(c, 'unbinned') for c in shared]))
(FIG / 'results.json').write_text(json.dumps(results, indent=2))

# 1) Wall time per stage: grouped horizontal bars, one row per stage.
fig, ax = plt.subplots(figsize=(7.5, 3.4))
y = np.arange(len(STAGES))
h = 0.36
for offset, run in [(-h / 2, 'old'), (h / 2, 'new')]:
    values = [results[run]['stages_s'][s] / 60 for s in STAGES]
    bars = ax.barh(y + offset, values, height=h - 0.04, color=COLORS[run], label=VERSIONS[run])
    for bar, value in zip(bars, values):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2, f'{value:.1f} min',
                va='center', fontsize=8.5, color=INK2)
ax.set_yticks(y, STAGES)
ax.invert_yaxis()
ax.set_xlabel('wall time (minutes)')
ax.xaxis.grid(True, color=GRID, linewidth=0.8)
ax.set_axisbelow(True)
speedup = results['old']['total_s'] / results['new']['total_s']
ax.set_title(f"Wall time per stage — total {results['old']['total_s'] / 60:.0f} → "
             f"{results['new']['total_s'] / 60:.0f} min ({speedup:.1f}× faster)",
             loc='left', fontsize=11, color=INK)
ax.legend(frameon=False, loc='lower right', fontsize=9)
fig.tight_layout()
fig.savefig(FIG / 'speed_stages.png', dpi=150)

# 2) Memory over the run (one measure, one axis).
fig, ax = plt.subplots(figsize=(7.5, 3.0))
for run in VERSIONS:
    t, rss = memory(run)
    ax.plot(t / 60, rss, color=COLORS[run], linewidth=2, label=f"{VERSIONS[run]} (peak {rss.max():.1f} GiB)")
ax.set_xlabel('elapsed (minutes)')
ax.set_ylabel('RSS of COMEBin + tools (GiB)')
ax.yaxis.grid(True, color=GRID, linewidth=0.8)
ax.set_axisbelow(True)
ax.set_ylim(bottom=0)
ax.set_title('Memory use over the run', loc='left', fontsize=11, color=INK)
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(FIG / 'memory.png', dpi=150)

# 3) Bin quality: completeness vs contamination per bin, MIMAG thresholds shaded.
fig, ax = plt.subplots(figsize=(7.5, 3.6))
for run, marker in [('old', 'o'), ('new', 's')]:
    comp, cont = checkm_quality(run)
    r = results[run]
    ax.scatter(comp, cont, s=36, marker=marker, facecolor=COLORS[run], edgecolor='#fcfcfb', linewidth=1.2,
               alpha=0.85, label=f"{VERSIONS[run]}: {r['bins']} bins, {r['hq']} HQ, {r['mq']} MQ")
ax.axhline(5, color=INK2, linewidth=0.8, linestyle=':')
ax.axhline(10, color=INK2, linewidth=0.8, linestyle=':')
ax.axvline(50, color=INK2, linewidth=0.8, linestyle=':')
ax.axvline(90, color=INK2, linewidth=0.8, linestyle=':')
ax.text(91, ax.get_ylim()[1] * 0.95 if ax.get_ylim()[1] > 12 else 11.5, 'HQ ≥90 / <5', fontsize=8, color=INK2, va='top')
ax.text(51, ax.get_ylim()[1] * 0.95 if ax.get_ylim()[1] > 12 else 11.5, 'MQ ≥50 / <10', fontsize=8, color=INK2, va='top')
ax.set_xlabel('CheckM completeness (%)')
ax.set_ylabel('CheckM contamination (%)')
ax.set_xlim(0, 101)
ax.set_ylim(bottom=0)
ax.set_title(f"Bin quality (CheckM, Bacteria marker set) — agreement ARI {results['ari_old_vs_new']:.2f}",
             loc='left', fontsize=11, color=INK)
ax.legend(frameon=False, fontsize=9, loc='upper left')
fig.tight_layout()
fig.savefig(FIG / 'quality.png', dpi=150)
print(json.dumps(results, indent=2))
