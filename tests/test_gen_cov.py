"""Vectorised per-view coverage statistics must match the original streaming algorithm."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'COMEBin'))
from data_aug.gen_cov import _summarize_contig, _validate_ignored_group  # noqa: E402


def reference_summary(rows, starts, ends):
    """COMEBin v1.1.0 implementation: weighted Welford update per row and view."""
    counts = np.zeros(len(starts))
    means = np.zeros(len(starts))
    m2 = np.zeros(len(starts))
    for depth_start, depth_end, depth in rows:
        for view in range(len(starts)):
            weight = max(0, min(depth_end, ends[view]) - max(depth_start, starts[view]))
            if weight <= 0:
                continue
            new_count = counts[view] + weight
            delta = depth - means[view]
            means[view] += delta * weight / new_count
            m2[view] += delta * delta * counts[view] * weight / new_count
            counts[view] = new_count
    return means, m2 / counts


def bedgraph(name, length, rng):
    cuts = np.unique(np.concatenate(([0, length], rng.integers(1, length, size=40))))
    depths = rng.integers(0, 50, size=len(cuts) - 1).astype(float) / 4
    rows = list(zip(cuts[:-1].tolist(), cuts[1:].tolist(), depths.tolist()))
    lines = [f'{name}\t{s}\t{e}\t{d:g}\n' for s, e, d in rows]
    return rows, lines


def test_matches_reference():
    rng = np.random.default_rng(1)
    for _ in range(50):
        length = int(rng.integers(1001, 5000))
        rows, lines = bedgraph('c1', length, rng)
        starts = np.concatenate(([0], rng.integers(0, length // 2, size=5)))
        ends = np.concatenate(([length], starts[1:] + rng.integers(500, length // 2, size=5)))
        ends = np.minimum(ends, length)
        means, variances = _summarize_contig('c1', iter(lines), starts, ends)
        ref_means, ref_variances = reference_summary(rows, starts, ends)
        np.testing.assert_allclose(means, ref_means, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(variances, ref_variances, rtol=1e-9, atol=1e-9)


@pytest.mark.parametrize('lines,message', [
    (['c1\t0\t10\t1\n', 'c1\t12\t20\t1\n'], 'Non-contiguous'),
    (['c1\t0\t10\t1\n'], 'Incomplete bedGraph span'),
    (['c1\t0\t20\t-1\n'], 'Invalid bedGraph interval'),
    (['c1\t0\tx\t1\n'], 'Malformed bedGraph values'),
    (['c2\t0\t20\t1\n'], 'Malformed bedGraph row'),
])
def test_rejects_bad_bedgraph(lines, message):
    with pytest.raises(ValueError, match=message):
        _summarize_contig('c1', iter(lines), np.array([0]), np.array([20]))


def test_ignored_group_must_be_contiguous():
    _validate_ignored_group('x', iter(['x\t0\t5\t1\n', 'x\t5\t9\t0\n']))
    with pytest.raises(ValueError, match='Non-contiguous'):
        _validate_ignored_group('x', iter(['x\t0\t5\t1\n', 'x\t6\t9\t0\n']))
