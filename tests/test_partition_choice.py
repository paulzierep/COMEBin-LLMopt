import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'COMEBin'))
gfr = pytest.importorskip('get_final_result')


class FakeMarkers:
    """Completeness = 60 % per marker contig in the bin; no contamination."""

    def bin_quality(self, names):
        return None, 60.0 * len(names), 0.0


def test_small_bins_do_not_count_toward_partition_score():
    # contigs: a0, a1 carry markers; lengths make bin 0 small (< 200 kb) in partition "split"
    labels_small = np.array([0, 0, 1])   # bin 0 = {a0, a1}: "120 %" complete but only 1 kb
    labels_big = np.array([0, 1, 1])     # bin 1 = {a1, x}: 60 % complete, 300 kb
    lengths = np.array([500, 500, 300000])
    marker_idx, marker_names = np.array([0, 1]), ('a0', 'a1')
    m = FakeMarkers()
    s_all, _ = gfr.score_marker_labels(labels_small, marker_idx, marker_names, m)
    assert s_all['num_9010'] == 1
    keep = gfr.retained_final_labels(labels_small, lengths)
    s_kept, q = gfr.score_marker_labels(labels_small, marker_idx, marker_names, m, keep)
    assert s_kept['sum'] == 0 and q == {}
    keep_big = gfr.retained_final_labels(labels_big, lengths)
    s_big, _ = gfr.score_marker_labels(labels_big, marker_idx, marker_names, m, keep_big)
    assert s_big['num_5010'] == 1
