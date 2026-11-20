"""Precomputed bedGraph coverage (bedtools genomecov -bga) is accepted instead of BAM files (upstream #28)."""
import gzip
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'COMEBin'))
from data_aug.gen_cov import (  # noqa: E402
    CoverageContext, _stream_bam_coverage, _summarize_groups, find_coverage_files, is_bedgraph)

LINES = ['c1\t0\t10\t1\n', 'c1\t10\t20\t3\n', 'skip\t0\t5\t0\n', 'c2\t0\t30\t2\n']


def context():
    return CoverageContext(
        contig_names=('c1', 'c2'), contig_index={'c1': 0, 'c2': 1},
        contig_lengths=np.array([20, 30]),
        starts=[np.array([0, 5]), np.array([0, 10])],
        ends=[np.array([20, 15]), np.array([30, 20])])


def summarize(source):
    mean, var, seen = np.zeros((2, 2)), np.zeros((2, 2)), np.zeros(2, dtype=bool)
    if isinstance(source, Path):
        _stream_bam_coverage(source, context(), mean, var, seen)
    else:
        _summarize_groups(iter(source), 'lines', context(), mean, var, seen)
    assert seen.all()
    return mean, var


@pytest.mark.parametrize('name,compress', [('s.bedgraph', False), ('s.bg', False), ('s.bg.gz', True)])
def test_bedgraph_file_matches_stream(tmp_path, name, compress):
    path = tmp_path / name
    text = 'track type=bedGraph\n' + ''.join(LINES)
    if compress:
        with gzip.open(path, 'wt') as handle:
            handle.write(text)
    else:
        path.write_text(text)
    assert is_bedgraph(path)
    for got, expected in zip(summarize(path), summarize(LINES)):
        np.testing.assert_array_equal(got, expected)
    mean, _ = summarize(path)
    np.testing.assert_allclose(mean[:, 0], [2.0, 2.0])  # c1, whole contig and 5..15


def test_unsorted_bedgraph_is_rejected(tmp_path):
    path = tmp_path / 's.bg'
    path.write_text(''.join(['c1\t0\t20\t1\n', 'c2\t0\t30\t2\n', 'c1\t0\t20\t3\n']))
    with pytest.raises(ValueError, match='repeated groups'):
        summarize(path)


def test_find_coverage_files(tmp_path):
    for name in ['b.bam', 'a.bedgraph', 'c.bg.gz', 'notes.txt', 'b.bam.bai']:
        (tmp_path / name).write_text('')
    assert [p.name for p in find_coverage_files(tmp_path)] == ['a.bedgraph', 'b.bam', 'c.bg.gz']
