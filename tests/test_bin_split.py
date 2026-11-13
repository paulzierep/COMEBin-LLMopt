import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'COMEBin'))

from bin_split import split_contaminated_bins  # noqa: E402


class FakeMarkers:
    """Each marker-bearing contig carries one marker gene; completeness = distinct / 4,
    contamination = extra copies / 4 (percent)."""

    def __init__(self, gene_of):
        self.gene_of = gene_of
        self.bac_markers_on_contigs = dict.fromkeys(gene_of)
        self.ar_markers_on_contigs = {}

    def bin_quality(self, names):
        genes = [self.gene_of[n] for n in names]
        distinct = len(set(genes))
        return None, 100.0 * distinct / 4, 100.0 * (len(genes) - distinct) / 4


def _two_genomes(rng):
    # genome A: contigs 0-7, genome B: contigs 8-15, far apart in embedding space
    emb = np.vstack([rng.normal([5, 0, 0], 0.1, (8, 3)), rng.normal([0, 5, 0], 0.1, (8, 3))]).astype(np.float32)
    ids = [f'c{i}' for i in range(16)]
    lengths = np.full(16, 50_000, dtype=np.int64)  # 400 kb per genome
    gene_of = {f'c{i}': f'g{i % 4}' for i in list(range(4)) + list(range(8, 12))}
    return emb, ids, lengths, FakeMarkers(gene_of)


def test_merged_genomes_are_split():
    emb, ids, lengths, markers = _two_genomes(np.random.default_rng(0))
    labels = np.zeros(16, dtype=np.uint32)  # both genomes in one bin: comp 100, cont 100
    new, n_split = split_contaminated_bins(labels, ids, lengths, emb, markers, {0}, 200_000)
    assert n_split == 1
    assert new.dtype == np.uint32
    assert len(set(new[:8])) == 1 and len(set(new[8:])) == 1 and new[0] != new[8]


def test_clean_bins_are_kept():
    emb, ids, lengths, markers = _two_genomes(np.random.default_rng(1))
    labels = np.array([0] * 8 + [1] * 8, dtype=np.uint32)
    new, n_split = split_contaminated_bins(labels, ids, lengths, emb, markers, {0, 1}, 200_000)
    assert n_split == 0
    assert np.array_equal(new, labels)


def test_unreported_bins_are_not_touched():
    emb, ids, lengths, markers = _two_genomes(np.random.default_rng(2))
    labels = np.zeros(16, dtype=np.uint32)
    new, n_split = split_contaminated_bins(labels, ids, lengths, emb, markers, set(), 200_000)
    assert n_split == 0 and np.array_equal(new, labels)
