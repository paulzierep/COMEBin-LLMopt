import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'COMEBin'))
from checkm2_select import final_labels, greedy_select  # noqa: E402


def fake_quality(truth):
    """Completeness = share of a genome's contigs; contamination = share of other genomes' contigs."""
    sizes = {g: sum(1 for t in truth.values() if t == g) for g in set(truth.values())}

    def q(bins):
        out = []
        for b in bins:
            counts = {}
            for c in b:
                counts[truth[c]] = counts.get(truth[c], 0) + 1
            g = max(counts, key=counts.get)
            out.append((100 * counts[g] / sizes[g], 100 * (len(b) - counts[g]) / sizes[g]))
        return out
    return q


def test_greedy_picks_clean_bins_across_partitions():
    truth = {f'a{i}': 'A' for i in range(4)} | {f'b{i}': 'B' for i in range(4)}
    lengths = {c: 100 for c in truth}
    # partition 1 merges A and B; partition 2 splits them cleanly but cuts B in half
    bins = [set(truth),
            {f'a{i}' for i in range(4)}, {'b0', 'b1'}, {'b2', 'b3'},
            {f'b{i}' for i in range(4)}]
    chosen, used = greedy_select(bins, lengths, fake_quality(truth), weight=4, min_bin_bp=0, min_comp=50)
    assert {frozenset(b) for b in chosen} == {frozenset(f'a{i}' for i in range(4)),
                                              frozenset(f'b{i}' for i in range(4))}
    assert used == set(truth)


def test_ties_go_to_the_larger_bin():
    truth = {'x1': 'X', 'x2': 'X', 'x3': 'X'}
    lengths = {'x1': 100, 'x2': 100, 'x3': 5}
    def q(bins):  # same score for both candidates
        return [(100.0, 0.0) for _ in bins]
    chosen, _ = greedy_select([{'x1', 'x2'}, {'x1', 'x2', 'x3'}], lengths, q, min_bin_bp=0)
    assert chosen[0] == {'x1', 'x2', 'x3'}


def test_stops_below_min_completeness():
    lengths = {'c': 100}
    chosen, used = greedy_select([{'c'}], lengths, lambda bins: [(40.0, 0.0)], min_bin_bp=0, min_comp=50)
    assert chosen == [] and used == set()


def test_final_labels_fill_from_base_partition():
    ids = ['a', 'b', 'c', 'd', 'e']
    base = np.array([7, 7, 3, 3, 9])
    labels = final_labels(ids, [{'a', 'c'}], {'a', 'c'}, base)
    assert labels[0] == labels[2] == 0          # selected bin
    assert labels[1] != labels[3]               # b and d keep separate base bins
    assert len({labels[1], labels[3], labels[4]}) == 3
    assert min(labels[[1, 3, 4]]) >= 1


def test_low_completeness_bins_are_not_rescored():
    truth = {f'a{i}': 'A' for i in range(4)} | {f'b{i}': 'B' for i in range(10)}
    lengths = {c: 100 for c in truth}
    q = fake_quality(truth)
    scored = []

    def counting(bins):
        scored.extend(frozenset(b) for b in bins)
        return q(bins)
    # {a0..a3, b0}: A complete, touched when A is picked; {b0, b1}: 20 % of B, pruned at once
    bins = [{f'a{i}' for i in range(4)}, {'a0', 'a1', 'a2', 'a3', 'b0'}, {'b0', 'b1'}]
    chosen, _ = greedy_select(bins, lengths, counting, weight=4, min_bin_bp=0, min_comp=50, prune_margin=5)
    assert chosen == [{f'a{i}' for i in range(4)}]
    assert scored.count(frozenset({'b0', 'b1'})) == 1
    assert frozenset({'b0'}) in scored  # the A+b0 bin minus A is rescored once
    chosen_all, _ = greedy_select([set(b) for b in [{f'a{i}' for i in range(4)}, {'b0', 'b1'}]], lengths, q,
                                  weight=4, min_bin_bp=0, min_comp=50, prune_margin=None)
    assert chosen_all == chosen
