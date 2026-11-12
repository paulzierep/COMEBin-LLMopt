"""CheckM2-guided final-bin selection for COMEBin-LLMopt.

COMEBin builds many Leiden partitions and by default keeps the single partition with the best
marker-gene score. This module instead picks bins across all partitions:

1. Annotate the assembly once: Prodigal (-p meta, parallel over contig chunks) and DIAMOND
   against the CheckM2 database with CheckM2's settings.
2. CheckM2 features are sums over a bin's proteins, so every candidate bin (each distinct bin of
   at least --min_bin_bp in any partition) gets CheckM2 completeness / contamination by summing
   per-contig features and applying CheckM2's own models.
3. Greedy: take the bin with the best completeness - weight * contamination (ties: larger bin),
   remove its contigs from all other candidates, rescore those, repeat while the best bin has
   completeness >= --min_comp.
4. Fill: contigs not yet used keep their bin from the base partition (the one COMEBin's marker
   score picked), so bins CheckM2 cannot score (viruses, plasmids) are kept.

Runs in a Python environment with CheckM2 installed; COMEBin calls it as a subprocess.

usage: checkm2_select.py --contigs F --ids IDS.txt --results LIST.txt --base BASE.labels.npy
                         --out OUT.labels.npy --cache DIR [--threads N] [--weight 4] [--db PATH]
Writes OUT.labels.npy: one int32 label per contig in IDS order. Contigs left unbinned get
labels of their own, so the usual minimum-bin-size filter decides about them.
"""
import argparse
import heapq
import logging
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

AA = 'ACDEFGHIKLMNPQRSTVWY'


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)


def greedy_select(bins, lengths, quality_many, weight=4.0, min_bin_bp=200000, min_comp=50.0, prune_margin=5.0):
    """Pick non-overlapping bins from overlapping candidates.

    bins: list of sets of contig IDs (modified in place); lengths: dict contig -> bp;
    quality_many: function(list of sets) -> list of (completeness, contamination) in percent.
    Returns the selected bins (list of sets) and the set of used contigs.
    Removing contigs from a bin does not raise its completeness, so a bin scored below
    min_comp - prune_margin is dropped and never rescored (prune_margin=None keeps all).
    """
    member = defaultdict(list)
    for i, b in enumerate(bins):
        for c in b:
            member[c].append(i)
    version = [0] * len(bins)
    heap = []

    def push(idx):
        for i, (comp, cont) in zip(idx, quality_many([bins[i] for i in idx])):
            if prune_margin is not None and comp < min_comp - prune_margin:
                bins[i] = set()
                continue
            size = sum(lengths[c] for c in bins[i])
            heapq.heappush(heap, (-(comp - weight * cont), -size, i, comp, version[i]))

    push(list(range(len(bins))))
    used, chosen = set(), []
    while heap:
        _, _, i, comp, ver = heapq.heappop(heap)
        if ver != version[i] or not bins[i]:
            continue
        if comp < min_comp:
            break
        chosen.append(set(bins[i]))
        used |= bins[i]
        touched = {j for c in bins[i] for j in member[c]}
        bins[i] = set()
        rescore = []
        for j in touched:
            if not bins[j]:
                continue
            bins[j] -= used
            version[j] += 1
            if bins[j] and sum(lengths[c] for c in bins[j]) >= min_bin_bp:
                rescore.append(j)
            else:
                bins[j] = set()
        if rescore:
            push(rescore)
    return chosen, used


def final_labels(contig_ids, chosen, used, base_labels):
    """Selected bins first, then base-partition bins of the unused contigs, then singletons."""
    labels = np.full(len(contig_ids), -1, dtype=np.int64)
    index = {c: k for k, c in enumerate(contig_ids)}
    for k, b in enumerate(chosen):
        for c in b:
            labels[index[c]] = k
    offset = len(chosen)
    base = np.asarray(base_labels, dtype=np.int64)
    free = labels < 0
    _, remap = np.unique(base[free], return_inverse=True)
    labels[free] = offset + remap
    return labels


class Annotation:
    """Per-contig CheckM2 base features from one annotation of the assembly."""

    def __init__(self, contigs, cache, threads, db=None):
        from checkm2 import keggData
        from checkm2.defaultValues import DefaultValues
        from checkm2.fileManager import DiamondDB
        cache = Path(cache)
        cache.mkdir(parents=True, exist_ok=True)
        feats = cache / 'contig_features.npz'
        if not feats.exists():
            faa = cache / 'proteins.faa'
            if not faa.exists():
                self._prodigal(contigs, cache, faa, threads)
            hits = cache / 'diamond.tsv'
            if not hits.exists():
                tmp = cache / 'diamond.tsv.tmp'
                subprocess.run(
                    ['diamond', 'blastp', '--outfmt', '6', '--max-target-seqs', '1', '--query', str(faa),
                     '-o', str(tmp), '--threads', str(threads), '--db', db or DiamondDB().get_DB_location(),
                     '--query-cover', str(DefaultValues.DIAMOND_QUERY_COVER),
                     '--subject-cover', str(DefaultValues.DIAMOND_SUBJECT_COVER),
                     '--id', str(DefaultValues.DIAMOND_PERCENT_ID), '--evalue', str(DefaultValues.DIAMOND_EVALUE),
                     '--block-size', '2', '--tmpdir', os.environ.get('TMPDIR', '/tmp'), '--quiet'], check=True)
                tmp.rename(hits)
            kos = list(keggData.KeggCalculator().return_default_values_from_category('KO_Genes').keys())
            self._features(faa, hits, kos, feats)
        d = np.load(feats, allow_pickle=False)
        self.index = {c: i for i, c in enumerate(d['contigs'].tolist())}
        self.meta, self.ko, self.kos = d['meta'], d['ko'], d['kos'].tolist()

    @staticmethod
    def _prodigal(contigs, cache, faa, threads):
        # prodigal -p meta treats every contig on its own, so chunks give the same genes.
        chunks = [cache / f'chunk_{i}.fa' for i in range(max(1, threads))]
        handles = [open(c, 'w') for c in chunks]
        sizes = [0] * len(chunks)
        for name, seq in sorted(read_fasta(contigs), key=lambda x: -len(x[1])):
            k = sizes.index(min(sizes))
            sizes[k] += len(seq)
            handles[k].write(f'>{name}\n{seq}\n')
        for h in handles:
            h.close()
        procs = [subprocess.Popen(['prodigal', '-p', 'meta', '-q', '-m', '-g', '11', '-a', f'{c}.faa',
                                   '-i', str(c), '-o', os.devnull, '-f', 'gff']) for c in chunks]
        failed = [p.wait() for p in procs]
        if any(failed):
            raise RuntimeError('prodigal failed')
        tmp = Path(f'{faa}.tmp')
        with open(tmp, 'w') as out:
            for c in chunks:
                with open(f'{c}.faa') as fh:
                    out.write(fh.read())
                os.remove(f'{c}.faa')
                os.remove(c)
        tmp.rename(faa)

    @staticmethod
    def _features(faa, hits, kos, out):
        ko_index = {k: i for i, k in enumerate(kos)}
        ko_of = {}
        with open(hits) as fh:
            for line in fh:
                f = line.split('\t')
                ko = f[1].split('~', 1)[1] if '~' in f[1] else None
                if ko in ko_index:
                    ko_of[f[0]] = ko_index[ko]
        meta = defaultdict(lambda: np.zeros(22))
        kcount = defaultdict(lambda: defaultdict(int))
        for pid, seq in read_fasta(faa):
            c = pid.rsplit('_', 1)[0]
            v = meta[c]
            for j, a in enumerate(AA):
                v[j] += seq.count(a)
            v[20] += len(seq)
            v[21] += 1
            if pid in ko_of:
                kcount[c][ko_of[pid]] += 1
        names = sorted(meta)
        K = np.zeros((len(names), len(kos)), dtype=np.float32)
        for i, c in enumerate(names):
            for j, n in kcount[c].items():
                K[i, j] = n
        tmp = Path(f'{out}.tmp.npz')
        np.savez_compressed(tmp, contigs=np.array(names), meta=np.array([meta[c] for c in names]),
                            ko=K, kos=np.array(kos))
        tmp.rename(out)


class Scorer:
    """CheckM2 completeness / contamination for many bins (CheckM2 'auto' model choice)."""

    def __init__(self, annotation, threads):
        self.n_scored = 0
        import pandas as pd
        from checkm2 import keggData, modelPostprocessing, modelProcessing
        self.pd, self.a = pd, annotation
        self.kegg = keggData.KeggCalculator()
        self.model = modelProcessing.modelProcessor(threads)
        self.post = modelPostprocessing.modelProcessor(threads)

    def _vectors(self, bins):
        a, pd = self.a, self.pd
        meta = np.zeros((len(bins), 22))
        ko = np.zeros((len(bins), len(a.kos)))
        for i, b in enumerate(bins):
            idx = [a.index[c] for c in b if c in a.index]
            if idx:
                meta[i] = a.meta[idx].sum(0)
                ko[i] = a.ko[idx].sum(0)
        ko_genes = pd.DataFrame(ko, columns=a.kos)
        ko_genes['Name'] = [str(i) for i in range(len(bins))]
        feats = pd.concat([ko_genes,
                           self.kegg.calculate_KO_group('KO_Pathways', ko_genes.copy()),
                           self.kegg.calculate_module_completeness(ko_genes.copy()),
                           self.kegg.calculate_KO_group('KO_Categories', ko_genes.copy())], axis=1)
        del feats['Name']
        return np.hstack([meta, feats.values.astype(float)])

    def __call__(self, bins, chunk=500):
        self.n_scored += len(bins)
        out = []
        for s in range(0, len(bins), chunk):
            X = self._vectors(bins[s:s + chunk])
            g_comp, g_cont = self.model.run_prediction_general(X)
            s_comp, scaled = self.model.run_prediction_specific(X, len(self.a.kos) + 22)
            comp, cont, _, _ = self.post.calculate_general_specific_ratio(X[:, 20], scaled, g_comp, g_cont, s_comp)
            out += list(zip(np.asarray(comp, float).tolist(), np.asarray(cont, float).tolist()))
        return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--contigs', required=True)
    p.add_argument('--ids', required=True, help='contig IDs, one per line, in label order')
    p.add_argument('--results', required=True, help='Leiden *.labels.npy files, one path per line')
    p.add_argument('--base', required=True, help='labels.npy of the partition picked by marker score')
    p.add_argument('--out', required=True)
    p.add_argument('--cache', required=True, help='directory for the annotation (reused if present)')
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--weight', type=float, default=4.0)
    p.add_argument('--min_bin_bp', type=int, default=200000)
    p.add_argument('--min_comp', type=float, default=50.0)
    p.add_argument('--prune_margin', type=float, default=5.0,
                   help='drop bins scored below min_comp minus this (negative: never drop)')
    p.add_argument('--db', help='CheckM2 DIAMOND database (default: CheckM2 configuration)')
    a = p.parse_args()
    # prodigal and diamond come with CheckM2: use the tools of this interpreter's environment first.
    os.environ['PATH'] = os.path.dirname(sys.executable) + os.pathsep + os.environ.get('PATH', '')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
    log = logging.getLogger('checkm2_select')

    ids = [line.strip() for line in open(a.ids)]
    lengths = {n: len(s) for n, s in read_fasta(a.contigs)}
    missing = [c for c in ids if c not in lengths]
    if missing:
        sys.exit(f'{len(missing)} contig IDs are not in {a.contigs}, e.g. {missing[0]}')

    log.info('Annotate the assembly for CheckM2 (cache: %s).', a.cache)
    annotation = Annotation(a.contigs, a.cache, a.threads, a.db)
    logging.getLogger().setLevel(logging.WARNING)  # CheckM2 logs every prediction batch
    scorer = Scorer(annotation, a.threads)
    logging.getLogger().setLevel(logging.INFO)

    seen, bins = set(), []
    for path in [line.strip() for line in open(a.results) if line.strip()]:
        groups = defaultdict(set)
        for lab, c in zip(np.load(path), ids):
            groups[int(lab)].add(c)
        for g in groups.values():
            key = frozenset(g)
            if key not in seen and sum(lengths[c] for c in g) >= a.min_bin_bp:
                seen.add(key)
                bins.append(set(g))
    log.info('Candidate bins (distinct, >= %d bp): %d.', a.min_bin_bp, len(bins))
    chosen, used = greedy_select(bins, lengths, scorer, a.weight, a.min_bin_bp, a.min_comp,
                                 a.prune_margin if a.prune_margin >= 0 else None)
    labels = final_labels(ids, chosen, used, np.load(a.base))
    tmp = Path(a.out).with_suffix('.tmp.npy')
    np.save(tmp, labels.astype(np.uint32))
    tmp.rename(a.out)
    log.info('Selected %d bins with CheckM2 (%d contigs, %d bin scorings); other contigs keep their '
             'base-partition bins.', len(chosen), len(used), scorer.n_scored)


if __name__ == '__main__':
    main()
