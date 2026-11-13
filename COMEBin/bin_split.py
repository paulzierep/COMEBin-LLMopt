"""Split final bins that marker genes show as contaminated (opt-in, run_comebin.sh -X).

A bin whose single-copy markers give contamination >= cont_min is re-clustered with
length-weighted k-means (k = 2..kmax) on its contigs' normalized embeddings. The split with the
best sum of (completeness - penalty * contamination) over its parts of at least the minimum
final bin size replaces the bin when it beats the unsplit bin. Smaller parts are then not
reported, like any other bin under that size.
"""
import numpy as np
from sklearn.cluster import KMeans

DEFAULT_CONT_MIN = 5.0
DEFAULT_PENALTY = 5.0
DEFAULT_KMAX = 4


def _quality(markers, marker_ids, names):
    hits = [n for n in names if n in marker_ids]
    if not hits:
        return 0.0, 0.0
    _, comp, cont = markers.bin_quality(hits)
    return comp, cont


def split_contaminated_bins(labels, contig_ids, contig_lengths, embeddings, markers, keep_labels,
                            minimum_size, cont_min=DEFAULT_CONT_MIN, penalty=DEFAULT_PENALTY,
                            kmax=DEFAULT_KMAX, seed=0, logger=None):
    """Return (new_labels, n_split). Only bins in keep_labels are considered."""
    labels = np.asarray(labels)
    lengths = np.asarray(contig_lengths, dtype=np.int64)
    emb = np.asarray(embeddings, dtype=np.float32)
    emb = emb / np.maximum(np.linalg.norm(emb, axis=1, keepdims=True), 1e-12)
    marker_ids = set(markers.bac_markers_on_contigs) | set(markers.ar_markers_on_contigs)
    out = labels.copy()  # keep the dtype (uint32 for Leiden labels)
    next_label = int(labels.max()) + 1 if len(labels) else 0
    n_split = 0

    def score(parts):
        total = 0.0
        for part in parts:
            if lengths[part].sum() < minimum_size:
                continue
            comp, cont = _quality(markers, marker_ids, [contig_ids[i] for i in part])
            total += comp - penalty * cont
        return total

    for label in sorted(int(x) for x in keep_labels):
        members = np.flatnonzero(labels == label)
        comp, cont = _quality(markers, marker_ids, [contig_ids[i] for i in members])
        if cont < cont_min or len(members) < 4:
            continue
        best_score, best_parts = score([members]), None
        for k in range(2, min(kmax, len(members) - 1) + 1):
            part_of = KMeans(k, n_init=5, random_state=seed).fit_predict(
                emb[members], sample_weight=lengths[members])
            parts = [members[part_of == t] for t in range(k)]
            s = score(parts)
            if s > best_score + 1e-9:
                best_score, best_parts = s, parts
        if best_parts is None:
            continue
        n_split += 1
        if logger is not None:
            logger.info('Split bin %d (completeness %.1f, contamination %.1f) into %d parts.',
                        label, comp, cont, len(best_parts))
        for part in best_parts[1:]:
            out[part] = next_label
            next_label += 1
    return out, n_split
