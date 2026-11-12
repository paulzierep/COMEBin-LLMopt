# Changelog

COMEBin-LLMopt versions. Each entry lists changes relative to the previous version;
1.0.0 lists the changes relative to COMEBin v1.1.0.

## Unreleased

### Better bins
- **`-S checkm2`: final bins chosen across all Leiden partitions with CheckM2 estimates** (opt-in; needs CheckM2 and its database, `-D`).
  - How it works:
    - COMEBin builds about 120 Leiden partitions and by default keeps the one with the best marker-gene score.
    - The new mode annotates the assembly once (Prodigal in parallel, plus CheckM2's DIAMOND search) and gives every distinct candidate bin CheckM2 completeness and contamination. CheckM2's own models run on per-bin sums of per-contig features; on real bins this matches `checkm2 predict` at r ≥ 0.995.
    - It then greedily picks the best bins (completeness − 4 × contamination; ties go to the larger bin), removing the picked contigs from the other candidates.
    - Contigs left over keep their bins from the marker-chosen partition, so viruses and plasmids that CheckM2 cannot score are kept.
  - CAMI II human, genomes recovered at > 90 % completeness / < 5 % contamination (gold standard, seeds 42 / 7 / 13): 20 / 21 / 21 → 25 / 24 / 26. AMBER ARI goes up (0.768 → 0.818 for seed 42); F1 (bp) is 0.01 lower, within seed noise.
  - Cost: about 11 % more run time on CAMI II human (36 instead of 32 min). On a 1.1 Gbp assembly the annotation takes about 12 min (32 threads), plus a few minutes of selection.

### Faster
- **Dropout:** the keep-rate scale is folded into the mask, so dropout costs one multiply instead of two in the forward and backward pass (about 3% faster step; outputs and gradients are bit-identical to 1.0.1).

## 1.0.1 (2026-11-11)

Same model, loss, clustering and bin selection as 1.0.0; only the dropout masks are drawn differently.

### Faster
- **Dropout:** on CPU, masks come from random bytes (7 per 64-bit draw) instead of `rand_like >= p`: 4.5× faster masks, 17% faster training step at 32 threads. The drop rate is rounded to 1/256 (p = 0.2 drops 0.199) and the output is scaled by the exact keep rate. With the same seed, results differ from 1.0.0 as with a different seed.
- **End to end:** CAMI II human takes 32–34 min instead of 38–39 min (seeds 42 and 7) and BATS medium about 14 min. Bin quality is unchanged within seed noise: AMBER genomes >90 % complete / <5 % contaminated 20–21 vs 21 and F1 (bp) 0.678/0.737 vs 0.678/0.750.

## 1.0.0 (2026-11-11)

Based on COMEBin v1.1.0.

### Renamed
- The tool is now **COMEBin-LLMopt**. Run it with `comebin-llmopt` (same options as `run_comebin.sh`, which still works). Modification notice in `NOTICE.md`.

### Faster
- **Training loss:** InfoNCE takes one log-sum-exp over the negatives per anchor instead of building the full logits matrix. The loss step is 4.5× faster and peaks 1.8 GB lower; loss, gradients and accuracy are unchanged.
- **Training batches:** batches are gathered from the in-memory views; no DataLoader worker processes are forked each epoch. Seeded batches are identical to v1.1.0.
- **Coverage features:** bedGraph statistics are vectorised, 1.7× faster.
- **Dropout:** masks are drawn with `rand_like >= p` instead of `bernoulli_` (same distribution, ~1.4× faster mask).
- **Final-bin selection:** the archaeal CheckM pass reuses the bacterial pass's gene calls.
- **End to end:** BATS medium (3,000 contigs, 32 CPU threads) takes 19.5 min instead of 38.5 min with near-identical bins (ARI 0.936; the difference is floating-point rounding).

### Added
- `-L INT`: bin only contigs of at least INT bp (default 1001, no extra filter); shorter contigs stay unbinned. For large assemblies (about 100k contigs), `-L 2500` keeps about a third of them. The automatic temperature then follows the N50 of the filtered contigs.

### Fixed
- Assemblies with no usable marker-gene seeds end with exit 0, an empty result and `NO_BINS.txt` instead of an error.
- Small assemblies (fewer contigs than `-m`) no longer crash; the neighbour count is capped at n − 1.
- CheckM failures stop the run instead of being ignored.
- Interrupted runs no longer reuse incomplete CheckM profiles, truncated FragGeneScan output or marker files older than the contigs.
- DataLoader `pymp-*` "Device or resource busy" errors on network file systems no longer occur.

### Added
- Per-stage wall time in the log.
- Unit tests in `tests/` (`python -m pytest tests`).
