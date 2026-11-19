# Changelog

COMEBin-LLMopt versions. Each entry lists changes relative to the previous version;
1.0.0 lists the changes relative to COMEBin v1.1.0.

## Unreleased

### New options
- **`-r INT`: number of training epochs** (default 200, unchanged; the cosine learning-rate schedule spans the given epochs). 100 or 300 epochs were not better than 200 on the gold-standard datasets (CAMI II human, seed 42: genomes > 90 % complete and < 5 % contaminated 28 at 200 epochs, 31 at 100, 25 at 300; a second dataset was worse with both), so the default stays.
- **`-O INT`: output (embedding) dimension** of the k-mer/coverage and coverage networks (default 128, unchanged). 64 and 256 were not better than 128 on the gold-standard datasets (CAMI II human, seed 42: genomes > 90 % complete and < 5 % contaminated 28 at 128, 27 at 64 and at 256; a second dataset was worse with both), so the default stays.
- **`-T DIR`: directory for temporary files** of CheckM, CheckM2, HMMER and Python (exported as `TMPDIR`, and passed to CheckM as `--tmpdir`; default `$TMPDIR`, else `/tmp`). The run stops at the start if the directory is not writable. On cluster nodes with a small `/tmp`, CheckM's `hmmfetch` otherwise fails with "No space left on device" on large assemblies (reported in COMEBin issue 61).

### Changed
- **gzip-compressed assemblies are accepted** (`-a contigs.fa.gz`, detected by content, not by name). They are decompressed once to `<output>/input_<name>` and all steps read that copy (COMEBin issue 16).
- **CheckM2 is checked before training.** Previously `-S auto` found out only at the end of a run that CheckM2 or its database was missing, and then only logged a warning. Now `run_comebin.sh` checks CheckM2 and its database at the start (`checkm2_check.py`) and prints which final-bin selection it will use. Without a working CheckM2, `-S auto` prints a prominent warning on stderr at the start and again at the end; `-S checkm2` stops at once instead of after training. `comebin_res/final_selection.txt` records the selection that made the final bins (including a CheckM2 failure during the run). Bins are unchanged.

## 1.1.0 (2026-11-17)

### New default (better bins)
- **`-S auto` is the new default:** the CheckM2 bin selection (`-S checkm2`) is used when CheckM2 is installed and works, otherwise the partition chosen by marker score, with a warning in the log. `-S checkm2` now fails when CheckM2 is missing or fails, instead of silently falling back.
  - CAMI II human (AMBER, seeds 42 / 7), CheckM2 selection against 1.0.x: genomes > 90 % complete and < 5 % contaminated 20 / 21 → 25 / 24, ARI 0.77 / 0.76 → 0.82 / 0.81. Run time 32 → 36 min, peak memory of the CheckM2 step about 10 GB.
- **Temperature 0.1 instead of 0.07** for assemblies with an N50 above 10 kb (the automatic choice when `-l` is not given; 0.15 for lower N50 is unchanged).
  - CAMI II human (AMBER, seeds 42 / 7 / 13, with `-S auto`): genomes > 90 % complete and < 5 % contaminated 25 / 24 / 26 → 28 / 26 / 30, > 50 % / < 10 % 32 / 31 / 32 → 35 / 32 / 36, ARI 0.82 → 0.84–0.85, bp F1 (seed 42) 0.667 → 0.725. Fewer, larger bins (100 → 75 with seed 42). Same run time.
  - Also better on a second gold-standard dataset with all three seeds.

### New options
- **`-X`: split contaminated bins** (opt-in). After the final bins are chosen, each bin whose single-copy marker genes show at least 5 % contamination is re-clustered with length-weighted k-means (k = 2–4) on its contigs' embeddings. The split is kept when it raises the sum of (completeness − 5 × contamination) over the parts of at least 200 kb; smaller parts are not reported. No extra run time.
  - CAMI II human (AMBER, seeds 42 / 7 / 13), with the CheckM2 selection: genomes > 90 % complete and < 5 % contaminated 25 / 24 / 26 → 26 / 27 / 26, ARI 0.82 / 0.81 / 0.82 → 0.87 / 0.84 / 0.87. Assigned sequence drops from 99 % to 94–98 %.
  - It stays opt-in because it can lower the average bin completeness: the bp-weighted AMBER F1 fell from 0.728 to 0.671 with seed 7, and splitting often cuts low-quality bins into several small ones.
  - BATS medium (CheckM2, 3 seeds): high/medium-quality bins unchanged (0 / 2), mean contamination 3.5–4.2 % → 1.4–1.7 %.

## 1.0.2 (2026-11-12)

### Fixes
- **No more empty output on small assemblies.** COMEBin picks the Leiden partition with the best marker-gene score, but counts bins of any size, while only bins of at least 200 kb are reported. On small assemblies the best-scoring partition could hold only small bins, so the run finished with no bins at all. Now, only when that partition would report no bin, the choice is made among the bins that are reported, with ties going to the partition that reports the most sequence.
  - On larger data the choice is the same as before (COMEBin's rule decides whenever its partition reports a bin).
  - Test assembly of 300 contigs (CAMI II human subset, 10 % of the reads, 2 threads): seeds 7 / 13 / 1 gave 0 / 0 / 0 bins before, now 6 / 8 / 8; seed 42 gives 6 (about 7 min each).

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
