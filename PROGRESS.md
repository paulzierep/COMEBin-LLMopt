# Progress: speed and bug fixes (`perf-robustness`)

Base: `comebin-optimizations-v11` (upstream v1.1.0 plus the sequential bin-numbering fix `95f5ea8`), the version used for the runs in `paulzierep/metagenomic-binning-benchmark`.

## 2026-11-03: laptop session (CPU only, minimal data)

Test data:
- **Minimal:** `input_single` (40 contigs), the tools-iuc test set.
- **Small:** a 4,000-contig random subset of the COMEBin demo assembly; see `benchmark/README.md`.

### Done

| Commit | What | Evidence |
|---|---|---|
| Faster coverage and InfoNCE, cap HNSW neighbours | vectorised bedGraph statistics; exact memory-light InfoNCE; `-m` capped at n − 1 for small assemblies | coverage identical to 1e-15, 12.9 → 7.4 s; loss/gradients/accuracy identical, 2.26 → 0.50 s per step, −1.8 GB; the 40-contig set with default `-m 100` crashed before and now runs |
| Fail on CheckM errors, resume-safe final selection | `subprocess` with checked exit codes; error instead of `sys.exit()` (status 0) on missing tables; archaeal pass reuses Prodigal genes; an incomplete `unitem_profile` is redone, not reused | rerun after deleting the profile gives an identical `comebin_res.tsv` |
| In-memory batch loader | no DataLoader workers; same RNG consumption as the seeded DataLoader | `tests/test_batch_loader.py`; first-epoch loss identical (5.9617581) and final bins identical on the minimal set |
| No truncated or stale marker intermediates | FragGeneScan writes to a temporary prefix and publishes on success; caches must be newer than their inputs; FragGeneScan stdout no longer overwrites its own `.out` file | fresh, reuse and stale (`touch contigs.fasta`) cases tested; seeds identical |
| Zero bins instead of crash without marker seeds (benchmark issue #28) | the `bin` and `get_result` stages write an empty result plus `NO_BINS.txt` and exit 0; exit 2 (contig-ID mismatch) still fails | random-sequence assembly now ends with exit 0 and the reason in `NO_BINS.txt` |
| Per-stage timing in `run_comebin.sh` | `[COMEBin] Stage '…' finished in N s (total M s)` | tools-iuc test parameters: 4 / 54 / 5 / 64 s |
| Tests | `tests/test_info_nce.py`, `tests/test_batch_loader.py`, `tests/test_gen_cov.py` (12 tests) | `python -m pytest tests` passes |

### Not finished / open

- **End-to-end old vs new run on the 4,000-contig set:** stopped at epoch 19 because the laptop ran out of free RAM (other desktop programs). The scripts are ready in `benchmark/`; run them on the VM.
- **CI:** `.github/workflows/tests.yml` exists locally but was not pushed, because the GitHub token lacks the `workflow` scope. Push it after `gh auth refresh -s workflow`.
- **Issue #28 follow-up:** the Galaxy wrapper's grep work-around for "no usable marker seed set" can be dropped once this branch is released.

### Next steps (ranked)

1. On the VM: full old vs new comparison (`benchmark/`) on the 4,000-contig set and on CAMI II human (4,900 contigs), with CheckM2 and AMBER as in the benchmark repo.
2. **Final-result selection cost:** CheckM runs Prodigal, HMM extraction and hmmsearch once per seed-k-means bin. That's about 1.5 s per bin even with gene reuse, and on small data this stage takes longer than training. Option: one pyrodigal (meta mode) + pyhmmer pass over all contigs with CheckM's marker HMMs. It changes gene calls, so add it as an option and benchmark it rather than make it the default.
3. **Newer sub-tools to benchmark as options:**
   - pyrodigal instead of FragGeneScan for marker seeds
   - pyhmmer instead of the hmmsearch subprocess
   - `mosdepth`/`samtools depth` instead of `bedtools genomecov`, since `bedtools` is now the coverage bottleneck; these count duplicates, secondary alignments and deletions differently
   - CheckM2 to rescore only the top-k Leiden candidates. Scoring all 120 candidates with CheckM2 is too expensive.
4. **Quality on CAMI human (phase 4):** start from the sweep winner (temperature 0.05, batch ≤ 512; benchmark issue #25).
