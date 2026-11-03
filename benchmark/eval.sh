#!/usr/bin/env bash
# Score the final bins of the "old" and "new" runs with CheckM v1 (Bacteria domain
# marker set; no pplacer, so it fits in a few GB of RAM) and draw the README figures.
# usage: BENCH_DIR=/path/to/bench eval.sh   (COMEBin conda environment active)
set -eu
here=$(dirname "$(realpath "$0")")
cd "${BENCH_DIR:?set BENCH_DIR}"
for v in old new; do
  rm -rf "$v/checkm_out"
  checkm taxonomy_wf -t "${THREADS:-6}" -x fa --tab_table -f "$v/checkm.tsv" domain Bacteria \
    "$v/out/comebin_res/comebin_res_bins" "$v/checkm_out" > "$v/checkm.log" 2>&1
done
python "$here/analyze.py" . "$here/figures"
