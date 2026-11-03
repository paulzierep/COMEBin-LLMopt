#!/usr/bin/env bash
# Time one COMEBin version on one dataset and sample its memory every 2 s.
#
# usage: BENCH_DIR=/path/to/bench run_one.sh <label> <comebin_checkout> <contigs.fasta> <bam_dir> [threads]
#   Writes $BENCH_DIR/<label>/{run.log,mem.tsv,out/}. The FASTA is copied into
#   $BENCH_DIR/<label>/data/ so each version generates its own marker-seed files.
#   Run versions one after another, never concurrently, so timings are comparable.
set -u
label=$1; code=$(realpath "$2"); fasta=$(realpath "$3"); bams=$(realpath "$4"); threads=${5:-6}
B=${BENCH_DIR:?set BENCH_DIR}/$label
mkdir -p "$B/data"; rm -rf "$B/out"; : > "$B/mem.tsv"
cp "$fasta" "$B/data/contigs.fasta"
( while sleep 2; do
    pids=$(pgrep -f "$code|checkm|FragGeneScan|hmmsearch|bedtools|prodigal" | paste -sd, -)
    [ -n "$pids" ] && printf '%s\t%s\n' "$(date +%s)" \
      "$(ps -o rss= -p "$pids" 2>/dev/null | awk '{s+=$1} END{print s}')" >> "$B/mem.tsv"
  done ) & sampler=$!
start=$(date +%s)
/usr/bin/time -v bash "$code/COMEBin/run_comebin.sh" -a "$B/data/contigs.fasta" -p "$bams" \
  -o "$B/out" -t "$threads" -d "${DEVICE:-cpu}" -s 1 > "$B/run.log" 2>&1
echo "exit=$? start=$start end=$(date +%s)" >> "$B/run.log"
kill "$sampler"
