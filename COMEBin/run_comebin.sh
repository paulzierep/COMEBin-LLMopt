#!/usr/bin/env bash

##############################################################################################################################################################
# COMEBin-LLMopt: run the binning pipeline after obtaining the bam files.
# COMEBin-LLMopt is a modified version of COMEBin v1.1.0 by Ziye Wang (zwang17@fudan.edu.cn),
# changed by Paul Zierep since 2026-09 (see NOTICE.md and CHANGELOG.md).
# Report problems with COMEBin-LLMopt at https://github.com/paulzierep/COMEBin-LLMopt/issues,
# not to the COMEBin authors.
##############################################################################################################################################################
VERSION="1.1.1"
BASED_ON="COMEBin 1.1.0"

help_message () {
  echo ""
  echo "COMEBin-LLMopt version: $VERSION (based on $BASED_ON)"
  echo "Usage: comebin-llmopt [options] -a contig_file -o output_dir -p bam_file_path"
	echo "Options:"
	echo ""
	echo "  -h              display this help message and exit"
	echo "  -V              display the COMEBin-LLMopt version and exit"
	echo "  -a STR          metagenomic assembly file"
	echo "  -o STR          output directory"
	echo "  -p STR          directory with one coverage file per sample: BAM (*.bam), or a bedGraph"
	echo "                  from 'bedtools genomecov -bga' (*.bedgraph, *.bg, optionally gzipped)"
	echo "  -n INT          number of views for contrastive multiple-view learning (default=6)"
	echo "  -t INT          stage worker/thread count where explicitly supported (default=5)"
	echo "  -w INT          concurrent Leiden worker processes (default=same as -t)"
	echo "  -m INT          HNSW neighbors retained per contig (default=100)"
	echo "  -d STR          NN training device: cuda, cuda:N, mps (Apple silicon), or cpu (default=cuda)"
	echo "  -s INT          random seed for reproducible runs (default=not set)"
	echo "  -E FLOAT        HMMER sequence E-value used when marker HMMs do not define TC cutoffs (default=1e-5)"
	echo "  -l FLOAT        temperature in loss function (default=0.1 for assemblies with an N50 > 10000, default=0.15 for others)"
	echo "  -e INT          embedding size for comebin network (default=2048)"
	echo "  -c INT          embedding size for coverage network (default=2048)"
	echo "  -O INT          output (embedding) dimension of the comebin and coverage networks (default=128)"
	echo "  -C STR          coverage transform before normalization: none or log (log(1+x)) (default=none)"
	echo "  -b INT          batch size for training process (default=1024)"
	echo "  -r INT          training epochs (default=200; cosine learning-rate schedule over the epochs)"
	echo "  -L INT          only bin contigs of at least INT bp; shorter ones stay unbinned (default=1001, i.e. no extra filter)"
	echo "  -S STR          final-bin selection: partition (best Leiden partition by marker score, as COMEBin),"
	echo "                  checkm2 (best bins across all partitions by CheckM2 estimates; fails without CheckM2), or"
	echo "                  auto (checkm2 when CheckM2 is installed and works, else partition) (default=auto)"
	echo "  -D STR          CheckM2 DIAMOND database for -S checkm2 (default=CheckM2's configured database)"
	echo "  -T STR          directory for temporary files of CheckM, CheckM2, HMMER and Python (default=\$TMPDIR, else /tmp)"
	echo "  -X              split final bins whose marker genes show >= 5% contamination (default=off; -x turns it off)"
	echo "";}

stage_start=${SECONDS}
pipeline_start=${SECONDS}
report_stage_time () {
    echo "[COMEBin-LLMopt] Stage '$1' finished in $((SECONDS - stage_start)) s (total $((SECONDS - pipeline_start)) s)."
    stage_start=${SECONDS}
}

report_stage_failure () {
    local stage_name=$1
    local stage_status=$2

    if (( stage_status >= 128 )); then
        echo "${stage_name} exited with status ${stage_status} (possible signal $((stage_status - 128)))." >&2
    else
        echo "${stage_name} exited with status ${stage_status}." >&2
    fi
    exit "${stage_status}"
}


########################################################################################################
########################     LOADING IN THE PARAMETERS AND RUNNING              ########################
########################################################################################################

num_threads=5
n_views=6
temperature=
emb_szs_forcov=2048
emb_szs=2048
batch_size=1024
epochs=200
out_dim=128
cov_transform=none
device=cuda
random_seed=
leiden_workers=
max_edges=100
hmm_evalue=1e-5
min_contig_length=

# Absolute paths without GNU-only options: BioContainers use BusyBox, whose realpath
# accepts no options. abs_path e PATH requires PATH to exist; abs_path m does not.
abs_path() {
  python -c 'import os, sys
path = os.path.realpath(os.path.expanduser(sys.argv[2]))
if sys.argv[1] == "e" and not os.path.exists(path):
    sys.exit(f"{sys.argv[2]}: No such file or directory")
print(path)' "$1" "$2"
}

split_bins=
while getopts hVXxa:o:p:n:t:w:m:d:s:E:l:e:c:b:r:L:S:D:O:C:T: OPT; do
 case ${OPT} in
  h) help_message
    exit 0
    ;;
  V) echo "COMEBin-LLMopt version: $VERSION (based on $BASED_ON)"
    exit 0
    ;;
  a) contig_file=$(abs_path e "${OPTARG}") || exit 1
    ;;
  o) output_dir=$(abs_path m "${OPTARG}") || exit 1
    ;;
  p) bam_file_path=$(abs_path e "${OPTARG}") || exit 1
    ;;
  n) n_views=${OPTARG}
    ;;
  t) num_threads=${OPTARG}
    ;;
  w) leiden_workers=${OPTARG}
    ;;
  m) max_edges=${OPTARG}
    ;;
  d) device=${OPTARG}
    ;;
  s) random_seed=${OPTARG}
    ;;
  E) hmm_evalue=${OPTARG}
    ;;
  l) temperature=${OPTARG}
    ;;
  e) emb_szs=${OPTARG}
    ;;
  c) emb_szs_forcov=${OPTARG}
    ;;
  b) batch_size=${OPTARG}
    ;;
  r) epochs=${OPTARG}
    ;;
  O) out_dim=${OPTARG}
    ;;
  C) cov_transform=${OPTARG}
    ;;
  S) final_selection=${OPTARG}
    ;;
  D) checkm2_db=$(abs_path e "${OPTARG}") || exit 1
    ;;
  L) min_contig_length=${OPTARG}
    ;;
  T) tmp_dir=$(abs_path m "${OPTARG}") || exit 1
    ;;
  X) split_bins=1
    ;;
  x) split_bins=
    ;;
  \?)
#    printf "[Usage] `date '+%F %T'` -i <INPUT_FILE> -o <OUTPUT_DIR> -o <P
#RODUCT_CODE> -s <SOFTWARE_VERSION> -t <TYPE>\n" >&2
    exit 1
 esac
done

export PYTHONUNBUFFERED=1
if [ -n "${tmp_dir:-}" ]; then
    mkdir -p -- "${tmp_dir}" || { echo "Cannot create the temporary directory (-T): ${tmp_dir}" >&2; exit 1; }
    export TMPDIR="${tmp_dir}"
fi
if [ ! -d "${TMPDIR:-/tmp}" ] || [ ! -w "${TMPDIR:-/tmp}" ]; then
    echo "The temporary directory is not writable: ${TMPDIR:-/tmp} (set it with -T DIR)" >&2
    exit 1
fi
echo "Temporary files: ${TMPDIR:-/tmp}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
seed_options=()
final_selection=${final_selection:-auto}
case "${final_selection}" in
  partition|checkm2|auto) ;;
  *) echo "Unknown -S value: ${final_selection} (use partition, checkm2 or auto)" >&2; exit 1 ;;
esac
selection_options=(--final_selection "${final_selection}")
if [ -n "${split_bins:-}" ]; then
    selection_options+=(--split_bins)
fi
if [ -n "${checkm2_db:-}" ]; then
    selection_options+=(--checkm2_db "${checkm2_db}")
fi
if [ -n "${random_seed}" ]; then
    export PYTHONHASHSEED="${random_seed}"
    export CUBLAS_WORKSPACE_CONFIG="${CUBLAS_WORKSPACE_CONFIG:-:4096:8}"
    seed_options=(--seed "${random_seed}")
    echo "Random seed: ${random_seed}"
fi
echo "Leiden workers: ${leiden_workers:-same as -t (${num_threads})}"
echo "HNSW neighbors per contig: ${max_edges}"
echo "HMMER sequence E-value: ${hmm_evalue}"

# check parameter
if [ -z "${contig_file}" -o -z "${output_dir}" -o -z "${bam_file_path}" ]; then
  help_message
  exit 1
fi

script_path=$(abs_path e "${BASH_SOURCE[0]}") || exit 1
script_dir=$(dirname -- "${script_path}")

if [[ -f "${script_dir}/main.py" ]]; then
    comebin_dir="${script_dir}"
elif [[ -f "${script_dir}/COMEBin/main.py" ]]; then
    comebin_dir="${script_dir}/COMEBin"
else
    echo "Cannot locate COMEBin/main.py relative to ${script_path}." >&2
    exit 1
fi

main_py="${comebin_dir}/main.py"

# Final-bin selection: check CheckM2 now, not after hours of training, and say what will be used.
if [[ "${final_selection}" == partition ]]; then
    echo "[COMEBin-LLMopt] Final-bin selection: partition with the best marker score (-S partition)."
else
    checkm2_check_options=()
    if [ -n "${checkm2_db:-}" ]; then
        checkm2_check_options+=(--db "${checkm2_db}")
    fi
    checkm2_message=$(python "${comebin_dir}/checkm2_check.py" "${checkm2_check_options[@]}" 2>&1)
    checkm2_usable=$?
    checkm2_message=$(printf '%s\n' "${checkm2_message}" | tail -n 1)
    if [[ ${checkm2_usable} -eq 0 ]]; then
        echo "[COMEBin-LLMopt] Final-bin selection: CheckM2 (-S ${final_selection}); ${checkm2_message}."
    elif [[ "${final_selection}" == checkm2 ]]; then
        echo "-S checkm2: ${checkm2_message}. Install CheckM2 and its database (or pass -D), or use -S partition." >&2
        exit 1
    else
        {
            echo "####################################################################################"
            echo "WARNING: -S auto cannot use CheckM2: ${checkm2_message}."
            echo "The final bins will be chosen by marker score (as -S partition), which recovered"
            echo "fewer high-quality genomes than the CheckM2 selection in our benchmarks."
            echo "Install CheckM2 and its database (or pass -D) for the CheckM2 selection, or pass"
            echo "-S partition to choose the marker-score selection explicitly."
            echo "####################################################################################"
        } >&2
        final_selection=partition
        selection_options[1]=partition
        selection_note="partition with the best marker score (-S auto without CheckM2: ${checkm2_message})"
    fi
fi

# A gzip-compressed assembly (detected by its magic bytes) is decompressed once into
# the output directory; all later steps read the plain copy.
if [ "$(python -c 'import sys; print(open(sys.argv[1], "rb").read(2).hex())' "${contig_file}")" = "1f8b" ]; then
    plain_name=$(basename -- "${contig_file}")
    plain_name=${plain_name%.gz}
    plain_name=${plain_name%.gzip}
    [[ "${plain_name}" == *.* ]] || plain_name="${plain_name}.fa"
    plain_contig_file="${output_dir}/input_${plain_name}"
    if [[ ! -s "${plain_contig_file}" ]]; then
        mkdir -p "${output_dir}" || exit 1
        gzip -cd -- "${contig_file}" > "${plain_contig_file}.part" \
            && mv -- "${plain_contig_file}.part" "${plain_contig_file}" \
            || { echo "Cannot decompress ${contig_file}" >&2; rm -f -- "${plain_contig_file}.part"; exit 1; }
    fi
    echo "Decompressed ${contig_file} to ${plain_contig_file}"
    contig_file="${plain_contig_file}"
fi

# COMEBin itself only uses contigs longer than 1000 bp. A larger minimum length
# is applied by binning a filtered copy of the assembly; coverage ignores BAM
# references that are not in it.
if [ -n "${min_contig_length}" ]; then
    if ! [[ "${min_contig_length}" =~ ^[0-9]+$ ]]; then
        echo "-L needs a positive integer, got '${min_contig_length}'." >&2
        exit 1
    fi
    if (( min_contig_length > 1001 )); then
        filtered_contig_file="${output_dir}/contigs_min${min_contig_length}.fa"
        if [[ ! -s "${filtered_contig_file}" ]]; then
            mkdir -p "${output_dir}" || exit 1
            python "${comebin_dir}/scripts/filter_min_length.py" \
                "${contig_file}" "${filtered_contig_file}" "${min_contig_length}" || exit 1
        fi
        contig_file="${filtered_contig_file}"
        echo "Minimum contig length: ${min_contig_length} bp (binning ${contig_file})"
    fi
fi

temperature_options=()
if [ -n "${temperature}" ]; then
    temperature_options=(--temperature "${temperature}")
    echo "Tau(temperature): ${temperature}"
else
    echo "Tau(temperature): automatic from complete-assembly N50 in feature_manifest.json"
fi


########################################################################################################
###### Get augmentation data
########################################################################################################
folder="${output_dir}/data_augmentation"
augmentation_status=0
feature_manifest="${folder}/feature_manifest.json"
coverage_manifest="${folder}/coverage_manifest.json"

if [[ -s "${feature_manifest}" && -s "${folder}/contig_ids.txt" &&
      -s "${folder}/contig_lengths.npy" && -s "${folder}/augmentation_intervals.npy" &&
      -s "${folder}/kmer_counts.npy" && -s "${coverage_manifest}" &&
      -s "${folder}/bam_ids.txt" && -s "${folder}/coverage_mean.npy" &&
      -s "${folder}/coverage_var.npy" ]]; then
    echo "No need to run data augmentation."
else
    echo "Running data augmentation."
    python "${main_py}" generate_aug_data --contig_file "${contig_file}" \
    --out_augdata_path "${output_dir}/data_augmentation" \
    --n_views "${n_views}" --bam_file_path "${bam_file_path}" --num_threads "${num_threads}" "${seed_options[@]}"
    augmentation_status=$?
fi

if [[ ${augmentation_status} -ne 0 ]]; then
    report_stage_failure "Data augmentation" "${augmentation_status}"
fi
report_stage_time "data augmentation"

########################################################################################################
###### Get representation (training process)
########################################################################################################
train_status=0
embedding_dir="${output_dir}/comebin_res"
embedding_manifest="${embedding_dir}/embedding_manifest.json"
embedding_ids="${embedding_dir}/embedding_ids.txt"
emb_file="${embedding_dir}/embeddings.npy"
covemb_file="${embedding_dir}/covembeddings.npy"
embedding_lengths="${embedding_dir}/embedding_lengths.npy"

if [[ -s "${embedding_manifest}" && -s "${embedding_ids}" &&
      -s "${emb_file}" && -s "${covemb_file}" &&
      -s "${embedding_lengths}" ]]; then
    echo "No need to run getting representation."
else
    echo "Running getting representation."
    python "${main_py}" train --data "${output_dir}/data_augmentation" \
    --emb_szs_forcov "${emb_szs_forcov}" --out_dim "${out_dim}" --out_dim_forcov "${out_dim}" --cov_transform "${cov_transform}" \
    --batch_size "${batch_size}" --epochs "${epochs}" --emb_szs "${emb_szs}" --n_views "${n_views}" \
    --output_path "${output_dir}/comebin_res" --earlystop --addvars --vars_sqrt \
    --device "${device}" --num_threads "${num_threads}" "${temperature_options[@]}" "${seed_options[@]}"
    train_status=$?
fi


if [[ ${train_status} -ne 0 ]]; then
    report_stage_failure "Training" "${train_status}"
fi
report_stage_time "training"

########################################################################################################
#### Clustering (run Leiden-based clustering methods and get the final result)
########################################################################################################
seed_file="${contig_file}.bacar_marker.2quarter_lencutoff_1001.seed"
leiden_options=(--max_edges "${max_edges}")
if [ -n "${leiden_workers}" ]; then
    leiden_options+=(--leiden_workers "${leiden_workers}")
fi
if [ -n "${random_seed}" ]; then
    leiden_options+=(--seed "${random_seed}")
fi
python "${main_py}" bin --contig_file "${contig_file}" \
--emb_file "${emb_file}" \
--output_path "${output_dir}/comebin_res" \
--seed_file "${seed_file}" --num_threads "${num_threads}" --hmm_evalue "${hmm_evalue}" "${leiden_options[@]}"
bin_status=$?

if [[ ${bin_status} -ne 0 ]]; then
    report_stage_failure "Clustering" "${bin_status}"
fi
report_stage_time "clustering"

python "${main_py}" get_result --contig_file "${contig_file}" \
--output_path "${output_dir}/comebin_res" \
--emb_file "${emb_file}" \
--seed_file "${seed_file}" --num_threads "${num_threads}" --max_edges "${max_edges}" \
--hmm_evalue "${hmm_evalue}" "${seed_options[@]}" "${selection_options[@]}"
result_status=$?

if [[ ${result_status} -ne 0 ]]; then
    report_stage_failure "Final-result selection" "${result_status}"
fi
report_stage_time "final-result selection"

selection_file="${output_dir}/comebin_res/final_selection.txt"
if [[ -n "${selection_note:-}" ]]; then
    echo "${selection_note}" > "${selection_file}"
    echo "WARNING: final bins were chosen by ${selection_note}." >&2
elif [[ -s "${selection_file}" ]]; then
    echo "[COMEBin-LLMopt] Final bins chosen by: $(cat "${selection_file}")"
fi
