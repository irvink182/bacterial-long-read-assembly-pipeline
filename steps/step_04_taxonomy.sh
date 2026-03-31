#!/usr/bin/env bash
set -euo pipefail

############################################
# PATHS AND CONFIG
############################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

CONFIG_FILE="${ROOT_DIR}/config/config.sh"
DB_FILE="${ROOT_DIR}/config/databases.config"
UTILS_FILE="${ROOT_DIR}/scripts/utils.sh"

[[ -f "${CONFIG_FILE}" ]] || { echo "[ERROR] Config file not found: ${CONFIG_FILE}" >&2; exit 1; }
[[ -f "${DB_FILE}" ]] || { echo "[ERROR] Database config file not found: ${DB_FILE}" >&2; exit 1; }
[[ -f "${UTILS_FILE}" ]] || { echo "[ERROR] Utils config file not found: ${UTILS_FILE}" >&2; exit 1; }

source "${CONFIG_FILE}"
source "${DB_FILE}"
source "${UTILS_FILE}"

############################################
# ARGUMENTS
############################################

SAMPLES=""
TRIMMER=""
RESUME="false"
FORCE="false"
export RESUME
export FORCE

# =========================
# Parse arguments
# =========================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --samples)
            SAMPLES="$2"
            shift 2
            ;;
        --trimmer)
            TRIMMER="$2"
            shift 2
            ;;
        --resume)
            RESUME="true"
            shift
            ;;
        --force)
            FORCE="true"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "${SAMPLES}" ]]; then
    die "missing --samples"
fi

require_file "${SAMPLES}" "samples file"

# =========================
# Info
# =========================
log_info "Starting taxonomy step"
[[ "${RESUME}" == "true" ]] && log_info "Resume mode enabled"

############################################
# DIRECTORIES
############################################

ASM_DIR="${ROOT_DIR}/${FINAL_ASM_DIR}"

SOURMASH_SIG_DIR="${ROOT_DIR}/${SOURMASH_SIG_DIR}"
SOURMASH_SEARCH_DIR="${ROOT_DIR}/${SOURMASH_SEARCH_DIR}"

SKANI_SEARCH_DIR="${ROOT_DIR}/${SKANI_SEARCH_DIR}"

MLST_INDIVIDUAL_DIR="${ROOT_DIR}/${MLST_INDIVIDUAL_DIR}"
MLST_SUMMARY_DIR="${ROOT_DIR}/${MLST_SUMMARY_DIR}"

TAXONOMY_MERGED_DIR="${ROOT_DIR}/${TAXONOMY_MERGED_DIR}"

# Inputs for merge
SYLPH_TABLE_DIR="${ROOT_DIR}/${SYLPH_TABLE_DIR}"
CHECKM2_INDIVIDUAL_DIR="${ROOT_DIR}/${CHECKM2_INDIVIDUAL_DIR}"

TAXONOMY_LOG_DIR="${ROOT_DIR}/${LOG_DIR}/taxonomy"

mkdir -p "${SOURMASH_SIG_DIR}"
mkdir -p "${SOURMASH_SEARCH_DIR}"
mkdir -p "${SKANI_SEARCH_DIR}"
mkdir -p "${MLST_INDIVIDUAL_DIR}"
mkdir -p "${MLST_SUMMARY_DIR}"
mkdir -p "${TAXONOMY_MERGED_DIR}"
mkdir -p "${TAXONOMY_LOG_DIR}"

############################################
# CHECK CONDA ENVS
############################################

if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_TAXONOMY}"; then
    echo "ERROR: conda environment not found: ${ENV_TAXONOMY}"
    exit 1
fi

if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_ANALYSIS}"; then
    echo "ERROR: conda environment not found: ${ENV_ANALYSIS}"
    exit 1
fi

############################################
# PARALLEL SETTINGS
############################################
TOTAL_CORES="${TOTAL_CORES:-32}"
THREADS_PER_JOB="${THREADS_PER_JOB:-8}"

PARALLEL_JOBS=$(( TOTAL_CORES / THREADS_PER_JOB ))
if [[ "${PARALLEL_JOBS}" -lt 1 ]]; then
    PARALLEL_JOBS=1
fi

# Keep explicit jobs if user defined them in config
SKETCH_JOBS="${SKETCH_JOBS:-${PARALLEL_JOBS}}"
SMASH_JOBS="${SMASH_JOBS:-${PARALLEL_JOBS}}"

SKANI_THREADS="${SKANI_THREADS:-${THREADS_PER_JOB}}"
SKANI_JOBS="${SKANI_JOBS:-${PARALLEL_JOBS}}"

MLST_JOBS="${MLST_JOBS:-${PARALLEL_JOBS}}"

############################################
# FUNCTIONS
############################################

resolve_asm() {
    local sample="$1"

    if [[ -s "${ASM_DIR}/${sample}.${FINAL_ASM_PREFIX}.fasta" ]]; then
        echo "${ASM_DIR}/${sample}.${FINAL_ASM_PREFIX}.fasta"
        return
    fi
    if [[ -s "${ASM_DIR}/${sample}.${FINAL_ASM_PREFIX}.fa" ]]; then
        echo "${ASM_DIR}/${sample}.${FINAL_ASM_PREFIX}.fa"
        return
    fi
    if [[ -s "${ASM_DIR}/${sample}.${FINAL_ASM_PREFIX}.fna" ]]; then
        echo "${ASM_DIR}/${sample}.${FINAL_ASM_PREFIX}.fna"
        return
    fi

    echo ""
}

export ENV_TAXONOMY ENV_ANALYSIS
export ASM_DIR FINAL_ASM_PREFIX
export SOURMASH_DB SKANI_DB
export SOURMASH_SIG_DIR SOURMASH_SEARCH_DIR SKANI_SEARCH_DIR MLST_INDIVIDUAL_DIR
export SOURMASH_KSIZE SOURMASH_SCALED SKANI_THREADS SKANI_TOPN SKANI_BOTH_MIN_AF
export TAXONOMY_LOG_DIR
export -f resolve_asm


############################################
# MAIN
############################################

echo "=================================================="
echo "STEP 04: TAXONOMY"
echo "Samples file        : ${SAMPLES}"
echo "Assembly dir        : ${ASM_DIR}"
echo "Sourmash sig dir    : ${SOURMASH_SIG_DIR}"
echo "Sourmash search dir : ${SOURMASH_SEARCH_DIR}"
echo "Skani search dir    : ${SKANI_SEARCH_DIR}"
echo "MLST dir            : ${MLST_INDIVIDUAL_DIR}"
echo "Merged dir          : ${TAXONOMY_MERGED_DIR}"
echo "ENV_TAXONOMY        : ${ENV_TAXONOMY}"
echo "ENV_ANALYSIS        : ${ENV_ANALYSIS}"
echo "PARALLEL_JOBS       : ${PARALLEL_JOBS}"
echo "=================================================="

echo "Checking first 5 samples resolve..."
head -n 6 "${SAMPLES}" | tail -n +2 | cut -f1 | while read -r s; do
    a=$(resolve_asm "$s")
    [[ -n "$a" ]] && echo "OK $s -> $a" || echo "MISSING $s"
done

############################################
# 1) SOURMASH SKETCH
############################################

echo "[1/4] sourmash sketch..."

if [[ "${RESUME}" == "true" && -n "$(ls -A ${SOURMASH_SIG_DIR}/*.sig.gz 2>/dev/null)" ]]; then
    log_info "Skipping sourmash sketch: signatures already exist"
else

tail -n +2 "${SAMPLES}" | \
parallel --colsep '\t' -j "${SKETCH_JOBS}" --halt now,fail=1 --linebuffer \
'
sample={1}
gsize={3}

asm=$(resolve_asm "$sample")
if [[ -z "$asm" ]]; then
    echo "Missing assembly for sample: $sample" >&2
    exit 2
fi

out_sig="'"${SOURMASH_SIG_DIR}"'/${sample}.sig.gz"
log="'"${TAXONOMY_LOG_DIR}"'/${sample}.sourmash.sketch.log"

echo "sample=$sample expected_genome_size=$gsize asm=$asm" > "$log"

conda run --no-capture-output -n "'"${ENV_TAXONOMY}"'" \
    sourmash sketch dna \
    -p k='"${SOURMASH_KSIZE}"',scaled='"${SOURMASH_SCALED}"',noabund \
    -o "$out_sig" \
    "$asm" >> "$log" 2>&1
'
fi

############################################
# 2) SOURMASH SEARCH
############################################

echo "[2/4] sourmash search..."

if [[ "${RESUME}" == "true" && -n "$(ls -A ${SOURMASH_SEARCH_DIR}/*.sourmash.search.csv 2>/dev/null)" ]]; then
    log_info "Skipping sourmash search: results already exist"
else

tail -n +2 "${SAMPLES}" | \
parallel --colsep '\t' -j "${SMASH_JOBS}" --halt now,fail=1 --linebuffer \
'
sample={1}
gsize={3}

sig="'"${SOURMASH_SIG_DIR}"'/${sample}.sig.gz"
out_csv="'"${SOURMASH_SEARCH_DIR}"'/${sample}.sourmash.search.csv"
log="'"${TAXONOMY_LOG_DIR}"'/${sample}.sourmash.search.log"

if [[ ! -s "$sig" ]]; then
    echo "Missing sig: $sig" >&2
    exit 2
fi

echo "sample=$sample expected_genome_size=$gsize sig=$sig" > "$log"

conda run --no-capture-output -n "'"${ENV_TAXONOMY}"'" \
    sourmash search \
    "$sig" "'"${SOURMASH_DB}"'" \
    --output "$out_csv" >> "$log" 2>&1
'
fi

############################################
# 3) SKANI SEARCH
############################################

echo "[3/4] skani search..."

if [[ "${RESUME}" == "true" && -n "$(ls -A ${SKANI_SEARCH_DIR}/*.skani.search.tsv 2>/dev/null)" ]]; then
    log_info "Skipping skani search: results already exist"
else

tail -n +2 "${SAMPLES}" | \
parallel --colsep '\t' -j "${SKANI_JOBS}" --halt now,fail=1 --linebuffer \
'
sample={1}
gsize={3}

asm=$(resolve_asm "$sample")
if [[ -z "$asm" ]]; then
    echo "Missing assembly for sample: $sample" >&2
    exit 2
fi

out_tsv="'"${SKANI_SEARCH_DIR}"'/${sample}.skani.search.tsv"
log="'"${TAXONOMY_LOG_DIR}"'/${sample}.skani.search.log"

echo "sample=$sample expected_genome_size=$gsize asm=$asm threads='"${SKANI_THREADS}"'" > "$log"

conda run --no-capture-output -n "'"${ENV_TAXONOMY}"'" \
    skani search \
    "$asm" \
    -d "'"${SKANI_DB}"'" \
    --both-min-af "'"${SKANI_BOTH_MIN_AF}"'" \
    -n "'"${SKANI_TOPN}"'" \
    -t "'"${SKANI_THREADS}"'" \
    -o "$out_tsv" >> "$log" 2>&1
'
fi

############################################
# 4) MLST
############################################

echo "[4/4] MLST..."

if [[ "${RESUME}" == "true" && -n "$(ls -A ${MLST_INDIVIDUAL_DIR}/*.mlst.tsv 2>/dev/null)" ]]; then
    log_info "Skipping MLST: results already exist"
else

tail -n +2 "${SAMPLES}" | \
parallel --colsep '\t' -j "${MLST_JOBS}" --halt now,fail=1 --linebuffer \
'
sample={1}
gsize={3}

asm=$(resolve_asm "$sample")
if [[ -z "$asm" ]]; then
    echo "Missing assembly for sample: $sample" >&2
    exit 2
fi

out_tsv="'"${MLST_INDIVIDUAL_DIR}"'/${sample}.mlst.tsv"
log="'"${TAXONOMY_LOG_DIR}"'/${sample}.mlst.log"

echo "sample=$sample expected_genome_size=$gsize asm=$asm" > "$log"

conda run --no-capture-output -n "'"${ENV_TAXONOMY}"'" \
    mlst \
    --label "$sample" \
    "$asm" > "$out_tsv" 2>> "$log"
'
fi

############################################
# MLST SUMMARY
############################################

echo "[INFO] Creating MLST summary table"

cd "${MLST_INDIVIDUAL_DIR}"

if compgen -G "*.mlst.tsv" > /dev/null; then
    cat *.mlst.tsv | awk 'BEGIN {print "Sample\tScheme\tST\tGeneA\tGeneB\tGeneC\tGeneD\tGeneE\tGeneF\tGeneG"}1' \
        > "${MLST_SUMMARY_DIR}/complete_mlst_table.tsv"
fi

cd "${ROOT_DIR}"

############################################
# 5) MERGE FINAL TAXONOMY TABLE
############################################

echo "[5/5] merge_taxonomy_reports.py..."

if should_skip_global "${FINAL_SPECIES_TABLE}"; then
    log_info "Skipping taxonomy merge"
else

conda run --no-capture-output -n "${ENV_ANALYSIS}" python3 "${MERGE_TAXONOMY_BIN}"\
    --asm-type "${ASM_TYPE_FOR_MERGE}" \
    --sample-list "${SAMPLES}" \
    --skani-files "${SKANI_SEARCH_DIR}/*.skani.search.tsv" \
    --sylph-files "${SYLPH_TABLE_DIR}/*.sylph_table.tsv" \
    --sourmash-files "${SOURMASH_SEARCH_DIR}/*.sourmash.search.csv" \
    --checkm2-files "${CHECKM2_INDIVIDUAL_DIR}/*/*.checkm2.quality_report.tsv" \
    --out "${FINAL_SPECIES_TABLE}" \
    > "${TAXONOMY_LOG_DIR}/merge_taxonomy.stdout.log" \
    2> "${TAXONOMY_LOG_DIR}/merge_taxonomy.stderr.log"
fi

echo "[INFO] Taxonomy module completed"
echo "[DONE] Outputs written to: ${ROOT_DIR}/${TAXONOMY_RESULTS_DIR}"