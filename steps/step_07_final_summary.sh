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

############################################
# STATUS FILE
############################################

STATUS_FILE="${ROOT_DIR}/results/pipeline_status.tsv"
init_status_file

export STATUS_FILE
export UTILS_FILE

############################################
# DIRECTORIES
############################################

FINAL_SUMMARY_DIR="${ROOT_DIR}/${FINAL_SUMMARY_DIR}"
FINAL_SUMMARY_LOG_DIR="${ROOT_DIR}/${LOG_DIR}/final_summary"

mkdir -p "${FINAL_SUMMARY_DIR}"
mkdir -p "${FINAL_SUMMARY_LOG_DIR}"

OUT_TSV="${FINAL_SUMMARY_DIR}/final_master_table.tsv"
OUT_XLSX="${FINAL_SUMMARY_DIR}/final_master_table.xlsx"

READS_QC_DIR="${ROOT_DIR}/results/reads/qc/nanostat_clean"
QUAST_DIR="${ROOT_DIR}/${QUAST_INDIVIDUAL_DIR}"
CHECKM2_FILE="${ROOT_DIR}/${CHECKM2_SUMMARY_DIR}/checkm2_summary.tsv"
TAXONOMY_FILE="${ROOT_DIR}/${TAXONOMY_MERGED_DIR}/final_species_table.tsv"
MLST_FILE="${ROOT_DIR}/${MLST_SUMMARY_DIR}/complete_mlst_table.tsv"
BAKTA_DIR="${ROOT_DIR}/${BAKTA_DIR}"
AMR_FILE="${ROOT_DIR}/${AMRFINDER_SUMMARY_DIR}/amrfinder.result.combined.tsv"
PLASMIDFINDER_FILE="${ROOT_DIR}/${PLASMIDFINDER_SUMMARY_DIR}/plasmidfinder_inc_types.tsv"
MOBTYPER_FILE="${ROOT_DIR}/${MOBSUITE_SUMMARY_DIR}/mobtyper.result.combined.tsv"

############################################
# CHECK ENV
############################################

if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_ANALYSIS}"; then
    echo "ERROR: conda environment not found: ${ENV_ANALYSIS}"
    exit 1
fi

############################################
# CHECK INPUT FILES
############################################

for f in \
    "${CHECKM2_FILE}" \
    "${TAXONOMY_FILE}" \
    "${MLST_FILE}" \
    "${AMR_FILE}" \
    "${PLASMIDFINDER_FILE}" \
    "${MOBTYPER_FILE}"
do
    if [[ ! -f "${f}" ]]; then
        log_error "Missing input file"
        status_fail "GLOBAL" "final_summary" "missing_input"
        exit 1
    fi
done

if [[ ! -d "${READS_QC_DIR}" ]]; then
    echo "ERROR: Reads nanostats directory not found: ${READS_QC_DIR}"
    exit 1
fi

if [[ ! -d "${QUAST_DIR}" ]]; then
    echo "ERROR: QUAST directory not found: ${QUAST_DIR}"
    exit 1
fi

if [[ ! -d "${BAKTA_DIR}" ]]; then
    echo "ERROR: Bakta directory not found: ${BAKTA_DIR}"
    exit 1
fi

for sample in "${CHECKM2_FILE}" "${TAXONOMY_FILE}" "${MLST_FILE}" "${AMR_FILE}" "${PLASMIDFINDER_FILE}" "${MOBTYPER_FILE}"; do
    if [[ ! -f "${sample}" ]]; then
        log_error "Missing critical input file: $sample"
        status_start "GLOBAL" "final_summary"
        status_fail "GLOBAL" "final_summary" "missing_input_$(basename "$sample")"
        exit 1
    fi
done

############################################
# RUN MERGE
############################################

echo "=================================================="
echo "STEP 07: FINAL SUMMARY"
echo "Samples file         : ${SAMPLES}"
echo "Reads QC file        : ${READS_QC_DIR}"
echo "QUAST directory      : ${QUAST_DIR}"
echo "CheckM2 file         : ${CHECKM2_FILE}"
echo "Taxonomy file        : ${TAXONOMY_FILE}"
echo "MLST file            : ${MLST_FILE}"
echo "Bakta directory      : ${BAKTA_DIR}"
echo "AMR file             : ${AMR_FILE}"
echo "PlasmidFinder file   : ${PLASMIDFINDER_FILE}"
echo "MOB-typer file       : ${MOBTYPER_FILE}"
echo "Output TSV           : ${OUT_TSV}"
echo "Output XLSX          : ${OUT_XLSX}"
echo "ENV_ANALYSIS         : ${ENV_ANALYSIS}"
echo "=================================================="

#Status global
status_start "GLOBAL" "final_summary"

set +e
conda run --no-capture-output -n "${ENV_ANALYSIS}" python3 \
    "${SCRIPT_MERGE_FINAL_REPORT}" \
    --samples "${SAMPLES}" \
    --reads-qc "${READS_QC_DIR}" \
    --quast "${QUAST_DIR}" \
    --checkm2 "${CHECKM2_FILE}" \
    --taxonomy "${TAXONOMY_FILE}" \
    --mlst "${MLST_FILE}" \
    --bakta "${BAKTA_DIR}" \
    --amr "${AMR_FILE}" \
    --plasmidfinder "${PLASMIDFINDER_FILE}" \
    --mobtyper "${MOBTYPER_FILE}" \
    --out "${OUT_TSV}" \
    --out-xlsx "${OUT_XLSX}" \
    > "${FINAL_SUMMARY_LOG_DIR}/merge_final_report.stdout.log" \
    2> "${FINAL_SUMMARY_LOG_DIR}/merge_final_report.stderr.log"
status=$?
set -e

if [[ $status -ne 0 ]]; then
    status_fail "GLOBAL" "final_summary" "merged_failed"
    exit 1
fi

#Check output
if ! check_file_not_empty "${OUT_TSV}"; then
    status_fail "GLOBAL" "final_summary" "empty_output"
    exit 1
fi

#Final log info
log_info "final_summary step completed for all samples"
status_ok "GLOBAL" "final_summary"

echo "[DONE] Final summary written to: ${FINAL_SUMMARY_DIR}"