#!/usr/bin/env bash
set -euo pipefail

############################################
# PATHS AND CONFIG
############################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

CONFIG_FILE="${ROOT_DIR}/config/config.sh"
DB_FILE="${ROOT_DIR}/config/databases.config"


[[ -f "${CONFIG_FILE}" ]] || { echo "[ERROR] Config file not found: ${CONFIG_FILE}" >&2; exit 1; }
[[ -f "${DB_FILE}" ]] || { echo "[ERROR] Database config file not found: ${DB_FILE}" >&2; exit 1; }

source "${CONFIG_FILE}"
source "${DB_FILE}"
source "${ROOT_DIR}/scripts/utils.sh"

############################################
# ARGUMENTS
############################################

SAMPLES=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --samples)
            SAMPLES="$2"
            shift 2
            ;;
        *)
            echo "ERROR: unknown argument: $1"
            exit 1
            ;;
    esac
done

if [[ -z "${SAMPLES}" ]]; then
    die "missing --samples"
fi

require_file "${SAMPLES}" "samples file"


############################################
# DIRECTORIES
############################################

FINAL_SUMMARY_DIR="${ROOT_DIR}/${FINAL_SUMMARY_DIR}"
FINAL_SUMMARY_LOG_DIR="${ROOT_DIR}/${LOG_DIR}/final_summary"

mkdir -p "${FINAL_SUMMARY_DIR}"
mkdir -p "${FINAL_SUMMARY_LOG_DIR}"

OUT_TSV="${FINAL_SUMMARY_DIR}/final_master_table.tsv"
OUT_XLSX="${FINAL_SUMMARY_DIR}/final_master_table.xlsx"

READS_QC_FILE="${ROOT_DIR}/results/reads/qc/multiqc/clean/multiqc_nanostats_clean_report.txt"
QUAST_FILE="${ROOT_DIR}/${QUAST_SUMMARY_DIR}/multiqc_quast_report.tsv"
CHECKM2_FILE="${ROOT_DIR}/${CHECKM2_SUMMARY_DIR}/checkm2_summary.tsv"
TAXONOMY_FILE="${ROOT_DIR}/${TAXONOMY_MERGED_DIR}/final_species_table.tsv"
MLST_FILE="${ROOT_DIR}/${MLST_SUMMARY_DIR}/complete_mlst_table.tsv"
BAKTA_FILE="${ROOT_DIR}/${ANNOTATION_SUMMARY_DIR}/multiqc_bakta_report.tsv"
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
    "${READS_QC_FILE}" \
    "${QUAST_FILE}" \
    "${CHECKM2_FILE}" \
    "${TAXONOMY_FILE}" \
    "${MLST_FILE}" \
    "${AMR_FILE}" \
    "${PLASMIDFINDER_FILE}" \
    "${MOBTYPER_FILE}"
do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required input file not found: ${f}"
        exit 1
    fi
done

if [[ ! -f "${BAKTA_FILE}" ]]; then
    echo "ERROR: Bakta file not found: ${BAKTA_FILE}"
    exit 1
fi

############################################
# RUN MERGE
############################################

echo "=================================================="
echo "STEP 07: FINAL SUMMARY"
echo "Samples file         : ${SAMPLES}"
echo "Reads QC file        : ${READS_QC_FILE}"
echo "QUAST file           : ${QUAST_FILE}"
echo "CheckM2 file         : ${CHECKM2_FILE}"
echo "Taxonomy file        : ${TAXONOMY_FILE}"
echo "MLST file            : ${MLST_FILE}"
echo "Bakta file            : ${BAKTA_FILE}"
echo "AMR file             : ${AMR_FILE}"
echo "PlasmidFinder file   : ${PLASMIDFINDER_FILE}"
echo "MOB-typer file       : ${MOBTYPER_FILE}"
echo "Output TSV           : ${OUT_TSV}"
echo "Output XLSX          : ${OUT_XLSX}"
echo "ENV_ANALYSIS         : ${ENV_ANALYSIS}"
echo "=================================================="

conda run --no-capture-output -n "${ENV_ANALYSIS}" python3 \
    "${SCRIPT_MERGE_FINAL_REPORT}" \
    --samples "${SAMPLES}" \
    --reads-qc "${READS_QC_FILE}" \
    --quast "${QUAST_FILE}" \
    --checkm2 "${CHECKM2_FILE}" \
    --taxonomy "${TAXONOMY_FILE}" \
    --mlst "${MLST_FILE}" \
    --bakta "${BAKTA_FILE}" \
    --amr "${AMR_FILE}" \
    --plasmidfinder "${PLASMIDFINDER_FILE}" \
    --mobtyper "${MOBTYPER_FILE}" \
    --out "${OUT_TSV}" \
    --out-xlsx "${OUT_XLSX}" \
    > "${FINAL_SUMMARY_LOG_DIR}/merge_final_report.stdout.log" \
    2> "${FINAL_SUMMARY_LOG_DIR}/merge_final_report.stderr.log"

echo "[DONE] Final summary written to: ${FINAL_SUMMARY_DIR}"