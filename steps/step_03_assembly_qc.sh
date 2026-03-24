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

INPUT_FASTA_DIR="${ROOT_DIR}/${FINAL_ASM_DIR}"

QUAST_DIR="${ROOT_DIR}/${QUAST_DIR}"
QUAST_SUMMARY_DIR="${ROOT_DIR}/${QUAST_SUMMARY_DIR}"

CHECKM2_DIR="${ROOT_DIR}/${CHECKM2_DIR}"
CHECKM2_INDIVIDUAL_DIR="${ROOT_DIR}/${CHECKM2_INDIVIDUAL_DIR}"
CHECKM2_SUMMARY_DIR="${ROOT_DIR}/${CHECKM2_SUMMARY_DIR}"

QC_LOG_DIR="${ROOT_DIR}/${LOG_DIR}/assembly_qc"

mkdir -p "${QUAST_DIR}"
mkdir -p "${QUAST_SUMMARY_DIR}"
mkdir -p "${CHECKM2_DIR}"
mkdir -p "${CHECKM2_INDIVIDUAL_DIR}"
mkdir -p "${CHECKM2_SUMMARY_DIR}"
mkdir -p "${QC_LOG_DIR}"

############################################
# CHECK CONDA ENV
############################################

if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_QUAST}"; then
    echo "ERROR: conda environment not found: ${ENV_QUAST}"
    exit 1
fi

if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_CHECKM2}"; then
    echo "ERROR: conda environment not found: ${ENV_CHECKM2}"
    exit 1
fi

############################################
# FUNCTIONS
############################################

run_quast() {
    local sample_id="$1"
    local assembly_fasta="$2"
    local outdir="$3"

    conda run -n "${ENV_QUAST}" quast.py \
        --glimmer \
        --threads "${QUAST_THREADS}" \
        --min-contig "${QUAST_MIN_CONTIG}" \
        -l "${sample_id}" \
        "${assembly_fasta}" \
        -o "${outdir}"
}

run_checkm2() {
    local sample_id="$1"
    local assembly_fasta="$2"
    local outdir="$3"

    conda run -n "${ENV_CHECKM2}" checkm2 predict \
        --database_path "${CHECKM2_DB}" \
        --allmodels \
        --threads "${CHECKM2_THREADS}" \
        --input "${assembly_fasta}" \
        --output-directory "${outdir}" \
        --force

    if [[ -f "${outdir}/quality_report.tsv" ]]; then
        mv "${outdir}/quality_report.tsv" \
           "${outdir}/${sample_id}.${CHECKM2_PREFIX}.quality_report.tsv"
    fi
}

############################################
# MAIN
############################################

echo "=================================================="
echo "STEP 03: ASSEMBLY QC"
echo "Samples file       : ${SAMPLES}"
echo "Input fasta dir    : ${INPUT_FASTA_DIR}"
echo "QUAST dir          : ${QUAST_DIR}"
echo "CheckM2 dir        : ${CHECKM2_DIR}"
echo "Conda env          : ${ENV_QUAST} & ${ENV_CHECKM2}"
echo "=================================================="

tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SAMPLE_ID ASM_TYPE EXPECTED_GENOME_SIZE LONG_READS SHORT_R1 SHORT_R2
do
    [[ -z "${SAMPLE_ID}" ]] && continue

    echo
    echo "--------------------------------------------------"
    echo "Processing sample: ${SAMPLE_ID}"
    echo "--------------------------------------------------"

    validate_asm_type "${ASM_TYPE}"


    if is_na "${EXPECTED_GENOME_SIZE}"; then
        log_warn "sample ${SAMPLE_ID} has empty expected_genome_size"
    else
        validate_genome_size "${EXPECTED_GENOME_SIZE}"
    fi

    ASSEMBLY_FASTA="${INPUT_FASTA_DIR}/${SAMPLE_ID}.${FINAL_ASM_PREFIX}.fasta"
    require_file "${ASSEMBLY_FASTA}" "assembly fasta for ${SAMPLE_ID}"

    SAMPLE_QUAST_DIR="${QUAST_DIR}/${SAMPLE_ID}"
    SAMPLE_CHECKM2_DIR="${CHECKM2_INDIVIDUAL_DIR}/${SAMPLE_ID}"

    SAMPLE_LOG="${QC_LOG_DIR}/${SAMPLE_ID}.assembly_qc.log"

    mkdir -p "${SAMPLE_QUAST_DIR}"
    mkdir -p "${SAMPLE_CHECKM2_DIR}"

    if [[ ! -f "${ASSEMBLY_FASTA}" ]]; then
        echo "ERROR: final assembly not found for ${SAMPLE_ID}: ${ASSEMBLY_FASTA}"
        exit 1
    fi

    {
        echo "[INFO] Sample: ${SAMPLE_ID}"
        echo "[INFO] Assembly: ${ASSEMBLY_FASTA}"

        echo "[INFO] QUAST: ${SAMPLE_ID}"
        run_quast "${SAMPLE_ID}" "${ASSEMBLY_FASTA}" "${SAMPLE_QUAST_DIR}"
        echo "[INFO] QUAST: ${SAMPLE_ID} finished"

        echo "[INFO] CHECKM2: ${SAMPLE_ID}"
        run_checkm2 "${SAMPLE_ID}" "${ASSEMBLY_FASTA}" "${SAMPLE_CHECKM2_DIR}"
        echo "[INFO] CHECKM2: ${SAMPLE_ID} finished"

    } > "${SAMPLE_LOG}" 2>&1

    echo "Done: ${SAMPLE_ID}"
done

############################################
# MULTIQC FOR QUAST
############################################

echo
echo "[INFO] Creating QUAST summary report"

cd "${QUAST_DIR}"

conda run -n "${ENV_QUAST}" multiqc . --force --cl-config "max_table_rows: 3000" \
    > "${QC_LOG_DIR}/multiqc_quast.stdout.log" \
    2> "${QC_LOG_DIR}/multiqc_quast.stderr.log" || true

if [[ -f "multiqc_report.html" ]]; then
    mv multiqc_report.html multiqc_quast_report.html
fi

if [[ -f "multiqc_data/multiqc_quast.txt" ]]; then
    mv multiqc_data/multiqc_quast.txt multiqc_data/multiqc_quast_report.tsv
    sed -i 's/\.0//g' multiqc_data/multiqc_quast_report.tsv
    cp multiqc_data/multiqc_quast_report.tsv "${QUAST_SUMMARY_DIR}/"
fi

if [[ -f "multiqc_quast_report.html" ]]; then
    cp multiqc_quast_report.html "${QUAST_SUMMARY_DIR}/"
fi

cd "${ROOT_DIR}"

############################################
# CHECKM2 SUMMARY
############################################

echo "[INFO] Creating CheckM2 summary report"

CHECKM2_SUMMARY_FILE="${CHECKM2_SUMMARY_DIR}/checkm2_summary.tsv"
FIRST=1
> "${CHECKM2_SUMMARY_FILE}"

tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SAMPLE_ID ASM_TYPE EXPECTED_GENOME_SIZE LONG_READS SHORT_R1 SHORT_R2
do
    [[ -z "${SAMPLE_ID}" ]] && continue

    REPORT="${CHECKM2_INDIVIDUAL_DIR}/${SAMPLE_ID}/${SAMPLE_ID}.${CHECKM2_PREFIX}.quality_report.tsv"

    if [[ -f "${REPORT}" ]]; then
        if [[ ${FIRST} -eq 1 ]]; then
            awk -v sample="${SAMPLE_ID}" 'BEGIN{FS=OFS="\t"} NR==1{print "Sample",$0} NR>1{print sample,$0}' "${REPORT}" \
                > "${CHECKM2_SUMMARY_FILE}"
            FIRST=0
        else
            awk -v sample="${SAMPLE_ID}" 'BEGIN{FS=OFS="\t"} NR>1{print sample,$0}' "${REPORT}" \
                >> "${CHECKM2_SUMMARY_FILE}"
        fi
    fi
done

echo "[INFO] ASSEMBLY QC module completed"