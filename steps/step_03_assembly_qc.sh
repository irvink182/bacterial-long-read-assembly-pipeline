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
log_info "Starting assembly_qc step"
[[ "${RESUME}" == "true" ]] && log_info "Resume mode enabled"

############################################
# DIRECTORIES
############################################

INPUT_FASTA_DIR="${ROOT_DIR}/${FINAL_ASM_DIR}"

QUAST_DIR="${ROOT_DIR}/${QUAST_DIR}"
QUAST_INDIVIDUAL_DIR="${ROOT_DIR}/${QUAST_INDIVIDUAL_DIR}"
QUAST_SUMMARY_DIR="${ROOT_DIR}/${QUAST_SUMMARY_DIR}"

CHECKM2_DIR="${ROOT_DIR}/${CHECKM2_DIR}"
CHECKM2_INDIVIDUAL_DIR="${ROOT_DIR}/${CHECKM2_INDIVIDUAL_DIR}"
CHECKM2_SUMMARY_DIR="${ROOT_DIR}/${CHECKM2_SUMMARY_DIR}"

QC_LOG_DIR="${ROOT_DIR}/${LOG_DIR}/assembly_qc"

mkdir -p "${QUAST_DIR}"
mkdir -p "${QUAST_INDIVIDUAL_DIR}"
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

while IFS=$'\t' read -r SAMPLE_ID ASM_TYPE EXPECTED_GENOME_SIZE LONG_READS SHORT_R1 SHORT_R2
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

    SAMPLE_QUAST_DIR="${QUAST_INDIVIDUAL_DIR}/${SAMPLE_ID}"
    SAMPLE_CHECKM2_DIR="${CHECKM2_INDIVIDUAL_DIR}/${SAMPLE_ID}"

    SAMPLE_LOG="${QC_LOG_DIR}/${SAMPLE_ID}.assembly_qc.log"

    # =========================
    # Define expected outputs
    # =========================
    QUAST_REPORT="${SAMPLE_QUAST_DIR}/report.tsv"
    CHECKM2_REPORT="${SAMPLE_CHECKM2_DIR}/${SAMPLE_ID}.${CHECKM2_PREFIX}.quality_report.tsv"

    RUN_QUAST=true
    RUN_CHECKM2=true

    # =========================
    # Resume logic (granular)
    # =========================
    if [[ "${RESUME}" == "true" && -s "${QUAST_REPORT}" ]]; then
        log_info "Skipping QUAST for ${SAMPLE_ID}"
        RUN_QUAST=false
    fi

    if [[ "${RESUME}" == "true" && -s "${CHECKM2_REPORT}" ]]; then
        log_info "Skipping CheckM2 for ${SAMPLE_ID}"
        RUN_CHECKM2=false
    fi

    # Si ambos existen → skip total
    if [[ "${RUN_QUAST}" == "false" && "${RUN_CHECKM2}" == "false" ]]; then
        log_info "Skipping ${SAMPLE_ID}: all QC outputs already exist"
        continue
    fi

    mkdir -p "${SAMPLE_QUAST_DIR}"
    mkdir -p "${SAMPLE_CHECKM2_DIR}"

    {
        echo "[INFO] Sample: ${SAMPLE_ID}"
        echo "[INFO] Assembly: ${ASSEMBLY_FASTA}"

        # =========================
        # QUAST
        # =========================
        if [[ "${RUN_QUAST}" == "true" ]]; then
            echo "[INFO] QUAST: ${SAMPLE_ID}"
            run_quast "${SAMPLE_ID}" "${ASSEMBLY_FASTA}" "${SAMPLE_QUAST_DIR}"
            echo "[INFO] QUAST: ${SAMPLE_ID} finished"
        else
            echo "[INFO] QUAST skipped"
        fi

        # =========================
        # CHECKM2
        # =========================
        if [[ "${RUN_CHECKM2}" == "true" ]]; then
            echo "[INFO] CHECKM2: ${SAMPLE_ID}"
            run_checkm2 "${SAMPLE_ID}" "${ASSEMBLY_FASTA}" "${SAMPLE_CHECKM2_DIR}"
            echo "[INFO] CHECKM2: ${SAMPLE_ID} finished"
        else
            echo "[INFO] CHECKM2 skipped"
        fi

    } > "${SAMPLE_LOG}" 2>&1

    echo "Done: ${SAMPLE_ID}"

done < <(tail -n +2 "${SAMPLES}")

# Count number of samples (excluding header)
N_SAMPLES=$(($(wc -l < "${SAMPLES}") - 1))
log_info "Detected ${N_SAMPLES} samples"

############################################
# MULTIQC FOR QUAST
############################################

############################################
# MULTIQC FOR QUAST
############################################

echo
echo "[INFO] Creating QUAST summary report"

QUAST_SUMMARY_REPORT="${QUAST_SUMMARY_DIR}/multiqc_quast_report.tsv"

# Skip if resume and already exists
if should_skip_global "${QUAST_SUMMARY_REPORT}"; then
    log_info "Skipping MultiQC (QUAST)"

# Skip if only 1 sample
elif [[ "${N_SAMPLES}" -lt 2 ]]; then
    log_warn "Only one sample detected → skipping MultiQC (QUAST)"

else
    log_info "Running MultiQC (QUAST)"

    cd "${QUAST_INDIVIDUAL_DIR}"

    conda run -n "${ENV_QUAST}" multiqc . \
        --force \
        --cl-config "max_table_rows: 3000" \
        > "${QC_LOG_DIR}/multiqc_quast.stdout.log" \
        2> "${QC_LOG_DIR}/multiqc_quast.stderr.log" || true

    # Handle outputs safely
    if [[ -f multiqc_report.html ]]; then
        mv multiqc_report.html multiqc_quast_report.html
        cp multiqc_quast_report.html "${QUAST_SUMMARY_DIR}/"
    else
        log_warn "MultiQC HTML report not generated"
    fi

    if [[ -f multiqc_data/multiqc_quast.txt ]]; then
        mv multiqc_data/multiqc_quast.txt multiqc_data/multiqc_quast_report.tsv
        sed -i 's/\.0//g' multiqc_data/multiqc_quast_report.tsv
        cp multiqc_data/multiqc_quast_report.tsv "${QUAST_SUMMARY_DIR}/"
    else
        log_warn "MultiQC QUAST table not generated"
    fi

    cd "${ROOT_DIR}"
fi

############################################
# CHECKM2 SUMMARY
############################################

echo "[INFO] Creating CheckM2 summary report"

CHECKM2_SUMMARY_FILE="${CHECKM2_SUMMARY_DIR}/checkm2_summary.tsv"
FIRST=1
> "${CHECKM2_SUMMARY_FILE}"

while IFS=$'\t' read -r SAMPLE_ID ASM_TYPE EXPECTED_GENOME_SIZE LONG_READS SHORT_R1 SHORT_R2
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
done < <(tail -n +2 "${SAMPLES}")

echo "[INFO] ASSEMBLY QC module completed"