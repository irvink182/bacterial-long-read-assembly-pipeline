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
log_info "Starting annotation step"
[[ "${RESUME}" == "true" ]] && log_info "Resume mode enabled"

############################################
# DIRECTORIES
############################################

INPUT_FASTA_DIR="${ROOT_DIR}/${FINAL_ASM_DIR}"

BAKTA_DIR="${ROOT_DIR}/${BAKTA_DIR}"
ANNOTATION_SUMMARY_DIR="${ROOT_DIR}/${ANNOTATION_SUMMARY_DIR}"

ANNOTATION_LOG_DIR="${ROOT_DIR}/${LOG_DIR}/annotation"

mkdir -p "${BAKTA_DIR}"
mkdir -p "${ANNOTATION_SUMMARY_DIR}"
mkdir -p "${ANNOTATION_LOG_DIR}"

############################################
# STATUS FILE
############################################

STATUS_FILE="${ROOT_DIR}/results/pipeline_status.tsv"
init_status_file

export STATUS_FILE
export UTILS_FILE

############################################
# CHECK CONDA ENV
############################################

if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_ANNOTATION}"; then
    echo "ERROR: conda environment not found: ${ENV_ANNOTATION}"
    exit 1
fi

############################################
# FUNCTIONS
############################################

run_bakta() {
    local sample_id="$1"
    local assembly_fasta="$2"
    local outdir="$3"

    conda run --no-capture-output -n "${ENV_ANNOTATION}" \
        bakta \
        --db "${BAKTA_DB}" \
        --force \
        --compliant \
        --threads "${BAKTA_THREADS}" \
        --prefix "${sample_id}" \
        --locus-tag "${sample_id}" \
        --keep-contig-headers \
        "${assembly_fasta}" \
        --output "${outdir}"
}

make_nosequence_gbff() {
    local sample_id="$1"
    local outdir="$2"

    local gbff="${outdir}/${sample_id}.gbff"
    local gbff_noseq="${outdir}/${sample_id}_nosequence.gbff"

    if [[ -f "${gbff}" ]]; then
        sed '/^##FASTA/Q' "${gbff}" > "${gbff_noseq}"
    fi
}

############################################
# MAIN
############################################

echo "=================================================="
echo "STEP 05: ANNOTATION"
echo "Samples file       : ${SAMPLES}"
echo "Input fasta dir    : ${INPUT_FASTA_DIR}"
echo "Bakta dir          : ${BAKTA_DIR}"
echo "Summary dir        : ${ANNOTATION_SUMMARY_DIR}"
echo "ENV_ANNOTATION     : ${ENV_ANNOTATION}"
echo "=================================================="

while IFS=$'\t' read -r SAMPLE_ID ASM_TYPE EXPECTED_GENOME_SIZE LONG_READS SHORT_R1 SHORT_R2
do
    [[ -z "${SAMPLE_ID}" ]] && continue

    echo
    echo "--------------------------------------------------"
    echo "Processing sample: ${SAMPLE_ID}"
    echo "--------------------------------------------------"

    # --- FILTERS FOR EMPTY/FAILED SAMPLES ---
    # # Checking if the sample failed in previous steps
    if grep -qE "${SAMPLE_ID}.*FAILED" "${STATUS_FILE}"; then
        log_warn "Skipping ${SAMPLE_ID}: Sample has FAILED status in previous steps."
        continue
    fi

    # Checking if the sample has assembly, if not, SKIPPED.
    ASSEMBLY_FASTA="${INPUT_FASTA_DIR}/${SAMPLE_ID}.${FINAL_ASM_PREFIX}.fasta"
    
    if [[ ! -f "${ASSEMBLY_FASTA}" || ! -s "${ASSEMBLY_FASTA}" ]]; then
        log_warn "Assembly not found for ${SAMPLE_ID}. Marking as SKIPPED in QC."
        status_skip "$SAMPLE_ID" "annotation" "missing_assembly"
        continue
    fi

    #Status per file
    status_start "$SAMPLE_ID" "annotation"
    
    validate_asm_type "${ASM_TYPE}"

    if is_na "${EXPECTED_GENOME_SIZE}"; then
        log_warn "sample ${SAMPLE_ID} has empty expected_genome_size"
    else
        validate_genome_size "${EXPECTED_GENOME_SIZE}"
    fi

    ASSEMBLY_FASTA="${INPUT_FASTA_DIR}/${SAMPLE_ID}.${FINAL_ASM_PREFIX}.fasta"
    
    SAMPLE_BAKTA_DIR="${BAKTA_DIR}/${SAMPLE_ID}"
    BAKTA_FILE="${SAMPLE_BAKTA_DIR}/${SAMPLE_ID}.gff3"
    SAMPLE_LOG="${ANNOTATION_LOG_DIR}/${SAMPLE_ID}.bakta.log"

    mkdir -p "${SAMPLE_BAKTA_DIR}"

    # =========================
    # Checking inputs
    # =========================

    #Checking input assembly
    if ! check_file_not_empty "${ASSEMBLY_FASTA}"; then
        log_warn "Assembly not found: ${SAMPLE_ID}"
        status_skip "$SAMPLE_ID" "annotation" "missing_assembly"
        continue
    fi

    # RESUME LOGIC
    if should_skip_sample "${BAKTA_FILE}"; then
        status_skip "$SAMPLE_ID" "annotation" "resume_existing_output"
        continue
    fi

    echo "[INFO] Sample: ${SAMPLE_ID}"
    echo "[INFO] Assembly: ${ASSEMBLY_FASTA}"
    
    # =========================
    # BAKTA
    # =========================
    #RUN BAKTA BLOCK

    {
        echo "[INFO] Running Bakta: ${SAMPLE_ID}"
        run_safe "bakta" run_bakta "${SAMPLE_ID}" "${ASSEMBLY_FASTA}" "${SAMPLE_BAKTA_DIR}"

        echo "[INFO] Creating no-sequence GBFF"
        run_safe "gbff" make_nosequence_gbff "${SAMPLE_ID}" "${SAMPLE_BAKTA_DIR}"

        echo "[INFO] Bakta annotation finished: ${SAMPLE_ID}"
    } >> "${SAMPLE_LOG}" 2>&1
    
    # Checking bakta output
    if ! check_file_not_empty "${BAKTA_FILE}"; then
        log_warn "Bakta output file empty for ${SAMPLE_ID}"
        status_fail "$SAMPLE_ID" "annotation" "bakta_failed"
        continue
    fi

    #FINAL LOG STATUS
    status_ok "$SAMPLE_ID" "annotation"

    echo "Done: ${SAMPLE_ID}"

done < <(tail -n +2 "${SAMPLES}")

############################################
# MULTIQC REPORT
############################################

echo
echo "[INFO] Creating Bakta MultiQC report"

BAKTA_SUMMARY_REPORT="${ANNOTATION_SUMMARY_DIR}/multiqc_bakta_report.tsv"
# RESUME LOGIC
if should_skip_global "${BAKTA_SUMMARY_REPORT}"; then
    log_info "Skipping MultiQC (BAKTA)"
else
    log_info "Running MultiQC (BAKTA)"
    cd "${BAKTA_DIR}"
    conda run --no-capture-output -n "${ENV_ANNOTATION}" \
    multiqc . --force --cl-config "max_table_rows: 3000" \
    > "${ANNOTATION_LOG_DIR}/multiqc_bakta.stdout.log" \
    2> "${ANNOTATION_LOG_DIR}/multiqc_bakta.stderr.log" || true

    if [[ -f multiqc_report.html ]]; then
        mv multiqc_report.html multiqc_bakta_report.html
        cp multiqc_bakta_report.html "${ANNOTATION_SUMMARY_DIR}/"
    else
        log_warn "MultiQC HTML report not generated"
    fi

    if [[ -f multiqc_data/multiqc_bakta.txt ]]; then
        mv multiqc_data/multiqc_bakta.txt multiqc_data/multiqc_bakta_report.tsv
        cp multiqc_data/multiqc_bakta_report.tsv "${ANNOTATION_SUMMARY_DIR}/"
    else
        log_warn "MultiQC table not generated"
    fi

    cd "${ROOT_DIR}"
fi

#Final log info
log_info "Annotation step completed for all samples"
awk '$2=="annotation"' "${STATUS_FILE}"
echo "[INFO] ANNOTATION module completed"
echo "[DONE] Outputs written to: ${ROOT_DIR}/${ANNOTATION_RESULTS_DIR}"