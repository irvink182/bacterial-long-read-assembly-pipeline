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

READS_DIR="${ROOT_DIR}/${CLEAN_READS_DIR}"
RAW_DIR="${ROOT_DIR}/${RAW_READS_DIR}"

FLYE_DIR="${ROOT_DIR}/${FLYE_DIR}"
MEDAKA_DIR="${ROOT_DIR}/${MEDAKA_DIR}"
FINAL_ASM_DIR="${ROOT_DIR}/${FINAL_ASM_DIR}"

ASSEMBLY_LOG_DIR="${ROOT_DIR}/${LOG_DIR}/assembly"

mkdir -p "${FLYE_DIR}"
mkdir -p "${MEDAKA_DIR}"
mkdir -p "${FINAL_ASM_DIR}"
mkdir -p "${ASSEMBLY_LOG_DIR}"

############################################
# CHECK CONDA ENV
############################################

if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_ASSEMBLY}"; then
    echo "ERROR: conda environment not found: ${ENV_ASSEMBLY}"
    exit 1
fi

############################################
# FUNCTIONS
############################################

run_flye() {
    local reads_fastq="$1"
    local genome_size="$2"
    local outdir="$3"

    conda run -n "${ENV_ASSEMBLY}" flye \
        --nano-hq "${reads_fastq}" \
        --threads "${FLYE_THREADS}" \
        --asm-coverage "${FLYE_ASM_COVERAGE}" \
        --genome-size "${genome_size}" \
        --iterations "${FLYE_ITERATIONS}" \
        --out-dir "${outdir}"
}

run_medaka() {
    local reads_fastq="$1"
    local draft_fasta="$2"
    local outdir="$3"

    conda run -n "${ENV_ASSEMBLY}" medaka_consensus \
        -b "${MEDAKA_BATCH_SIZE}" \
        -t "${MEDAKA_THREADS}" \
        -m "${MEDAKA_MODEL}" \
        -i "${reads_fastq}" \
        -d "${draft_fasta}" \
        -o "${outdir}"
}

sort_and_rename_contigs() {
    local input_fasta="$1"
    local output_fasta="$2"
    local sample_id="$3"
    local tmp_sorted="$4"

    conda run -n "${ENV_ASSEMBLY}" seqkit sort \
        -l -r "${input_fasta}" \
        > "${tmp_sorted}"

    conda run -n "${ENV_ASSEMBLY}" seqkit replace \
        -p ".*" \
        -r "${sample_id}_contig_{nr}" \
        "${tmp_sorted}" \
        > "${output_fasta}"
}

run_pypolca() {
    local assembly_fasta="$1"
    local short_r1="$2"
    local short_r2="$3"
    local sample_id="$4"
    local output_dir="$5"

    mkdir -p "${output_dir}"

    conda run -n "${ENV_ASSEMBLY}" pypolca run \
        --assembly "${assembly_fasta}" \
        --reads1 "${short_r1}" \
        --reads2 "${short_r2}" \
        --prefix "${sample_id}" \
        --threads "${PYPOLCA_THREADS}" \
        --careful \
        --force \
        --output "${output_dir}"
}

############################################
# MAIN
############################################

echo "=================================================="
echo "STEP 02: ASSEMBLY"
echo "Samples file     : ${SAMPLES}"
echo "Clean reads dir  : ${READS_DIR}"
echo "Flye dir         : ${FLYE_DIR}"
echo "Medaka dir       : ${MEDAKA_DIR}"
echo "Final asm dir    : ${FINAL_ASM_DIR}"
echo "Conda env        : ${ENV_ASSEMBLY}"
echo "TRIMMER          : ${TRIMMER}"
echo "=================================================="

tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SAMPLE_ID ASM_TYPE EXPECTED_GENOME_SIZE LONG_READS SHORT_R1 SHORT_R2
do
    [[ -z "${SAMPLE_ID}" ]] && continue

    echo
    echo "--------------------------------------------------"
    echo "Processing sample: ${SAMPLE_ID}"
    echo "asm_type             : ${ASM_TYPE}"
    echo "expected_genome_size : ${EXPECTED_GENOME_SIZE}"
    echo "--------------------------------------------------"

    validate_asm_type "${ASM_TYPE}"

    if is_na "${EXPECTED_GENOME_SIZE}"; then
        log_warn "sample ${SAMPLE_ID} has empty expected_genome_size"
    else
        validate_genome_size "${EXPECTED_GENOME_SIZE}"
    fi

    if [[ "${ASM_TYPE}" != "long" && "${ASM_TYPE}" != "hybrid" ]]; then
        echo "WARNING: sample ${SAMPLE_ID} has asm_type=${ASM_TYPE}. Skipping in step_02_assembly.sh"
        continue
    fi

    if [[ -z "${EXPECTED_GENOME_SIZE}" ]]; then
        echo "WARNING: sample ${SAMPLE_ID} has empty expected_genome_size. Using default: ${FLYE_GENOME_SIZE_DEFAULT}"
        EXPECTED_GENOME_SIZE="${FLYE_GENOME_SIZE_DEFAULT}"
    else
        validate_genome_size "${EXPECTED_GENOME_SIZE}"
    fi

    if [[ "${ASM_TYPE}" == "hybrid" ]]; then
        SHORT_R1_PATH="${RAW_DIR}/${SHORT_R1}"
        SHORT_R2_PATH="${RAW_DIR}/${SHORT_R2}"

        require_file "${SHORT_R1_PATH}" "short_r1 fo ${SAMPLE_ID}"
        require_file "${SHORT_R2_PATH}" "short_r2 fo ${SAMPLE_ID}"
        log_info "Sample ${SAMPLE_ID} configured for hybrid assembly"
    else
        log_info "Sample ${SAMPLE_ID} configured for only long-read assembly"
    fi


    CLEAN_FASTQ="${READS_DIR}/${SAMPLE_ID}.${TRIMMER}.${FILTLONG_PREFIX}.fastq.gz"
    require_file "${CLEAN_FASTQ}" "clean FASTQ for ${SAMPLE_ID}"

    SAMPLE_FLYE_DIR="${FLYE_DIR}/${SAMPLE_ID}"
    SAMPLE_MEDAKA_DIR="${MEDAKA_DIR}/${SAMPLE_ID}"

    SAMPLE_LOG="${ASSEMBLY_LOG_DIR}/${SAMPLE_ID}.assembly.log"

    mkdir -p "${SAMPLE_FLYE_DIR}"
    mkdir -p "${SAMPLE_MEDAKA_DIR}"

    require_file "${CLEAN_FASTQ}" clean FASTQ for ${SAMPLE_ID}

    {
        echo "[INFO] Sample: ${SAMPLE_ID}"
        echo "[INFO] Clean FASTQ: ${CLEAN_FASTQ}"
        echo "[INFO] Genome size: ${EXPECTED_GENOME_SIZE}"

        ############################################
        # FLYE
        ############################################
        echo "[INFO] Running Flye"
        run_flye "${CLEAN_FASTQ}" "${EXPECTED_GENOME_SIZE}" "${SAMPLE_FLYE_DIR}"

        require_file "${SAMPLE_FLYE_DIR}/assembly.fasta" "Flye assembly.fasta for ${SAMPLE_ID}"

        cd "${SAMPLE_FLYE_DIR}"

        mv assembly.fasta "${SAMPLE_ID}.${FLYE_PREFIX}.fasta"
        mv assembly_graph.gfa "${SAMPLE_ID}.${FLYE_PREFIX}.gfa"
        mv flye.log "${SAMPLE_ID}.${FLYE_PREFIX}.log"
        mv assembly_info.txt "${SAMPLE_ID}.${FLYE_PREFIX}.info.txt"

        FLYE_FASTA="${SAMPLE_FLYE_DIR}/${SAMPLE_ID}.${FLYE_PREFIX}.fasta"

        echo "[INFO] Flye finished"

        ############################################
        # MEDAKA
        ############################################
        echo "[INFO] Running Medaka"
        run_medaka "${CLEAN_FASTQ}" "${FLYE_FASTA}" "${SAMPLE_MEDAKA_DIR}"

        require_file "${SAMPLE_MEDAKA_DIR}/consensus.fasta" "Medaka consensus.fasta for ${SAMPLE_ID}"

        echo "[INFO] Sorting and renaming contigs"
        TMP_SORTED="${SAMPLE_MEDAKA_DIR}/consensus.sorted.fasta"
        RENAMED_FASTA="${SAMPLE_MEDAKA_DIR}/${SAMPLE_ID}.${MEDAKA_PREFIX}.fasta"

        sort_and_rename_contigs \
            "${SAMPLE_MEDAKA_DIR}/consensus.fasta" \
            "${RENAMED_FASTA}" \
            "${SAMPLE_ID}" \
            "${TMP_SORTED}"

        rm -f "${TMP_SORTED}"
        rm -f "${SAMPLE_MEDAKA_DIR}/consensus.fasta"

        require_file "${RENAMED_FASTA}" "renamed Medaka fasta for ${SAMPLE_ID}"

        FINAL_INPUT_FASTA="${RENAMED_FASTA}"

        ############################################
        # HYBRID ASSEMBLY
        ############################################

        if [[ "${ASM_TYPE}" == "hybrid" ]]; then
            SHORT_R1_PATH="${RAW_DIR}/${SHORT_R1}"
            SHORT_R2_PATH="${RAW_DIR}/${SHORT_R2}"
            require_file "${SHORT_R1_PATH}" "short_r1 for ${SAMPLE_ID}"
            require_file "${SHORT_R2_PATH}" "short_r2 for ${SAMPLE_ID}"

            SAMPLE_PYPOLCA_DIR="${SAMPLE_MEDAKA_DIR}/pypolca"

            echo "[INFO] Running pypolca"
            run_pypolca "${RENAMED_FASTA}" "${SHORT_R1_PATH}" "${SHORT_R2_PATH}" "${SAMPLE_ID}" "${SAMPLE_PYPOLCA_DIR}"

            PYPOLCA_FASTA="${SAMPLE_PYPOLCA_DIR}/${SAMPLE_ID}_corrected.fasta"
            require_file "${PYPOLCA_FASTA}" "pypolca polished fasta for ${SAMPLE_ID}"

            FINAL_INPUT_FASTA="${PYPOLCA_FASTA}"
        fi

        ############################################
        # FINAL ASSEMBLY
        ############################################
        cp "${FINAL_INPUT_FASTA}" "${FINAL_ASM_DIR}/${SAMPLE_ID}.${FINAL_ASM_PREFIX}.fasta"

        require_file "${FINAL_ASM_DIR}/${SAMPLE_ID}.${FINAL_ASM_PREFIX}.fasta" "final assembly for ${SAMPLE_ID}"

        echo "[INFO] Final assembly written to: ${FINAL_ASM_DIR}/${SAMPLE_ID}.${FINAL_ASM_PREFIX}.fasta"
        echo "[INFO] Sample ${SAMPLE_ID} completed successfully"
    } > "${SAMPLE_LOG}" 2>&1

    echo "Done: ${SAMPLE_ID}"
done

echo
echo "[INFO] ASSEMBLY module completed"