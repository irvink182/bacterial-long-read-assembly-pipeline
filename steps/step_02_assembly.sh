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

#Helper scripts
SCRIPT_CONTIGS_STATS="${ROOT_DIR}/scripts/generate_contigs_stats.py"

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
log_info "Starting assembly step"
[[ "${RESUME}" == "true" ]] && log_info "Resume mode enabled"

############################################
# DIRECTORIES
############################################

READS_DIR="${ROOT_DIR}/${CLEAN_READS_DIR}"
RAW_DIR="${ROOT_DIR}/${RAW_READS_DIR}"

FLYE_DIR="${ROOT_DIR}/${FLYE_DIR}"
MEDAKA_DIR="${ROOT_DIR}/${MEDAKA_DIR}"
DNAAPLER_DIR="${ROOT_DIR}/${DNAAPLER_DIR}"
COVERAGE_DIR="${ROOT_DIR}/${COVERAGE_DIR}"
FINAL_ASM_DIR="${ROOT_DIR}/${FINAL_ASM_DIR}"

ASSEMBLY_LOG_DIR="${ROOT_DIR}/${LOG_DIR}/assembly"

mkdir -p "${FLYE_DIR}"
mkdir -p "${MEDAKA_DIR}"
mkdir -p "${COVERAGE_DIR}"
mkdir -p "${DNAAPLER_DIR}"
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
        -f \
        -i "${reads_fastq}" \
        -d "${draft_fasta}" \
        -o "${outdir}"
}

run_mosdepth() {
    local input_fasta="$1"
    local reads_fastq="$2"
    local tmp_bam="$3"
    local output_bam="$4"
    local cov_prefix="$5"
    
    conda run -n "${ENV_ASSEMBLY}" minimap2 -ax \
        map-ont -t "${MINIMAP_THREADS}" \
        "${input_fasta}" "${reads_fastq}" \
        > "${tmp_bam}"

    conda run -n "${ENV_ASSEMBLY}" samtools sort \
        "${tmp_bam}" -t "${SAMTOOLS_THREADS}" -o "${output_bam}"

    conda run -n "${ENV_ASSEMBLY}" samtools index \
        "${output_bam}"

    conda run -n "${ENV_ASSEMBLY}" mosdepth \
        --by 1000 \
        "${cov_prefix}" "${output_bam}"
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

run_dnaapler() {
    local input_fasta="$1"
    local output_dir="$2"

    conda run -n "${ENV_ASSEMBLY}" dnaapler "${DNAAPLER_MODE}" \
        --force \
        --threads "${DNAAPLER_THREADS}" \
        --input "${input_fasta}" \
        --output "${output_dir}"
}

clean_dnaapler_headers() {
    local input_fasta="$1"
    local output_fasta="$2"

    conda run -n "${ENV_ASSEMBLY}" seqkit replace \
        -p "^(\S+).*" \
        -r "\$1" \
        "${input_fasta}" \
        > "${output_fasta}"
}

extract_dnaapler_metadata() {
    local input_fasta="$1"
    local output_tsv="$2"

    grep "^>" "${input_fasta}" | sed 's/>//' | \
    awk '{
        id=$1
        rotated=$2
        gene=$3
        print id"\t"rotated"\t"gene
    }' > "${output_tsv}"
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
echo "Dnaapler dir     : ${DNAAPLER_DIR}"
echo "Coverage dir     : ${COVERAGE_DIR}"
echo "Final asm dir    : ${FINAL_ASM_DIR}"
echo "Conda env        : ${ENV_ASSEMBLY}"
echo "TRIMMER          : ${TRIMMER}"
echo "=================================================="

while IFS=$'\t' read -r SAMPLE_ID ASM_TYPE EXPECTED_GENOME_SIZE LONG_READS SHORT_R1 SHORT_R2
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

    FINAL_ASSEMBLY="${FINAL_ASM_DIR}/${SAMPLE_ID}.${FINAL_ASM_PREFIX}.fasta"

    COV_PREFIX="${COVERAGE_DIR}/${SAMPLE_ID}.coverage"
    TMP_BAM="$COVERAGE_DIR/${SAMPLE_ID}.tmp.bam"
    BAM="${COVERAGE_DIR}/${SAMPLE_ID}.bam"

    # RESUME LOGIC
    should_skip_sample "${FINAL_ASSEMBLY}" && continue

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

        ############################################
        # CONTIG COVERAGE
        ############################################
        echo "[INFO] Calculating contig coverage"

        if [[ "${RESUME}" == "true" && -f "${COV_PREFIX}.mosdepth.summary.txt" ]]; then
            echo "[SKIP] coverage ${SAMPLE_ID}"
        else
            run_mosdepth \
            "${SAMPLE_MEDAKA_DIR}/consensus.fasta" \
            "${CLEAN_FASTQ}" \
            "${TMP_BAM}" \
            "${BAM}" \
            "${COV_PREFIX}"

        fi

        rm -f "${TMP_BAM}"
        rm -f "${BAM}"
        rm -f "${BAM}.bai"

        ############################################
        # FILTER CONTIGS BY COVERAGE
        ############################################

        FILTERED_FASTA="${SAMPLE_MEDAKA_DIR}/consensus.filtered.fasta"

        if [[ "${FILTER_CONTIGS}" == "true" ]]; then
            echo "[INFO] Filtering contigs by coverage"

            SUMMARY="${COV_PREFIX}.mosdepth.summary.txt"

            require_file "${SUMMARY}" "mosdepth summary for ${SAMPLE_ID}"

            # Get chromosome coverage (largest contig)

            CHROM_COV=$(awk '
            NR > 1 && $1 !~ /total/ {
                print $4, $2
            }' "${SUMMARY}" | sort -k2,2nr | head -n1 | awk '{print $1}')

            echo "[INFO] Chromosome coverage: ${CHROM_COV}"

            MIN_COV=$(awk -v c="$CHROM_COV" -v f="$MIN_COV_FACTOR" 'BEGIN{print c*f}')

            echo "[INFO] Min coverage threshold: ${MIN_COV}"

            # Select contigs
            awk -v mincov="$MIN_COV" -v minlen="$MIN_CONTIG_LENGTH" '
            NR > 1 && $1 !~ /total/ {
                if ($4 >= mincov && $2 >= minlen) print $1
                }' "${SUMMARY}" > "${SAMPLE_MEDAKA_DIR}/contigs.keep.txt"

            # Filter FASTA
            conda run -n "${ENV_ASSEMBLY}" seqkit grep \
            -f "${SAMPLE_MEDAKA_DIR}/contigs.keep.txt" \
            "${SAMPLE_MEDAKA_DIR}/consensus.fasta" \
            > "${FILTERED_FASTA}"

            FINAL_INPUT_FASTA="${FILTERED_FASTA}"
        else
            FINAL_INPUT_FASTA="${SAMPLE_MEDAKA_DIR}/consensus.fasta"
        fi

        ############################################
        # SORTING AND RENAMING CONTIGS
        ############################################

        echo "[INFO] Sorting and renaming contigs"
        TMP_SORTED="${SAMPLE_MEDAKA_DIR}/consensus.sorted.fasta"
        RENAMED_FASTA="${SAMPLE_MEDAKA_DIR}/${SAMPLE_ID}.${MEDAKA_PREFIX}.fasta"

        sort_and_rename_contigs \
            "${FINAL_INPUT_FASTA}" \
            "${RENAMED_FASTA}" \
            "${SAMPLE_ID}" \
            "${TMP_SORTED}"

        # Making contig names
        paste \
        <(grep "^>" "$TMP_SORTED" | sed 's/>//') \
        <(grep "^>" "$RENAMED_FASTA" | sed 's/>//') \
        > "${SAMPLE_MEDAKA_DIR}/contig_name_map.tsv"

        rm -f "${TMP_SORTED}"
        # rm -f "${SAMPLE_MEDAKA_DIR}/consensus.fasta"

        require_file "${RENAMED_FASTA}" "renamed Medaka fasta for ${SAMPLE_ID}"

        FINAL_INPUT_FASTA="${RENAMED_FASTA}"


        ############################################
        # DNAAPLER
        ############################################

        if [[ "${RUN_DNAAPLER}" == "true" ]]; then
            echo "[INFO] Running dnaapler"

            DNAAPLER_SAMPLE_DIR="${DNAAPLER_DIR}/${SAMPLE_ID}"
            DNAAPLER_OUT="${DNAAPLER_SAMPLE_DIR}/dnaapler_reoriented.fasta"
            DNAAPLER_CLEAN="${DNAAPLER_SAMPLE_DIR}/${SAMPLE_ID}.dnaapler.clean.fasta"
            DNAAPLER_META="${DNAAPLER_SAMPLE_DIR}/${SAMPLE_ID}.dnaapler.metadata.tsv"

            if [[ "${RESUME}" == "true" && -f "${DNAAPLER_CLEAN}" ]]; then
                echo "[SKIP] dnaapler ${SAMPLE_ID}"
            else
                run_dnaapler "${RENAMED_FASTA}" "${DNAAPLER_SAMPLE_DIR}"

                require_file "${DNAAPLER_OUT}" "dnaapler output for ${SAMPLE_ID}"
                extract_dnaapler_metadata "${DNAAPLER_OUT}" "${DNAAPLER_META}"
                clean_dnaapler_headers "${DNAAPLER_OUT}" "${DNAAPLER_CLEAN}"
            fi

            rm -rf "${DNAAPLER_SAMPLE_DIR}/logs" 2>/dev/null || true
            rm -f "${DNAAPLER_SAMPLE_DIR}/*_output.txt"
            rm -f "${DNAAPLER_SAMPLE_DIR}/dnaapler_*.log"

            FINAL_INPUT_FASTA="${DNAAPLER_CLEAN}"
        else
            FINAL_INPUT_FASTA="${RENAMED_FASTA}"
        fi

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
            run_pypolca "${FINAL_INPUT_FASTA}" "${SHORT_R1_PATH}" "${SHORT_R2_PATH}" "${SAMPLE_ID}" "${SAMPLE_PYPOLCA_DIR}"

            PYPOLCA_FASTA="${SAMPLE_PYPOLCA_DIR}/${SAMPLE_ID}_corrected.fasta"
            require_file "${PYPOLCA_FASTA}" "pypolca polished fasta for ${SAMPLE_ID}"

            FINAL_INPUT_FASTA="${PYPOLCA_FASTA}"
        fi


        ############################################
        # GENERATE CONTIG STATS
        ############################################

        if [[ "${RUN_DNAAPLER}" != "true" || ! -f "${DNAAPLER_META}" ]]; then
            echo "[INFO] No dnaapler metadata, using empty file"
            DNAAPLER_META="${SAMPLE_MEDAKA_DIR}/empty_dnaapler.tsv"
            : > "${DNAAPLER_META}"
        fi

        if [[ "${GENERATE_CONTIG_STATS}" == "true" ]]; then

            echo "[INFO] Generating contig stats for ${SAMPLE_ID}"

            if [[ -z "${MIN_COV:-}" ]]; then
                MIN_COV=0
            fi

            conda run -n "${ENV_ANALYSIS}" python3 \
            "${SCRIPT_CONTIGS_STATS}" \
            --sample "${SAMPLE_ID}" \
            --map "${SAMPLE_MEDAKA_DIR}/contig_name_map.tsv" \
            --meta "${DNAAPLER_META}" \
            --mosdepth "${COV_PREFIX}.mosdepth.summary.txt" \
            --min-cov "${MIN_COV}" \
            --out "${COVERAGE_DIR}/${SAMPLE_ID}.contig.stats.tsv"
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
done < <(tail -n +2 "${SAMPLES}")

# generate a merge file with all the contig stats
OUT_ALL="${COVERAGE_DIR}/all_samples.contigs.stats.tsv"

files=("${COVERAGE_DIR}"/*.contig.stats.tsv)

if [[ ! -e "${files[0]}" ]]; then
    echo "[WARNING] No contig stats files found"
else
    first_file="${files[0]}"
    head -n1 "$first_file" > "$OUT_ALL"
    tail -n +2 -q "${files[@]}" >> "$OUT_ALL"
fi

echo
echo "[INFO] ASSEMBLY module completed"