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
# STATUS FILE
############################################

STATUS_FILE="${ROOT_DIR}/results/pipeline_status.tsv"
init_status_file

export STATUS_FILE
export UTILS_FILE

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

    echo "--------------------------------------------------"
    echo "Processing sample: ${SAMPLE_ID}"
    echo "asm_type             : ${ASM_TYPE}"
    echo "expected_genome_size : ${EXPECTED_GENOME_SIZE}"
    echo "--------------------------------------------------"

    # --- FILTER FOR EMPTY/FAILED SAMPLES ---
    # Checking if the sample failed in trimming step
    if grep -qE "${SAMPLE_ID}.*trimming.*FAILED" "${STATUS_FILE}"; then
        log_warn "Skipping ${SAMPLE_ID}: already FAILED in trimming step."
        continue
    fi

    #Status per file
    status_start "$SAMPLE_ID" "assembly"

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


    # =========================
    # Define expected outputs
    # =========================
    CLEAN_FASTQ="${READS_DIR}/${SAMPLE_ID}.${TRIMMER}.${FILTLONG_PREFIX}.fastq.gz"

    SAMPLE_FLYE_DIR="${FLYE_DIR}/${SAMPLE_ID}"
    SAMPLE_MEDAKA_DIR="${MEDAKA_DIR}/${SAMPLE_ID}"

    SAMPLE_LOG="${ASSEMBLY_LOG_DIR}/${SAMPLE_ID}.assembly.log"

    FINAL_ASSEMBLY="${FINAL_ASM_DIR}/${SAMPLE_ID}.${FINAL_ASM_PREFIX}.fasta"

    COV_PREFIX="${COVERAGE_DIR}/${SAMPLE_ID}.coverage"
    TMP_BAM="$COVERAGE_DIR/${SAMPLE_ID}.tmp.bam"
    BAM="${COVERAGE_DIR}/${SAMPLE_ID}.bam"

    # =========================
    # Checking inputs
    # =========================
    # BEFORE RUNNING ANYTHING
    MIN_READS=15000

    # Check if nanostat report exist
    CLEAN_STATS="${QC_CLEAN_DIR}/${SAMPLE_ID}.${TRIMMER}.clean.txt"

    if [[ ! -f "${CLEAN_STATS}" ]]; then
        log_warn "Skipping ${SAMPLE_ID}: No QC stats found (likely failed in step 01)"
        continue
    fi

    #Checking number of reads before
    if ! check_min_reads_nanostat "${CLEAN_STATS}" "${MIN_READS}"; then
        log_warn "Skipping ${SAMPLE_ID} (low reads: < ${MIN_READS})"
        status_skip "$SAMPLE_ID" "assembly" "low_reads"
        continue
    fi

    # RESUME LOGIC
    if should_skip_sample "${FINAL_ASSEMBLY}"; then
        status_skip "$SAMPLE_ID" "assembly" "resume_existing_output"
        continue
    fi

    mkdir -p "${SAMPLE_FLYE_DIR}"
    mkdir -p "${SAMPLE_MEDAKA_DIR}"

    echo "[INFO] Sample: ${SAMPLE_ID}"
    echo "[INFO] Clean FASTQ: ${CLEAN_FASTQ}"
    echo "[INFO] Genome size: ${EXPECTED_GENOME_SIZE}"

    #RUN FLYE ASSEMBLY BLOCK
    {   
        echo "[INFO] Running Flye"
        # Running Flye (safe mode)
        run_safe "flye" run_flye "${CLEAN_FASTQ}" "${EXPECTED_GENOME_SIZE}" "${SAMPLE_FLYE_DIR}"

        # Define output
        FLYE_ASM="${SAMPLE_FLYE_DIR}/assembly.fasta"

        # Rename outputs
        cd "${SAMPLE_FLYE_DIR}"

        mv assembly.fasta "${SAMPLE_ID}.${FLYE_PREFIX}.fasta"
        mv assembly_graph.gfa "${SAMPLE_ID}.${FLYE_PREFIX}.gfa"
        mv flye.log "${SAMPLE_ID}.${FLYE_PREFIX}.log"
        mv assembly_info.txt "${SAMPLE_ID}.${FLYE_PREFIX}.info.txt"

        echo "[INFO] Flye finished"
    } >> "${SAMPLE_LOG}" 2>&1

    FLYE_FASTA="${SAMPLE_FLYE_DIR}/${SAMPLE_ID}.${FLYE_PREFIX}.fasta"

    # Checking flye output
    if ! check_file_not_empty "${FLYE_FASTA}"; then
        log_warn "Flye output file empty for ${SAMPLE_ID}"
        status_fail "$SAMPLE_ID" "assembly" "flye_empty_assembly"
        continue
    fi

    #RUN MEDAKA ASSEMBLY BLOCK
    {    
        echo "[INFO] Running Medaka"
        # Running Medaka (safe mode)
        run_safe "medaka" run_medaka "${CLEAN_FASTQ}" "${FLYE_FASTA}" "${SAMPLE_MEDAKA_DIR}"

        # Define output
        MEDAKA_ASM="${SAMPLE_MEDAKA_DIR}/consensus.fasta"

        echo "[INFO] Medaka finished"
    } >> "${SAMPLE_LOG}" 2>&1

    # Checking medaka output
    if ! check_file_not_empty "${MEDAKA_ASM}"; then
        log_warn "Medaka output file empty for ${SAMPLE_ID}"
        status_fail "$SAMPLE_ID" "assembly" "medaka_empty_output"
        continue
    fi

    #RUN CONTIG COVERAGE BLOCK
    {  
        echo "[INFO] Calculating contig coverage"

        if [[ "${RESUME}" == "true" && -f "${COV_PREFIX}.mosdepth.summary.txt" ]]; then
            echo "[SKIP] coverage ${SAMPLE_ID}"
        else
            # Running Mosdepth (safe mode)
            run_safe "mosdepth" run_mosdepth \
            "${SAMPLE_MEDAKA_DIR}/consensus.fasta" \
            "${CLEAN_FASTQ}" \
            "${TMP_BAM}" \
            "${BAM}" \
            "${COV_PREFIX}"
            
            # Removing temporary files (BAM)
            rm -f "${TMP_BAM}" "${BAM}" "${BAM}.bai"
        fi

        echo "[INFO] Contig coverage finished"
    } >> "${SAMPLE_LOG}" 2>&1

    # Checking mosdepth output
    if ! check_file_not_empty "${COV_PREFIX}.mosdepth.summary.txt"; then
        log_warn "Mosdepth output empty for ${SAMPLE_ID}"
    fi

    #RUN FILTER CONTIGS BLOCK
    
    FILTERED_FASTA="${SAMPLE_MEDAKA_DIR}/consensus.filtered.fasta"
    FINAL_INPUT_FASTA="${SAMPLE_MEDAKA_DIR}/consensus.fasta"
    SUMMARY="${COV_PREFIX}.mosdepth.summary.txt"

    # Checking input for filtering
    if ! check_file_not_empty "${SUMMARY}"; then
        log_warn "Missing mosdepth summary for ${SAMPLE_ID}, skipping filtering"
        status_skip "$SAMPLE_ID" "assembly" "missing_coverage"
        continue
    fi

    {
        echo "[INFO] Filtering contigs by coverage"

        if [[ "${FILTER_CONTIGS}" == "true" ]]; then

            # GET CHROM COVERAGE
            CHROM_COV=$(awk '
            NR > 1 && $1 !~ /total/ {
                print $4, $2
            }' "${SUMMARY}" | sort -k2,2nr | head -n1 | awk '{print $1}')

            echo "[INFO] Chromosome coverage: ${CHROM_COV}"
            MIN_COV=$(awk -v c="$CHROM_COV" -v f="$MIN_COV_FACTOR" 'BEGIN{print c*f}')
            echo "[INFO] Min coverage threshold: ${MIN_COV}"

            # SELECT CONTIGS
            awk -v mincov="$MIN_COV" -v minlen="$MIN_CONTIG_LENGTH" '
            NR > 1 && $1 !~ /total/ {
                if ($4 >= mincov && $2 >= minlen) print $1
                }' "${SUMMARY}" > "${SAMPLE_MEDAKA_DIR}/contigs.keep.txt"

            # FILTER FASTA
            run_safe "filter_contigs" \
                conda run -n "${ENV_ASSEMBLY}" seqkit grep \
                -f "${SAMPLE_MEDAKA_DIR}/contigs.keep.txt" \
                "${SAMPLE_MEDAKA_DIR}/consensus.fasta" \
                > "${FILTERED_FASTA}"
        else
            FINAL_INPUT_FASTA="${SAMPLE_MEDAKA_DIR}/consensus.fasta"
        fi

        echo "[INFO] Filter contigs finished"
    } >> "${SAMPLE_LOG}" 2>&1

    # Checking filtering output
    if [[ ! -s "${SAMPLE_MEDAKA_DIR}/contigs.keep.txt" ]]; then
        log_warn "No contigs passed filtering for ${SAMPLE_ID}, using unfiltered assembly"
        status_skip "$SAMPLE_ID" "assembly" "no_contigs_after_filtering"
        continue
    fi

    if ! check_file_not_empty "${FILTERED_FASTA}"; then
        log_warn "Filtering failed for ${SAMPLE_ID}"
        FINAL_INPUT_FASTA="${SAMPLE_MEDAKA_DIR}/consensus.fasta"
    fi

    #RUN SORTING AND RENAMING CONTIGS BLOCK
    
    TMP_SORTED="${SAMPLE_MEDAKA_DIR}/consensus.sorted.fasta"
    RENAMED_FASTA="${SAMPLE_MEDAKA_DIR}/${SAMPLE_ID}.${MEDAKA_PREFIX}.fasta"

    # Checking input
    if ! check_file_not_empty "${FINAL_INPUT_FASTA}"; then
        log_warn "Input FASTA missing before renaming for ${SAMPLE_ID}"
        continue
    fi

    {
        echo "[INFO] Sorting and renaming contigs"

        # RUN RENAME (SAFE)
        run_safe "rename_contigs" sort_and_rename_contigs \
            "${FINAL_INPUT_FASTA}" \
            "${RENAMED_FASTA}" \
            "${SAMPLE_ID}" \
            "${TMP_SORTED}"

        # GENERATE MAPPING
        if [[ $(grep -c "^>" "$TMP_SORTED") -ne $(grep -c "^>" "$RENAMED_FASTA") ]]; then
            log_warn "Mismatch in contig counts for ${SAMPLE_ID}, skipping mapping"
        else
            paste \
            <(grep "^>" "$TMP_SORTED" | sed 's/>//') \
            <(grep "^>" "$RENAMED_FASTA" | sed 's/>//') \
            > "${SAMPLE_MEDAKA_DIR}/contig_name_map.tsv"
        fi

        # CLEANUP
        rm -f "${TMP_SORTED}" 2>/dev/null || true

        echo "[INFO] sorting and renaming contigs finished"
    } >> "${SAMPLE_LOG}" 2>&1

    FINAL_INPUT_FASTA="${RENAMED_FASTA}"

    # Checking output
    if ! check_file_not_empty "${RENAMED_FASTA}"; then
        log_warn "Renamed FASTA missing for ${SAMPLE_ID}"
        continue
    fi

    #RUN DNAAPLER BLOCK

    DNAAPLER_SAMPLE_DIR="${DNAAPLER_DIR}/${SAMPLE_ID}"
    DNAAPLER_OUT="${DNAAPLER_SAMPLE_DIR}/dnaapler_reoriented.fasta"
    DNAAPLER_CLEAN="${DNAAPLER_SAMPLE_DIR}/${SAMPLE_ID}.dnaapler.clean.fasta"
    DNAAPLER_META="${DNAAPLER_SAMPLE_DIR}/${SAMPLE_ID}.dnaapler.metadata.tsv"

    {
        if [[ "${RUN_DNAAPLER}" == "true" ]]; then
            echo "[INFO] Running dnaapler"

            if [[ "${RESUME}" == "true" && -f "${DNAAPLER_CLEAN}" ]]; then
                echo "[SKIP] dnaapler ${SAMPLE_ID}"
            else
                if ! run_safe "dnaapler" run_dnaapler "${RENAMED_FASTA}" "${DNAAPLER_SAMPLE_DIR}"; then
                    log_warn "dnaapler failed for ${SAMPLE_ID}, continue to next file"
                    FINAL_INPUT_FASTA="${RENAMED_FASTA}"
                else

                    #Check output
                    if ! check_file_not_empty "${DNAAPLER_OUT}"; then
                        log_warn "Dnaapler output empty for ${SAMPLE_ID}"
                        FINAL_INPUT_FASTA="${RENAMED_FASTA}"
                    else
                        extract_dnaapler_metadata "${DNAAPLER_OUT}" "${DNAAPLER_META}"
                        clean_dnaapler_headers "${DNAAPLER_OUT}" "${DNAAPLER_CLEAN}"
                    fi
                fi
            fi
        fi
        
        # removing temporary files
        rm -rf "${DNAAPLER_SAMPLE_DIR}/logs" 2>/dev/null || true
        rm -f "${DNAAPLER_SAMPLE_DIR}"/*_output.txt 2>/dev/null || true
        rm -f "${DNAAPLER_SAMPLE_DIR}"/dnaapler_*.log 2>/dev/null || true

        echo "[INFO] Dnaapler finished"
    } >>"${SAMPLE_LOG}" 2>&1

    #Final input
    if [[ -f "${DNAAPLER_CLEAN}" ]]; then
        FINAL_INPUT_FASTA="${DNAAPLER_CLEAN}"
    else
        FINAL_INPUT_FASTA="${RENAMED_FASTA}"
    fi


    #RUN HYBRID ASSEMBLY POLISH BLOCK

    {
        if [[ "${ASM_TYPE}" == "hybrid" ]]; then
            SHORT_R1_PATH="${RAW_DIR}/${SHORT_R1}"
            SHORT_R2_PATH="${RAW_DIR}/${SHORT_R2}"
            SAMPLE_PYPOLCA_DIR="${SAMPLE_MEDAKA_DIR}/pypolca"

            # CHECK SHORT READS
            if ! check_file_not_empty "${SHORT_R1_PATH}" || ! check_file_not_empty "${SHORT_R2_PATH}"; then
                log_warn "Missing short reads for ${SAMPLE_ID}, skipping hybrid polishing"
                # fallback → keep current assembly
            else

                echo "[INFO] Running pypolca"
                # RUN SAFE
                if ! run_safe "pypolca" run_pypolca \
                    "${FINAL_INPUT_FASTA}" \
                    "${SHORT_R1_PATH}" \
                    "${SHORT_R2_PATH}" \
                    "${SAMPLE_ID}" \
                    "${SAMPLE_PYPOLCA_DIR}"
                then
                    log_warn "pypolca failed for ${SAMPLE_ID}, using unpolished assembly"
                else
                    PYPOLCA_FASTA="${SAMPLE_PYPOLCA_DIR}/${SAMPLE_ID}_corrected.fasta"

                    # CHECK OUTPUT
                    if ! check_file_not_empty "${PYPOLCA_FASTA}"; then
                        log_warn "pypolca output empty for ${SAMPLE_ID}, skipping"
                    else
                        FINAL_INPUT_FASTA="${PYPOLCA_FASTA}"
                    fi
                fi
            fi
        fi
        echo "[INFO] hybrid assembly polish finished"            
    } >> "${SAMPLE_LOG}" 2>&1

    #RUN GENERATE CONTIG STATS BLOCK

    # Checking input files
    #Check input contig_name
    if ! check_file_not_empty "${SAMPLE_MEDAKA_DIR}/contig_name_map.tsv"; then
        log_warn "Missing contig map for ${SAMPLE_ID}"
        continue
    fi

    #Check empty file coverage
    if ! check_file_not_empty "${COV_PREFIX}.mosdepth.summary.txt"; then
        log_warn "Missing mosdepth summary for ${SAMPLE_ID}"
        continue
    fi

    # Minimum coverage value
    if [[ -z "${MIN_COV:-}" ]]; then
        log_warn "MIN_COV not set, using default 1"
        MIN_COV=1
    fi

    {
        if [[ "${RUN_DNAAPLER}" != "true" || ! -f "${DNAAPLER_META}" ]]; then
            echo "[INFO] No dnaapler metadata, using empty file"
            DNAAPLER_META="${SAMPLE_MEDAKA_DIR}/empty_dnaapler.tsv"
            [[ ! -f "$DNAAPLER_META" ]] && : > "${DNAAPLER_META}"
        fi

        if [[ "${GENERATE_CONTIG_STATS}" == "true" ]]; then

            echo "[INFO] Generating contig stats for ${SAMPLE_ID}"

            #Run contig stats python
            if ! run_safe "contig_stats" \
                conda run -n "${ENV_ANALYSIS}" python3 \
                "${SCRIPT_CONTIGS_STATS}" \
                --sample "${SAMPLE_ID}" \
                --map "${SAMPLE_MEDAKA_DIR}/contig_name_map.tsv" \
                --meta "${DNAAPLER_META}" \
                --mosdepth "${COV_PREFIX}.mosdepth.summary.txt" \
                --min-cov "${MIN_COV}" \
                --out "${COVERAGE_DIR}/${SAMPLE_ID}.contig.stats.tsv"
            then
                log_warn "contig stats failed for ${SAMPLE_ID}"
            fi
        fi

        echo "[INFO] Sample ${SAMPLE_ID} completed successfully"
    } >> "${SAMPLE_LOG}" 2>&1

    ############################################
    # FINAL ASSEMBLY
    ############################################

    FINAL_OUT="${FINAL_ASM_DIR}/${SAMPLE_ID}.${FINAL_ASM_PREFIX}.fasta"

    echo "[INFO] Final assembly written to: ${FINAL_OUT}"
    
    # CHECK INPUT
    if ! check_file_not_empty "${FINAL_INPUT_FASTA}"; then
        log_warn "Final input FASTA missing for ${SAMPLE_ID}, skipping final assembly"
        continue
    fi

    # COPY (SAFE)
    if ! run_safe "final_copy" cp "${FINAL_INPUT_FASTA}" "${FINAL_OUT}"; then
        log_warn "Failed to write final assembly for ${SAMPLE_ID}"
        continue
    fi

    # VALIDATE OUTPUT
    if ! check_file_not_empty "${FINAL_OUT}"; then
        log_warn "Final assembly empty for ${SAMPLE_ID}"
        continue
    fi

    #FINAL LOG STATUS
    status_ok "$SAMPLE_ID" "assembly"

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

#Final log info
log_info "Assembly step completed for all samples"
awk '$2=="assembly"' "${STATUS_FILE}"

echo
echo "[INFO] ASSEMBLY module completed"