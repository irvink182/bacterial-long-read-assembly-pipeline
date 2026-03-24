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
# INPUT / OUTPUT DIRS
############################################

RAW_DIR="${ROOT_DIR}/${RAW_READS_DIR}"

RESULT_DIR="${ROOT_DIR}/results/reads"
PORECHOP_DIR="${RESULT_DIR}/porechop_reads"
CLEAN_DIR="${RESULT_DIR}/clean_reads"

QC_DIR="${RESULT_DIR}/qc"
QC_RAW_DIR="${QC_DIR}/nanostat_raw"
QC_CLEAN_DIR="${QC_DIR}/nanostat_clean"
MULTIQC_DIR="${QC_DIR}/multiqc"

SYLPH_DIR="${RESULT_DIR}/sylph"
SYLPH_SKETCH_DIR="${SYLPH_DIR}/sylph_sketches"
SYLPH_TABLE_DIR="${SYLPH_DIR}/sylph_tables"

LOG_DIR="${ROOT_DIR}/logs/trimming"

mkdir -p "${PORECHOP_DIR}"
mkdir -p "${CLEAN_DIR}"
mkdir -p "${QC_RAW_DIR}"
mkdir -p "${QC_CLEAN_DIR}"
mkdir -p "${FASTPLONG_QC_DIR}"
mkdir -p "${FILTLONG_QC_DIR}"
mkdir -p "${PORECHOP_QC_DIR}"
mkdir -p "${QC_RAW_DIR}"
mkdir -p "${QC_CLEAN_DIR}"
mkdir -p "${MULTIQC_RAW_DIR}"
mkdir -p "${MULTIQC_CLEAN_DIR}"
#mkdir -p "${MULTIQC_DIR}"
mkdir -p "${SYLPH_DIR}"
mkdir -p "${SYLPH_SKETCH_DIR}"
mkdir -p "${SYLPH_TABLE_DIR}"
mkdir -p "${LOG_DIR}"

############################################
# CHECK CONDA ENV
############################################

if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_TRIMMING}"; then
    echo "ERROR: conda environment not found: ${ENV_TRIMMING}"
    exit 1
fi

############################################
# FUNCTIONS
############################################

run_nanostat() {
    local input_fastq="$1"
    local sample_name="$2"
    local outdir="$3"

    conda run -n "${ENV_TRIMMING}" NanoStat \
        --fastq "${input_fastq}" \
        --threads "${THREADS}" \
        -n "${sample_name}.txt" \
        --outdir "${outdir}"
}

run_sylph_sketch() {
    local clean_fastq="$1"
    local sample_id="$2"
    local sketch_dir="$3"

    conda run -n "${ENV_SYLPH}" sylph sketch \
        -t "${SYLPH_THREADS}" \
        "${clean_fastq}" \
        -S "${sample_id}" \
        -d "${sketch_dir}"
}

run_sylph_profile_individual() {
    local sample_id="$1"
    local sketch_dir="$2"
    local table_dir="$3"

    conda run -n "${ENV_SYLPH}" sylph profile \
        -t "${SYLPH_THREADS}" -u \
        "${sketch_dir}/${sample_id}.sylsp" \
        "${SYLPHDB1}" \
        "${SYLPHDB2}" \
        -o "${table_dir}/${sample_id}.sylph_table.tsv"
}

run_sylph_profile_general() {
    local sketch_dir="$1"
    local out_dir="$2"

    conda run -n "${ENV_SYLPH}" bash -c "
        set -euo pipefail
        sylph profile \
            -t '${SYLPH_THREADS}' -u \
            '${sketch_dir}'/*.sylsp \
            '${SYLPHDB1}' \
            '${SYLPHDB2}' \
            -o '${out_dir}/sylph_table.tsv'
    "
}

run_sylph_taxprof() {
    local out_dir="$1"

    conda run -n "${ENV_SYLPH}" sylph-tax taxprof \
        "${out_dir}/sylph_table.tsv" \
        -t "${SYLPHTAXDB1}" "${SYLPHTAXDB2}" \
        -o "${out_dir}/taxprof-"
}

run_sylph_merge_tables() {
    local out_dir="$1"

    conda run -n "${ENV_SYLPH}" bash -c "
        set -euo pipefail

        sylph-tax merge \
            '${out_dir}'/*.sylphmpa \
            --column relative_abundance \
            -o '${out_dir}/merged_abundance_file.tsv'

        grep -E 's__|clade_name' '${out_dir}/merged_abundance_file.tsv' \
        | grep -v 't__' \
        | sed 's/^.*|//g' \
        | sed 's/SRS[0-9]*-//g' \
        > '${out_dir}/merged_abundance_table_species.tsv'

        datamash transpose \
            < '${out_dir}/merged_abundance_file.tsv' \
            | sed 's/[ \t]*$//' \
            > '${out_dir}/merged_abundance_file.transposed.tsv'

        datamash transpose \
            < '${out_dir}/merged_abundance_table_species.tsv' \
            | sed 's/[ \t]*$//' \
            > '${out_dir}/merged_abundance_table_species.transposed.tsv'
    "
}

run_sylph_move_mpa() {
    local out_dir="$1"

    mkdir -p "${out_dir}/sylphmpa_files"
    mv "${out_dir}"/*.sylphmpa "${out_dir}/sylphmpa_files/" 2>/dev/null || true
}

run_sylph_plot() {
    local out_dir="$1"

    if [[ -n "${RSCRIPT_PLOT:-}" && -f "${RSCRIPT_PLOT}" ]]; then
        conda run -n "${ENV_ANALYSIS}" Rscript "${RSCRIPT_PLOT}" \
            input="${out_dir}/merged_abundance_table_species.tsv" \
            output="${out_dir}/abundance_species_plot.png" \
            title="Species Abundance" \
            min_pct=2
    else
        echo "[WARNING] RSCRIPT_PLOT not defined or file not found. Skipping plot."
    fi
}

run_fastplong() {
    local raw_fastq="$1"
    local clean_fastq="$2"
    local sample_id="$3"

    conda run -n "${ENV_TRIMMING}" fastplong \
        -i "${raw_fastq}" \
        -o "${clean_fastq}" \
        --thread "${THREADS}" \
        --length_required "${FASTPLONG_MIN_LENGTH}" \
        --json "${FASTPLONG_QC_DIR}/${sample_id}.fastplong.json" \
        --html "${FASTPLONG_QC_DIR}/${sample_id}.fastplong.html"
}

run_porechop_filtlong() {
    local raw_fastq="$1"
    local porechop_fastq="$2"
    local clean_fastq="$3"
    local sample_id="$4"

    echo "[INFO] Running porechop"

    if [[ "${PORECHOP_DISCARD_MIDDLE}" == "true" ]]; then
        conda run -n "${ENV_TRIMMING}" porechop \
            --discard_middle \
            -t "${THREADS}" \
            -i "${raw_fastq}" \
            -o "${porechop_fastq}" \
            > "${PORECHOP_QC_DIR}/${sample_id}.porechop.stdout.log" \
            2> "${PORECHOP_QC_DIR}/${sample_id}.porechop.stderr.log"
    else
        conda run -n "${ENV_TRIMMING}" porechop \
            -t "${THREADS}" \
            -i "${raw_fastq}" \
            -o "${porechop_fastq}" \
            > "${PORECHOP_QC_DIR}/${sample_id}.porechop.stdout.log" \
            2> "${PORECHOP_QC_DIR}/${sample_id}.porechop.stderr.log"
    fi

    if [[ ! -s "${porechop_fastq}" ]]; then
        echo "ERROR: porechop output was not generated for ${sample_id}"
        exit 1
    fi

    echo "[INFO] Running filtlong"

    conda run -n "${ENV_TRIMMING}" bash -c "
        set -euo pipefail
        filtlong \
            --min_length '${FILTLONG_MIN_LENGTH}' \
            --min_mean_q '${FILTLONG_MIN_MEAN_Q}' \
            --keep_percent '${FILTLONG_KEEP_PERCENT}' \
            '${porechop_fastq}' \
            2> "${FILTLONG_QC_DIR}/${sample_id}.filtlong.log" \
        | gzip > '${clean_fastq}'
    "
}

############################################
# MAIN
############################################

echo "=================================================="
echo "STEP 01: TRIMMING"
echo "Samples file      : ${SAMPLES}"
echo "Raw reads dir     : ${RAW_DIR}"
echo "Results dir       : ${RESULT_DIR}"
echo "Trimming strategy : ${TRIMMER}"
echo "Conda env         : ${ENV_TRIMMING}"
echo "=================================================="

tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SAMPLE_ID ASM_TYPE EXPECTED_GENOME_SIZE LONG_READS SHORT_R1 SHORT_R2
do
    [[ -z "${SAMPLE_ID}" ]] && continue

    echo
    echo "--------------------------------------------------"
    echo "Processing sample: ${SAMPLE_ID}"
    echo "asm_type              : ${ASM_TYPE}"
    echo "expected_genome_size  : ${EXPECTED_GENOME_SIZE}"
    echo "--------------------------------------------------"

    validate_asm_type "${ASM_TYPE}"

    if is_na "${EXPECTED_GENOME_SIZE}"; then
        log_warn "sample ${SAMPLE_ID} has empty expected_genome_size"
    else
        validate_genome_size "${EXPECTED_GENOME_SIZE}"
    fi


    if [[ "${ASM_TYPE}" == "hybrid" ]]; then
        log_info "Sample ${SAMPLE_ID} configured for hybrid assembly"
    else
        log_info "Sample ${SAMPLE_ID} configured for only long-read assembly"
    fi


    RAW_FASTQ="${RAW_DIR}/${LONG_READS}"
    require_file "${RAW_FASTQ}" "raw FASTQ for ${SAMPLE_ID}"

    PORECHOP_FASTQ="${PORECHOP_DIR}/${SAMPLE_ID}.${PORECHOP_PREFIX}.fastq.gz"
    CLEAN_FASTQ="${CLEAN_DIR}/${SAMPLE_ID}.${TRIMMER}.${FILTLONG_PREFIX}.fastq.gz"

    RAW_NANOSTAT="${QC_RAW_DIR}/${SAMPLE_ID}.nanostat.txt"
    CLEAN_NANOSTAT="${QC_CLEAN_DIR}/${SAMPLE_ID}.${TRIMMER}.nanostat.txt"

    SAMPLE_LOG="${LOG_DIR}/${SAMPLE_ID}.${TRIMMER}.trimming.log"

    {
        echo "[INFO] Sample: ${SAMPLE_ID}"
        echo "[INFO] Raw FASTQ: ${RAW_FASTQ}"
        echo "[INFO] Trimmer: ${TRIMMER}"
        echo "[INFO] Output clean FASTQ: ${CLEAN_FASTQ}"

        echo "[INFO] Running NanoStat on raw reads"
        run_nanostat "${RAW_FASTQ}" "${SAMPLE_ID}.raw" "${QC_RAW_DIR}"

        echo "[INFO] Running trimming"
        if [[ "${TRIMMER}" == "fastplong" ]]; then
            run_fastplong "${RAW_FASTQ}" "${CLEAN_FASTQ}" "${SAMPLE_ID}"

        elif [[ "${TRIMMER}" == "porechop_filtlong" ]]; then
            run_porechop_filtlong "${RAW_FASTQ}" "${PORECHOP_FASTQ}" "${CLEAN_FASTQ}" "${SAMPLE_ID}"

        else
            echo "ERROR: unknown TRIMMER value: ${TRIMMER}"
            exit 1
        fi

        if [[ ! -s "${CLEAN_FASTQ}" ]]; then
            echo "ERROR: clean FASTQ was not generated for ${SAMPLE_ID}"
            exit 1
        fi

        echo "[INFO] Running NanoStat on clean reads"
        run_nanostat "${CLEAN_FASTQ}" "${SAMPLE_ID}.${TRIMMER}.${FILTLONG_PREFIX}" "${QC_CLEAN_DIR}"

        echo "[INFO] Running sylph skecth"
        run_sylph_sketch "${CLEAN_FASTQ}" "${SAMPLE_ID}" "${SYLPH_SKETCH_DIR}"

        echo "[INFO] Running individual profile"
        run_sylph_profile_individual "${SAMPLE_ID}" "${SYLPH_SKETCH_DIR}" "${SYLPH_TABLE_DIR}"

        echo "[INFO] Sample ${SAMPLE_ID} completed successfully"
    } > "${SAMPLE_LOG}" 2>&1

    echo "Done: ${SAMPLE_ID}"
done

#SYLPH GENERAL REPORT AND PLOT ABUNDANCE TABLE
echo
echo "[INFO] Running sylph general profile"
run_sylph_profile_general "${SYLPH_SKETCH_DIR}" "${SYLPH_DIR}"

echo "[INFO] Running sylph-tax taxprof"
run_sylph_taxprof "${SYLPH_DIR}"

echo "[INFO] Merging sylph abundance tables"
run_sylph_merge_tables "${SYLPH_DIR}"

echo "[INFO] Moving .sylphmpa files"
run_sylph_move_mpa "${SYLPH_DIR}"

echo "[INFO] Plotting species abundance table"
run_sylph_plot "${SYLPH_DIR}"

#MULTIQC NANOSTAT GENERAL REPORT
echo
echo "[INFO] Running MultiQC for NanoStat raw reports"

conda run -n "${ENV_TRIMMING}" multiqc "${QC_RAW_DIR}" \
    --force \
    --cl-config "max_table_rows: 3000" \
    --outdir "${MULTIQC_RAW_DIR}" \
    > "${LOG_DIR}/multiqc_nanostat_raw.stdout.log" \
    2> "${LOG_DIR}/multiqc_nanostat_raw.stderr.log" || true

if [[ -f "${MULTIQC_RAW_DIR}/multiqc_report.html" ]]; then
    mv "${MULTIQC_RAW_DIR}/multiqc_report.html" \
       "${MULTIQC_RAW_DIR}/multiqc_nanostats_raw_report.html"
fi

if [[ -f "${MULTIQC_RAW_DIR}/multiqc_data/multiqc_nanostat.txt" ]]; then
    mv "${MULTIQC_RAW_DIR}/multiqc_data/multiqc_nanostat.txt" \
       "${MULTIQC_RAW_DIR}/multiqc_data/multiqc_nanostats_raw_report.txt"
    sed -i 's/\\.0//g' "${MULTIQC_RAW_DIR}/multiqc_data/multiqc_nanostats_raw_report.txt"
    cp "${MULTIQC_RAW_DIR}/multiqc_data/multiqc_nanostats_raw_report.txt" \
       "${MULTIQC_RAW_DIR}/"
fi

echo "[INFO] Running MultiQC for NanoStat clean reports"

conda run -n "${ENV_TRIMMING}" multiqc "${QC_CLEAN_DIR}" \
    --force \
    --cl-config "max_table_rows: 3000" \
    --outdir "${MULTIQC_CLEAN_DIR}" \
    > "${LOG_DIR}/multiqc_nanostat_clean.stdout.log" \
    2> "${LOG_DIR}/multiqc_nanostat_clean.stderr.log" || true

if [[ -f "${MULTIQC_CLEAN_DIR}/multiqc_report.html" ]]; then
    mv "${MULTIQC_CLEAN_DIR}/multiqc_report.html" \
       "${MULTIQC_CLEAN_DIR}/multiqc_nanostats_clean_report.html"
fi

if [[ -f "${MULTIQC_CLEAN_DIR}/multiqc_data/multiqc_nanostat.txt" ]]; then
    mv "${MULTIQC_CLEAN_DIR}/multiqc_data/multiqc_nanostat.txt" \
       "${MULTIQC_CLEAN_DIR}/multiqc_data/multiqc_nanostats_clean_report.txt"
    sed -i 's/\\.0//g' "${MULTIQC_CLEAN_DIR}/multiqc_data/multiqc_nanostats_clean_report.txt"
    cp "${MULTIQC_CLEAN_DIR}/multiqc_data/multiqc_nanostats_clean_report.txt" \
       "${MULTIQC_CLEAN_DIR}/"
fi

cd "${ROOT_DIR}"

echo "[INFO] TRIMMING module completed"