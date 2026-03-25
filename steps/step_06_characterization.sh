#!/usr/bin/env bash
set -euo pipefail

############################################
# PATHS AND CONFIG
############################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
SCRIPTS_DIR="${ROOT_DIR}/scripts"
SCRIPT_SUMMARY="${SCRIPTS_DIR}/amrfinder_unique_symbols_by_type_and_class.py"
SCRIPT_MATRIX="${SCRIPTS_DIR}/amrfinder_presence_absence_by_type.py"
SCRIPT_HEATMAP="${SCRIPTS_DIR}/plot_amrfinder_heatmap_amr_virulence.py"
SCRIPT_PLASMIDFINDER_SUMMARY="${SCRIPTS_DIR}/summarize_plasmidfinder.py"

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
BAKTA_DIR="${ROOT_DIR}/${BAKTA_DIR}"

AMRFINDER_INDIVIDUAL_DIR="${ROOT_DIR}/${AMRFINDER_INDIVIDUAL_DIR}"
AMRFINDER_SUMMARY_DIR="${ROOT_DIR}/${AMRFINDER_SUMMARY_DIR}"

PLASMIDFINDER_INDIVIDUAL_DIR="${ROOT_DIR}/${PLASMIDFINDER_INDIVIDUAL_DIR}"
PLASMIDFINDER_SUMMARY_DIR="${ROOT_DIR}/${PLASMIDFINDER_SUMMARY_DIR}"

MOBSUITE_INDIVIDUAL_DIR="${ROOT_DIR}/${MOBSUITE_INDIVIDUAL_DIR}"
MOBSUITE_SUMMARY_DIR="${ROOT_DIR}/${MOBSUITE_SUMMARY_DIR}"

CHARACTERIZATION_LOG_DIR="${ROOT_DIR}/${LOG_DIR}/characterization"

mkdir -p "${AMRFINDER_INDIVIDUAL_DIR}"
mkdir -p "${AMRFINDER_SUMMARY_DIR}"
mkdir -p "${PLASMIDFINDER_INDIVIDUAL_DIR}"
mkdir -p "${PLASMIDFINDER_SUMMARY_DIR}"
mkdir -p "${MOBSUITE_INDIVIDUAL_DIR}"
mkdir -p "${MOBSUITE_SUMMARY_DIR}"
mkdir -p "${CHARACTERIZATION_LOG_DIR}"

############################################
# CHECK CONDA ENVS
############################################

if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_CHARACTERIZATION}"; then
    echo "ERROR: conda environment not found: ${ENV_CHARACTERIZATION}"
    exit 1
fi

if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_ANALYSIS}"; then
    echo "ERROR: conda environment not found: ${ENV_ANALYSIS}"
    exit 1
fi

############################################
# FUNCTIONS
############################################

run_amrfinder() {
    local sample_id="$1"
    local assembly_fasta="$2"
    local bakta_sample_dir="$3"
    local outfile="$4"

    conda run --no-capture-output -n "${ENV_CHARACTERIZATION}" \
        amrfinder \
        --plus \
        --threads "${AMRFINDER_THREADS}" \
        --name "${sample_id}" \
        --database "${AMRFINDER_DB}" \
        -n "${assembly_fasta}" \
        -p "${bakta_sample_dir}/${sample_id}.faa" \
        -g "${bakta_sample_dir}/${sample_id}.gff3" \
        -a bakta \
        --coverage_min "${AMRFINDER_COVERAGE_MIN}" \
        > "${outfile}"
}

run_plasmidfinder() {
    local sample_id="$1"
    local assembly_fasta="$2"
    local outdir="$3"

    mkdir -p "${outdir}"

    conda run --no-capture-output -n "${ENV_CHARACTERIZATION}" \
        plasmidfinder.py \
        -p "${PLASMIDFINDER_DB}" \
        -x \
        -t 0.90 \
        -l 0.85 \
        -i "${assembly_fasta}" \
        -o "${outdir}"
}

run_mobtyper() {
    local sample_id="$1"
    local assembly_fasta="$2"
    local outfile="$3"

    conda run --no-capture-output -n "${ENV_CHARACTERIZATION}" \
        mob_typer \
        --database_directory "${MOBSUITE_DB}" \
        -n "${MOBSUITE_THREADS}" \
        --multi \
        --infile "${assembly_fasta}" \
        --sample_id "${sample_id}" \
        --out_file "${outfile}"
}

export ENV_CHARACTERIZATION ENV_ANALYSIS
export INPUT_FASTA_DIR BAKTA_DIR
export AMRFINDER_INDIVIDUAL_DIR AMRFINDER_SUMMARY_DIR PLASMIDFINDER_INDIVIDUAL_DIR PLASMIDFINDER_SUMMARY_DIR
export MOBSUITE_INDIVIDUAL_DIR MOBSUITE_SUMMARY_DIR
export CHARACTERIZATION_LOG_DIR
export -f run_amrfinder run_plasmidfinder run_mobtyper

############################################
# MAIN
############################################

echo "=================================================="
echo "STEP 06: CHARACTERIZATION"
echo "Samples file              : ${SAMPLES}"
echo "Input fasta dir           : ${INPUT_FASTA_DIR}"
echo "Bakta dir                 : ${BAKTA_DIR}"
echo "AMRFinder individual dir  : ${AMRFINDER_INDIVIDUAL_DIR}"
echo "PlasmidFinder indiv dir   : ${PLASMIDFINDER_INDIVIDUAL_DIR}"
echo "MOBsuite individual dir   : ${MOBSUITE_INDIVIDUAL_DIR}"
echo "=================================================="

echo "[INFO] Step A: GENOMIC CHARACTERIZATION"

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
    BAKTA_SAMPLE_DIR="${BAKTA_DIR}/${SAMPLE_ID}"
    require_file "${ASSEMBLY_FASTA}" "assembly FASTA for ${SAMPLE_ID}"
    require_dir "${BAKTA_SAMPLE_DIR}" "Bakta annotation dir for ${SAMPLE_ID}"

    AMRFINDER_OUTFILE="${AMRFINDER_INDIVIDUAL_DIR}/${SAMPLE_ID}.${AMRFINDER_PREFIX}.result.tsv"
    PLASMIDFINDER_SAMPLE_DIR="${PLASMIDFINDER_INDIVIDUAL_DIR}/${SAMPLE_ID}"
    MOBSUITE_OUTFILE="${MOBSUITE_INDIVIDUAL_DIR}/${SAMPLE_ID}.${MOBSUITE_PREFIX}.result.tsv"

    SAMPLE_LOG="${CHARACTERIZATION_LOG_DIR}/${SAMPLE_ID}.characterization.log"

    if [[ ! -f "${ASSEMBLY_FASTA}" ]]; then
        echo "ERROR: final assembly not found for ${SAMPLE_ID}: ${ASSEMBLY_FASTA}"
        exit 1
    fi

    if [[ ! -d "${BAKTA_SAMPLE_DIR}" ]]; then
        echo "ERROR: Bakta annotation dir not found for ${SAMPLE_ID}: ${BAKTA_SAMPLE_DIR}"
        exit 1
    fi

    {
        echo "[INFO] Sample: ${SAMPLE_ID}"
        echo "[INFO] Assembly: ${ASSEMBLY_FASTA}"
        echo "[INFO] Bakta dir: ${BAKTA_SAMPLE_DIR}"

        echo "[INFO] AMRFinderPlus: ${SAMPLE_ID}"
        run_amrfinder "${SAMPLE_ID}" "${ASSEMBLY_FASTA}" "${BAKTA_SAMPLE_DIR}" "${AMRFINDER_OUTFILE}"
        echo "[DONE] AMRFinderPlus: ${SAMPLE_ID}"

        echo "[INFO] PlasmidFinder: ${SAMPLE_ID}"
        run_plasmidfinder "${SAMPLE_ID}" "${ASSEMBLY_FASTA}" "${PLASMIDFINDER_SAMPLE_DIR}"
        echo "[DONE] PlasmidFinder: ${SAMPLE_ID}"

        echo "[INFO] MOB-typer: ${SAMPLE_ID}"
        run_mobtyper "${SAMPLE_ID}" "${ASSEMBLY_FASTA}" "${MOBSUITE_OUTFILE}"
        echo "[DONE] MOB-typer: ${SAMPLE_ID}"

    } > "${SAMPLE_LOG}" 2>&1

    echo "Done: ${SAMPLE_ID}"
done

############################################
# STEP B: AMRFINDER COMBINED TABLE
############################################

echo
echo "[INFO] Step B: combine -> amrfinder.result.combined.tsv"

cd "${AMRFINDER_INDIVIDUAL_DIR}"

if compgen -G "*.${AMRFINDER_PREFIX}.result.tsv" > /dev/null; then
    first=$(ls *."${AMRFINDER_PREFIX}".result.tsv | head -1)
    head -1 "${first}" > "${AMRFINDER_SUMMARY_DIR}/amrfinder.result.combined.tsv"
    tail -n +2 -q *."${AMRFINDER_PREFIX}".result.tsv >> "${AMRFINDER_SUMMARY_DIR}/amrfinder.result.combined.tsv"
    sed -i "s/'//g" "${AMRFINDER_SUMMARY_DIR}/amrfinder.result.combined.tsv"
fi

############################################
# STEP C: AMRFINDER SUMMARY TABLE
############################################

echo "[INFO] Step C: summary table (unique symbols by type/class)"

conda run --no-capture-output -n "${ENV_ANALYSIS}" python3 "${SCRIPT_SUMMARY}" \
    -i "${AMRFINDER_SUMMARY_DIR}/amrfinder.result.combined.tsv" \
    -o "${AMRFINDER_SUMMARY_DIR}/amrfinder.summary.unique_symbols.tsv" \
    --include-class-counts \
    --include-type-lists \
    --include-class-lists

############################################
# STEP D: PRESENCE/ABSENCE MATRICES
############################################

echo "[INFO] Step D: presence/absence matrices per Type"

conda run --no-capture-output -n "${ENV_ANALYSIS}" python3 "${SCRIPT_MATRIX}" \
    -i "${AMRFINDER_SUMMARY_DIR}/amrfinder.result.combined.tsv" \
    -o "${AMRFINDER_SUMMARY_DIR}/amrfinder" \
    --types AMR,STRESS,VIRULENCE

############################################
# STEP E: HEATMAP
############################################

echo "[INFO] Step E: heatmap AMR+VIRULENCE (PNG+SVG)"

conda run --no-capture-output -n "${ENV_ANALYSIS}" python3 "${SCRIPT_HEATMAP}" \
    --amr "${AMRFINDER_SUMMARY_DIR}/amrfinder.AMR.presence_absence.tsv" \
    --virulence "${AMRFINDER_SUMMARY_DIR}/amrfinder.VIRULENCE.presence_absence.tsv" \
    --top-amr 30 \
    --top-vir 30 \
    --order alphabetical \
    -o "${AMRFINDER_SUMMARY_DIR}/heatmap_AMR_VIRULENCE_top30"

############################################
# STEP F: PLASMIDFINDER SUMMARY
############################################

echo "[INFO] Step F: summarize plasmidfinder reports"

conda run --no-capture-output -n "${ENV_ANALYSIS}" python3 "${SCRIPT_PLASMIDFINDER_SUMMARY}" \
    "${PLASMIDFINDER_INDIVIDUAL_DIR}" \
    "${PLASMIDFINDER_SUMMARY_DIR}/plasmidfinder_inc_types.tsv"

############################################
# STEP G: MOB-TYPER COMBINED TABLE
############################################

echo "[INFO] Step G: combine -> mobtyper.result.combined.tsv"

cd "${MOBSUITE_INDIVIDUAL_DIR}"

if compgen -G "*.${MOBSUITE_PREFIX}.result.tsv" > /dev/null; then
    second=$(ls *."${MOBSUITE_PREFIX}".result.tsv | head -1)
    head -1 "${second}" > "${MOBSUITE_SUMMARY_DIR}/mobtyper.result.combined.tsv"
    tail -n +2 -q *."${MOBSUITE_PREFIX}".result.tsv >> "${MOBSUITE_SUMMARY_DIR}/mobtyper.result.combined.tsv"
fi

cd "${ROOT_DIR}"

echo "[INFO] CHARACTERIZATION module completed"
echo "[DONE] Outputs written to: ${ROOT_DIR}/${CHARACTERIZATION_RESULTS_DIR}"