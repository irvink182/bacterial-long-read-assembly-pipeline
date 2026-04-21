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

#HELPER SCRIPTS
SCRIPT_SUMMARY="${ROOT_DIR}/scripts/amrfinder_unique_symbols_by_type_and_class.py"
SCRIPT_MATRIX="${ROOT_DIR}/scripts/amrfinder_presence_absence_by_type.py"
SCRIPT_HEATMAP="${ROOT_DIR}/scripts/plot_amrfinder_heatmap_amr_virulence.py"
SCRIPT_PLASMIDFINDER_SUMMARY="${ROOT_DIR}/scripts/summarize_plasmidfinder.py"


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
log_info "Starting characterization step"
[[ "${RESUME}" == "true" ]] && log_info "Resume mode enabled"

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

MLST_INDIVIDUAL_DIR="${ROOT_DIR}/${MLST_INDIVIDUAL_DIR}"
MLST_SUMMARY_DIR="${ROOT_DIR}/${MLST_SUMMARY_DIR}"

CHARACTERIZATION_LOG_DIR="${ROOT_DIR}/${LOG_DIR}/characterization"

mkdir -p "${AMRFINDER_INDIVIDUAL_DIR}"
mkdir -p "${AMRFINDER_SUMMARY_DIR}"
mkdir -p "${PLASMIDFINDER_INDIVIDUAL_DIR}"
mkdir -p "${PLASMIDFINDER_SUMMARY_DIR}"
mkdir -p "${MOBSUITE_INDIVIDUAL_DIR}"
mkdir -p "${MOBSUITE_SUMMARY_DIR}"
mkdir -p "${CHARACTERIZATION_LOG_DIR}"

############################################
# STATUS FILE
############################################

STATUS_FILE="${ROOT_DIR}/results/pipeline_status.tsv"
init_status_file

export STATUS_FILE
export UTILS_FILE

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
VALID_SAMPLES_FILE=$(mktemp)

tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r SID TYPE SIZE LONG R1 R2
do
    [[ -z "${SID}" ]] && continue

    # --- FILTERS FOR EMPTY/FAILED SAMPLES ---
    # # Checking if the sample failed in previous steps
    if grep -qE "${SID}.*(trimming|assembly|annotation).*(FAILED|SKIPPED)" "${STATUS_FILE}"; then
        continue
    fi

    # INPUT FILES PATHS
    ASSEMBLY_FASTA="${INPUT_FASTA_DIR}/${SID}.${FINAL_ASM_PREFIX}.fasta"
    BAKTA_SAMPLE_DIR="${BAKTA_DIR}/${SID}"
    BAKTA_FILE="${BAKTA_SAMPLE_DIR}/${SID}.gff3"

    # Checking if the sample has assembly, if not, SKIPPED.
    if [[ -s "${ASSEMBLY_FASTA}" && -s "${BAKTA_FILE}" ]]; then
        echo -e "${SID}\t${TYPE}\t${SIZE}\t${LONG}\t${R1}\t${R2}" >> "${VALID_SAMPLES_FILE}"
        status_start "$SID" "characterization"
    else
        status_skip "${SID}" "characterization" "missing_input_files"
        log_warn "Required input files missing for ${SID}. Skipping characterization step."
    fi
done

#MAIN LOOP
while IFS=$'\t' read -r SAMPLE_ID ASM_TYPE EXPECTED_GENOME_SIZE LONG_READS SHORT_R1 SHORT_R2
do

    [[ -z "${SAMPLE_ID}" ]] && continue

    echo "--------------------------------------------------"
    echo "Processing sample: ${SAMPLE_ID}"
    echo "--------------------------------------------------"

    # PATHS TO CHECK INPUT FILES
    ASSEMBLY_FASTA="${INPUT_FASTA_DIR}/${SAMPLE_ID}.${FINAL_ASM_PREFIX}.fasta"
    BAKTA_SAMPLE_DIR="${BAKTA_DIR}/${SAMPLE_ID}"
    BAKTA_FILE="${BAKTA_SAMPLE_DIR}/${SAMPLE_ID}.gff3"

    SAMPLE_LOG="${CHARACTERIZATION_LOG_DIR}/${SAMPLE_ID}.characterization.log"

    #OUTPUT PATHS
    AMRFINDER_OUTFILE="${AMRFINDER_INDIVIDUAL_DIR}/${SAMPLE_ID}.${AMRFINDER_PREFIX}.result.tsv"
    PLASMIDFINDER_SAMPLE_DIR="${PLASMIDFINDER_INDIVIDUAL_DIR}/${SAMPLE_ID}"
    MOBSUITE_OUTFILE="${MOBSUITE_INDIVIDUAL_DIR}/${SAMPLE_ID}.${MOBSUITE_PREFIX}.result.tsv"

    validate_asm_type "${ASM_TYPE}"

    if is_na "${EXPECTED_GENOME_SIZE}"; then
        log_warn "sample ${SAMPLE_ID} has empty expected_genome_size"
    else
        validate_genome_size "${EXPECTED_GENOME_SIZE}"
    fi

    # =========================
    # Checking inputs
    # =========================

    #Checking input assembly
    if ! check_file_not_empty "${ASSEMBLY_FASTA}"; then
        log_warn "Assembly not found: ${SAMPLE_ID}"
        status_skip "$SAMPLE_ID" "characterization" "missing_assembly"
        continue
    fi

    #Checking input annotation
    if ! check_file_not_empty "${BAKTA_FILE}"; then
        log_warn "Annotation not found: ${SAMPLE_ID}"
        status_skip "$SAMPLE_ID" "characterization" "missing_bakta"
        continue
    fi


    # =========================
    # Define outputs
    # =========================

    AMRFINDER_OUTFILE="${AMRFINDER_INDIVIDUAL_DIR}/${SAMPLE_ID}.${AMRFINDER_PREFIX}.result.tsv"
    PLASMIDFINDER_OUT="${PLASMIDFINDER_SAMPLE_DIR}/results.txt"
    MOBSUITE_OUTFILE="${MOBSUITE_INDIVIDUAL_DIR}/${SAMPLE_ID}.${MOBSUITE_PREFIX}.result.tsv"
    MLST_OUT="${MLST_INDIVIDUAL_DIR}/${SAMPLE_ID}.mlst.tsv"

    RUN_AMR=true
    RUN_PLASMID=true
    RUN_MOBSUITE=true
    RUN_MLST=true

    # =========================
    # Resume logic
    # =========================

    if [[ "${RESUME}" == "true" && -s "${AMRFINDER_OUTFILE}" ]]; then
        log_info "Skipping AMRFinder for ${SAMPLE_ID}"
        RUN_AMR=false
    fi

    if [[ "${RESUME}" == "true" && -s "${PLASMIDFINDER_OUT}" ]]; then
        log_info "Skipping PlasmidFinder for ${SAMPLE_ID}"
        RUN_PLASMID=false
    fi

    if [[ "${RESUME}" == "true" && -s "${MOBSUITE_OUTFILE}" ]]; then
        log_info "Skipping MobTyper for ${SAMPLE_ID}"
        RUN_MOBSUITE=false
    fi

    if [[ "${RESUME}" == "true" && -s "${MLST_OUT}" ]]; then
        log_info "Skipping MLST for ${SAMPLE_ID}"
        RUN_MLST=false
    fi

    # skip total
    if [[ "${RUN_AMR}" == "false" && "${RUN_PLASMID}" == "false" && "${RUN_MOBSUITE}" == "false" && "${RUN_MLST}" == "false" ]]; then
        log_info "Skipping ${SAMPLE_ID}: all characterization outputs exist"
        status_skip "$SAMPLE_ID" "characterization" "resume_existing_output"
        continue
    fi

    echo "[INFO] Sample: ${SAMPLE_ID}"
    echo "[INFO] Assembly: ${ASSEMBLY_FASTA}"
    echo "[INFO] Bakta dir: ${BAKTA_SAMPLE_DIR}"

    mkdir -p "${PLASMIDFINDER_SAMPLE_DIR}"

    # =========================
    # AMRFinder
    # =========================
    #RUN AMRFinder BLOCK
    {
        if [[ "${RUN_AMR}" == "true" ]]; then
            echo "[INFO] Running AMRFinder"
            run_safe "amrfinder" run_amrfinder "${SAMPLE_ID}" "${ASSEMBLY_FASTA}" "${BAKTA_SAMPLE_DIR}" "${AMRFINDER_OUTFILE}"
            echo "[DONE] AMRFinderPlus: ${SAMPLE_ID}"
        else
            echo "[INFO] AMRFinder skipped: ${SAMPLE_ID}"
        fi
    } >> "${SAMPLE_LOG}" 2>&1

    if [[ "${RUN_AMR}" == "true" ]] && ! check_file_not_empty "${AMRFINDER_OUTFILE}"; then
        status_fail "$SAMPLE_ID" "characterization" "amrfinder_failed"
        continue
    fi

    # =========================
    # PlasmidFinder
    # =========================
    #RUN Plasmidfinder BLOCK
    {
        if [[ "${RUN_PLASMID}" == "true" ]]; then
            echo "[INFO] Running PlasmidFinder"
            run_safe "plasmidfinder" run_plasmidfinder "${SAMPLE_ID}" "${ASSEMBLY_FASTA}" "${PLASMIDFINDER_SAMPLE_DIR}"
            echo "[DONE] PlasmidFinder: ${SAMPLE_ID}"
        else
            echo "[INFO] PlasmidFinder skipped: ${SAMPLE_ID}"
        fi
    } >> "${SAMPLE_LOG}" 2>&1

    if [[ "${RUN_PLASMID}" == "true" ]] && ! check_file_not_empty "${PLASMIDFINDER_OUT}"; then
        status_fail "$SAMPLE_ID" "characterization" "plasmidfinder_failed"
        continue
    fi

    # =========================
    # Mobsuite - MobTyper
    # =========================
    #RUN Mobsuite BLOCK
    {
        if [[ "${RUN_MOBSUITE}" == "true" ]]; then
            echo "[INFO] Running MOB-typer"
            run_mobtyper "${SAMPLE_ID}" "${ASSEMBLY_FASTA}" "${MOBSUITE_OUTFILE}"
            echo "[DONE] MOB-typer: ${SAMPLE_ID}"
        else
            echo "[INFO] MOB-typer skipped: ${SAMPLE_ID}"
        fi
    } >> "${SAMPLE_LOG}" 2>&1

    if [[ "${RUN_MOBSUITE}" == "true" ]] && ! check_file_not_empty "${MOBSUITE_OUTFILE}"; then
        status_fail "$SAMPLE_ID" "characterization" "mobtyper_failed"
        continue
    fi

    echo "Done: ${SAMPLE_ID}"
done < "${VALID_SAMPLES_FILE}"

############################################
# FINAL STATUS PER SAMPLE
############################################

echo "[FINAL] Evaluating characterization status per sample..."

cat "${VALID_SAMPLES_FILE}" | cut -f1 | while read -r sample; do
    # Checking only OK samples
    if grep -qE "^${sample}.*characterization.*STARTED" "${STATUS_FILE}" && \
       ! grep -qE "^${sample}.*characterization.*FAILED" "${STATUS_FILE}"; then
        status_ok "$sample" "characterization"
    fi
done

rm "${VALID_SAMPLES_FILE}"

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
    --amr "${AMRFINDER_SUMMARY_DIR}/amrfinder.AMR.tsv" \
    --virulence "${AMRFINDER_SUMMARY_DIR}/amrfinder.VIRULENCE.tsv" \
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

#Final log info
log_info "Characterization step completed for all samples"
awk '$2=="characterization"' "${STATUS_FILE}"
echo "[INFO] CHARACTERIZATION module completed"
echo "[DONE] Outputs written to: ${ROOT_DIR}/${CHARACTERIZATION_RESULTS_DIR}"