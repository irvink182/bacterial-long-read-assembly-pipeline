#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${SCRIPT_DIR}"

STEPS_DIR="${ROOT_DIR}/steps"
SCRIPTS_DIR="${ROOT_DIR}/scripts"
CONFIG_DIR="${ROOT_DIR}/config"

CONFIG_FILE="${ROOT_DIR}/config/config.sh"
DB_FILE="${ROOT_DIR}/config/databases.config"
UTILS_FILE="${ROOT_DIR}/scripts/utils.sh"

[[ -f "${CONFIG_FILE}" ]] || { echo "[ERROR] Config file not found: ${CONFIG_FILE}" >&2; exit 1; }
[[ -f "${DB_FILE}" ]] || { echo "[ERROR] Database config file not found: ${DB_FILE}" >&2; exit 1; }

# =========================
# Defaults
# =========================

STEP=""
SAMPLES=""
TRIMMER_OVERRIDE=""
FILTER_CONTIGS="true"
RUN_DNAAPLER="true"
GENERATE_CONTIG_STATS="true"
MIN_COV_FACTOR="0.3"
RESUME="false"
FORCE="false"

usage() {
    echo "Usage:"
    echo "  bash run_pipeline.sh --step <step_name> --samples <samples.tsv> [options]"
    echo
    echo "Steps:"
    echo "  trimming"
    echo "  assembly"
    echo "  assembly_qc"
    echo "  taxonomy"
    echo "  annotation"
    echo "  characterization"
    echo "  final_summary"
    echo "  cleanup"
    echo "  all"
    echo
    echo "Optional overrides:"
    echo "  --trimmer fastplong|porechop_filtlong"
    echo "  --no-filter-contigs     Disable contig filtering"
    echo "  --no-dnaapler           Disable dnaapler reorientation"
    echo "  --no-contig-stats       Disable contig stats generation"
    echo "  --min-cov-factor FLOAT  Coverage threshold factor (default: 0.2)"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --step)
            STEP="$2"
            shift 2
            ;;
        --samples)
            SAMPLES="$2"
            shift 2
            ;;
        --trimmer)
            TRIMMER_OVERRIDE="$2"
            shift 2
            ;;
        --no-filter-contigs)
            FILTER_CONTIGS="false"
            shift
            ;;
        --no-dnaapler)
            RUN_DNAAPLER="false"
            shift
            ;;
        --no-contig-stats)
            GENERATE_CONTIG_STATS="false"
            shift
            ;;
        --min-cov-factor)
            MIN_COV_FACTOR="$2"
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
        *)
            echo "Unknown argument: $1"
            usage
            exit 1
            ;;
    esac
done

source "${CONFIG_FILE}"
source "${DB_FILE}"
export FILTER_CONTIGS
export RUN_DNAAPLER
export GENERATE_CONTIG_STATS
export MIN_COV_FACTOR
export RESUME
export FORCE

if [[ "${FORCE}" == "true" ]]; then
    RESUME="false"
fi

echo "[DEBUG run_pipeline] RESUME=${RESUME}"

[[ -z "${STEP}" ]] && usage
[[ -z "${SAMPLES}" ]] && usage
[[ ! -f "${SAMPLES}" ]] && { echo "Samples file not found: ${SAMPLES}"; exit 1; }

# Apply overrides
if [[ -n "${TRIMMER_OVERRIDE}" ]]; then
    export TRIMMER="${TRIMMER_OVERRIDE}"
fi

# Validate trimmer
if [[ "${TRIMMER}" != "fastplong" && "${TRIMMER}" != "porechop_filtlong" ]]; then
    echo "ERROR: invalid --trimmer value: ${TRIMMER}"
    echo "Allowed values: fastplong | porechop_filtlong"
    exit 1
fi

echo "=============================================="
echo "PIPELINE CONFIGURATION"
echo "STEP         : ${STEP}"
echo "SAMPLES      : ${SAMPLES}"
echo "TRIMMER      : ${TRIMMER:-not_set}"
echo "THREADS      : ${THREADS:-not_set}"
echo "ENV_TRIMMING : ${ENV_TRIMMING:-not_set}"
echo "ENV_ASSEMBLY : ${ENV_ASSEMBLY:-not_set}"
echo "=============================================="

run_step() {
    local step_name="$1"
    local script_path="$2"

    if [[ ! -f "${script_path}" ]]; then
        echo "ERROR: script not found for step '${step_name}': ${script_path}"
        exit 1
    fi

    echo "=================================================="
    echo "Running step: ${step_name}"
    echo "Samples file: ${SAMPLES}"
    echo "Script       : ${script_path}"
    echo "=================================================="

    STEP_ARGS=(--samples "${SAMPLES}")
    [[ "${RESUME}" == "true" ]] && STEP_ARGS+=(--resume)

    if [[ -n "${TRIMMER:-}" ]]; then
        STEP_ARGS+=(--trimmer "${TRIMMER}")
    fi

    bash "${script_path}" "${STEP_ARGS[@]}"
}

case "${STEP}" in
    trimming)
        run_step "trimming" "${STEPS_DIR}/step_01_trimming.sh"
        ;;
    assembly)
        run_step "assembly" "${STEPS_DIR}/step_02_assembly.sh"
        ;;
    assembly_qc)
        run_step "assembly_qc" "${STEPS_DIR}/step_03_assembly_qc.sh"
        ;;
    taxonomy)
        run_step "taxonomy" "${STEPS_DIR}/step_04_taxonomy.sh"
        ;;
    annotation)
        run_step "annotation" "${STEPS_DIR}/step_05_annotation.sh"
        ;;
    characterization)
        run_step "characterization" "${STEPS_DIR}/step_06_characterization.sh"
        ;;
    final_summary)
        run_step "final_summary" "${STEPS_DIR}/step_07_final_summary.sh"
        ;;
    cleanup)
        run_step "cleanup" "${STEPS_DIR}/step_08_cleanup.sh"
        ;;
    all)
        run_step "trimming" "${STEPS_DIR}/step_01_trimming.sh"
        run_step "assembly" "${STEPS_DIR}/step_02_assembly.sh"
        run_step "assembly_qc" "${STEPS_DIR}/step_03_assembly_qc.sh"
        run_step "taxonomy" "${STEPS_DIR}/step_04_taxonomy.sh"
        run_step "annotation" "${STEPS_DIR}/step_05_annotation.sh"
        run_step "characterization" "${STEPS_DIR}/step_06_characterization.sh"
        run_step "final_summary" "${STEPS_DIR}/step_07_final_summary.sh"
        run_step "cleanup" "${STEPS_DIR}/step_08_cleanup.sh"
        ;;
    *)
        echo "Invalid step: ${STEP}"
        usage
        ;;
esac