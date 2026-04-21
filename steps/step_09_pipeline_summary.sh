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

############################################
# STATUS FILE
############################################

STATUS_FILE="${ROOT_DIR}/results/pipeline_status.tsv"
init_status_file

export STATUS_FILE
export UTILS_FILE

############################################
# DIRECTORIES
############################################

RUN_SUMMARY_LOG_DIR="${ROOT_DIR}/${LOG_DIR}/run_summary"

mkdir -p "${RUN_SUMMARY_LOG_DIR}"

OUT_TSV="${ROOT_DIR}/${RESULTS_DIR}/pipeline_status_report.tsv"
OUT_EXCEL="${ROOT_DIR}/${RESULTS_DIR}/pipeline_status_report.xlsx"
OUT_HTML="${ROOT_DIR}/${RESULTS_DIR}/pipeline_status_report.html"

############################################
# CHECK ENV
############################################

if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_ANALYSIS}"; then
    echo "ERROR: conda environment not found: ${ENV_ANALYSIS}"
    exit 1
fi

############################################
# RUN MAIN
############################################

echo "=================================================="
echo "STEP 09: RUN SUMMARY"
echo "Samples file         : ${SAMPLES}"
echo "Input file           : ${STATUS_FILE}"
echo "Output TSV           : ${OUT_TSV}"
echo "Output HTML          : ${OUT_HTML}"
echo "ENV_ANALYSIS         : ${ENV_ANALYSIS}"
echo "=================================================="

status_start "GLOBAL" "run_summary"

log_info "Generating run summary report"

set +e
conda run --no-capture-output -n "${ENV_ANALYSIS}" python3 \
    "${SCRIPT_PIPELINE_SUMMARY}" \
    --status "${STATUS_FILE}" \
    --html "${OUT_HTML}" \
    --out_tsv "${OUT_TSV}" \
    --out_excel "${OUT_EXCEL}" \
    --pipeline-name "${PIPELINE_NAME}" \
    --pipeline-version "${PIPELINE_VERSION}" \
    --author "${PIPELINE_AUTHOR}" \
    --institution "${PIPELINE_INSTITUTION}" \
    > "${RUN_SUMMARY_LOG_DIR}/run_summary_report.stdout.log" \
    2> "${RUN_SUMMARY_LOG_DIR}/run_summary_report.stderr.log"
status=$?
set -e

if [[ $status -ne 0 ]]; then
    status_fail "GLOBAL" "run_summary" "summary_failed"
    exit 1
fi

#Check output
if ! check_file_not_empty "${OUT_HTML}" || ! check_file_not_empty "${OUT_TSV}" || ! check_file_not_empty "${OUT_EXCEL}"; then
    status_fail "GLOBAL" "run_summary" "empty_output"
    exit 1
fi

#Final log info
log_info "Run summary completed"
status_ok "GLOBAL" "run_summary"

echo "[DONE] Run summary report written to: ${OUT_HTML}"