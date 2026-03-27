#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${SCRIPT_DIR}"

# =========================
# Defaults
# =========================
SAMPLES=""
TRIMMER="porechop_filtlong"
STEP=""
RESUME="false"
FORCE="false"
RUN_CLEANUP="true"

# =========================
# Usage
# =========================
usage() {
    cat <<EOF
Usage: $0 --samples FILE [options]

Options:
  --samples FILE       Path to sample sheet TSV (required)
  --trimmer NAME       fastplong | porechop_filtlong (default:porechop_filtlong)
  --resume             Skip existing outputs
  --force              Ignore resume and re-run everything
  --no-cleanup         Do not run cleanup step
  -h, --help           Show this help message
EOF
}

# =========================
# Parse arguments
# =========================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --samples)
            [[ $# -lt 2 ]] && { echo "[ERROR] --samples requires a value"; exit 1; }
            SAMPLES="$2"
            shift 2
            ;;
        --trimmer)
            [[ $# -lt 2 ]] && { echo "[ERROR] --trimmer requires a value"; exit 1; }
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
        --no-cleanup)
            RUN_CLEANUP="false"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            usage
            exit 1
            ;;
    esac
done

# =========================
# Post-parse validation
# =========================

# required args
[[ -z "${SAMPLES}" ]] && { echo "[ERROR] --samples is required"; usage; exit 1; }

# file exists
[[ ! -f "${SAMPLES}" ]] && { echo "[ERROR] Sample sheet not found: ${SAMPLES}"; exit 1; }

# validate trimmer
if [[ "${TRIMMER}" != "fastplong" && "${TRIMMER}" != "porechop_filtlong" ]]; then
    echo "[ERROR] Invalid --trimmer value: ${TRIMMER}"
    exit 1
fi

# force overrides resume
if [[ "${FORCE}" == "true" ]]; then
    RESUME="false"
fi

#Running pipeline
echo "=========================================="
echo " Running complete bacterial LR pipeline"
echo " Samples : ${SAMPLES}"
echo " Trimmer : ${TRIMMER}"
echo "=========================================="

COMMON_ARGS=(--samples "${SAMPLES}")
[[ "${RESUME}" == "true" ]] && COMMON_ARGS+=(--resume)
[[ "${FORCE}" == "true" ]] && COMMON_ARGS+=(--force)

bash "${ROOT_DIR}/run_pipeline.sh" --step trimming "${COMMON_ARGS[@]}" --trimmer "${TRIMMER}"
bash "${ROOT_DIR}/run_pipeline.sh" --step assembly "${COMMON_ARGS[@]}" --trimmer "${TRIMMER}"
bash "${ROOT_DIR}/run_pipeline.sh" --step assembly_qc "${COMMON_ARGS[@]}"
bash "${ROOT_DIR}/run_pipeline.sh" --step taxonomy "${COMMON_ARGS[@]}"
bash "${ROOT_DIR}/run_pipeline.sh" --step annotation "${COMMON_ARGS[@]}"
bash "${ROOT_DIR}/run_pipeline.sh" --step characterization "${COMMON_ARGS[@]}"
bash "${ROOT_DIR}/run_pipeline.sh" --step final_summary "${COMMON_ARGS[@]}"
if [[ "${RUN_CLEANUP}" == "true" ]]; then
    bash "${ROOT_DIR}/run_pipeline.sh" --step cleanup "${COMMON_ARGS[@]}" --trimmer "${TRIMMER}"
else
    echo "[INFO] Skipping cleanup step (--no-cleanup enabled)"
fi

echo
echo "=========================================="
echo " Pipeline completed successfully"
echo "=========================================="