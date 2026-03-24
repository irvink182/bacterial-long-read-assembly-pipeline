#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${SCRIPT_DIR}"

SAMPLES="samplesheet.tsv"
TRIMMER="porechop_filtlong"

usage() {
    cat <<EOF
Usage: $0 [options]

Options:
  --samples FILE       Path to sample sheet TSV
  --trimmer NAME       Trimming method: fastplong or porechop_filtlong
  -h, --help           Show this help message
EOF
}

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

if [[ ! -f "${SAMPLES}" ]]; then
    echo "[ERROR] Sample sheet not found: ${SAMPLES}" >&2
    exit 1
fi

if [[ "${TRIMMER}" != "fastplong" && "${TRIMMER}" != "porechop_filtlong" ]]; then
    echo "[ERROR] Invalid --trimmer value: ${TRIMMER}" >&2
    exit 1
fi

echo "=========================================="
echo " Running complete bacterial LR pipeline"
echo " Samples : ${SAMPLES}"
echo " Trimmer : ${TRIMMER}"
echo "=========================================="

bash "${ROOT_DIR}/run_pipeline.sh" --step trimming --samples "${SAMPLES}" --trimmer "${TRIMMER}"
bash "${ROOT_DIR}/run_pipeline.sh" --step assembly --samples "${SAMPLES}" --trimmer "${TRIMMER}"
bash "${ROOT_DIR}/run_pipeline.sh" --step assembly_qc --samples "${SAMPLES}"
bash "${ROOT_DIR}/run_pipeline.sh" --step taxonomy --samples "${SAMPLES}"
bash "${ROOT_DIR}/run_pipeline.sh" --step annotation --samples "${SAMPLES}"
bash "${ROOT_DIR}/run_pipeline.sh" --step characterization --samples "${SAMPLES}"
bash "${ROOT_DIR}/run_pipeline.sh" --step final_summary --samples "${SAMPLES}"

echo
echo "=========================================="
echo " Pipeline completed successfully"
echo "=========================================="