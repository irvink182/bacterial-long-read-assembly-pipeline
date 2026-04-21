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
DRY_RUN="false"

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
        --dry-run)
            DRY_RUN="true"
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

[[ -z "${SAMPLES}" ]] && { echo "[ERROR] missing --samples"; exit 1; }

require_file "${SAMPLES}" "samples file"

export RESUME
export FORCE

############################################
# STATUS FILE
############################################

STATUS_FILE="${ROOT_DIR}/results/pipeline_status.tsv"
init_status_file

export STATUS_FILE
export UTILS_FILE

# =========================
# Info
# =========================
log_info "Starting cleanup step"
[[ "${RESUME}" == "true" ]] && log_info "Resume mode enabled"


############################################
# DIRECTORIES
############################################

# Reads
FASTPLONG_QC_DIR="${RESULTS_DIR}/reads/qc/fastplong"
PORECHOP_DIR="${RESULTS_DIR}/reads/porechop_reads"
PORECHOP_QC_DIR="${RESULTS_DIR}/reads/qc/porechop"
FILTLONG_QC_DIR="${RESULTS_DIR}/reads/qc/filtlong"
#clean read directory, only in aggresive mode?
#CLEAN_DIR="${RESULTS_DIR}/reads/clean_reads"

# Assemblies
FLYE_DIR="${RESULTS_DIR}/assemblies/flye_asm"
MEDAKA_DIR="${RESULTS_DIR}/assemblies/medaka_asm"
COVERAGE_DIR="${RESULTS_DIR}/assemblies/coverage_asm"

############################################
# FUNCTIONS
############################################

remove_path() {
    local path="$1"

    # if varible is empty, abort for safety
    [[ -z "${path}" || "${path}" == "/" ]] && return 0
    [[ ! -e "$path" ]] && return 0

    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "[DRY-RUN] Would remove: $path"
    else
        echo "[INFO] Removing: $path"
        rm -rf "$path"
    fi
}

############################################
# CLEANUP FUNCTIONS
############################################

cleanup_reads_by_trimmer() {
    case "${TRIMMER}" in
        fastplong)
            log_info "Cleanup for fastplong"

            remove_path "${FASTPLONG_QC_DIR}"
            # opcional:
            # remove_path "${CLEAN_DIR}"

            ;;

        porechop_filtlong)
            log_info "Cleanup for porechop_filtlong"

            remove_path "${PORECHOP_DIR}"
            remove_path "${PORECHOP_QC_DIR}"
            remove_path "${FILTLONG_QC_DIR}"
            # opcional:
            # remove_path "${CLEAN_DIR}"

            ;;

        *)
            log_warn "Unknown trimmer: ${TRIMMER}, skipping read cleanup"
            ;;
    esac
}

cleanup_assemblies() {
    log_info "Cleaning assembly intermediate files"

    remove_path "${FLYE_DIR}"
    remove_path "${MEDAKA_DIR}"
}

cleanup_coverage_assemblies() {
    log_info "Cleaning assembly coverage intermediate files"

    find "${COVERAGE_DIR}" -type f \
        \( -name "*.mosdepth.*" -o -name "*.bed.gz*" \) 2>/dev/null | while read -r file; do

        if [[ "${DRY_RUN}" == "true" ]]; then
            echo "[DRY-RUN] Would remove: ${file}"
        else
            echo "[INFO] Removing: ${file}"
            rm -f "${file}"
        fi

    done     

}

cleanup_plasmidfinder_tmp() {
    log_info "Cleaning plasmidfinder tmp directories"

    find "${RESULTS_DIR}/characterization/plasmidfinder" \
        -type d -name "tmp" 2>/dev/null | while read -r tmp_dir; do

        if [[ "${DRY_RUN}" == "true" ]]; then
            echo "[DRY-RUN] Would remove: ${tmp_dir}"
        else
            echo "[INFO] Removing: ${tmp_dir}"
            rm -rf "${tmp_dir}"
        fi

    done
}

############################################
# RUN CLEANUP
############################################

echo "=================================================="
echo "STEP 08: CLEANUP"
echo "=================================================="

[[ "${FORCE}" == "true" ]] && log_info "Force mode enabled"
[[ "${DRY_RUN}" == "true" ]] && log_info "Dry-run mode enabled"

#Status global
status_start "GLOBAL" "cleanup"

# Load status file to know which samples were included in this run
# We use 'trimming' step as a reference because is the first step in this pipeline
cat "${STATUS_FILE}" | grep "trimming" | cut -f1 | sort -u | while read -r SID; do
    
    log_info "Cleaning up data for sample: ${SID}"
    status_start "${SID}" "cleanup"

    # Removing directories
    cleanup_failed=0
    
    # Cleanup
    cleanup_reads_by_trimmer "${SID}" || cleanup_failed=1
    cleanup_plasmidfinder_tmp "${SID}" || cleanup_failed=1

    if [[ $cleanup_failed -eq 0 ]]; then
        status_ok "${SID}" "cleanup"
    else
        status_fail "${SID}" "cleanup" "partial_cleanup_error"
    fi
done

# Delting shared directories for all samples
cleanup_assemblies 
cleanup_coverage_assemblies


#Final log info
log_info "cleanup step completed"
status_ok "GLOBAL" "cleanup"
echo "[INFO] Cleanup completed"