#!/usr/bin/env bash

# ============================================================
# utils.sh
# Funciones comunes para el pipeline de genómica bacteriana
# ============================================================

set -o pipefail

# ------------------------------------------------------------
# Logging
# ------------------------------------------------------------

log_info() {
    echo "[INFO]  $(date '+%Y-%m-%d %H:%M:%S')  $*"
}

log_warn() {
    echo "[WARN]  $(date '+%Y-%m-%d %H:%M:%S')  $*" >&2
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S')  $*" >&2
}

die() {
    log_error "$*"
    exit 1
}

# ------------------------------------------------------------
# NA handling
# ------------------------------------------------------------

is_na() {
    [[ -z "${1:-}" || "$1" == "NA" || "$1" == "None" || "$1" == "none" ]]
}

normalize_na() {
    local v="${1:-}"
    if is_na "$v"; then
        echo ""
    else
        echo "$v"
    fi
}

# ------------------------------------------------------------
# File validation
# ------------------------------------------------------------

require_file() {
    local f="$1"
    local label="${2:-file}"

    if is_na "$f"; then
        die "Missing ${label}"
    fi

    if [[ ! -f "$f" ]]; then
        die "${label} not found: $f"
    fi
}

optional_file() {
    local f="${1:-}"
    [[ -n "$f" && "$f" != "NA" && -f "$f" ]]
}

require_dir() {
    local d="$1"
    local label="${2:-directory}"

    if is_na "$d"; then
        die "Missing ${label}"
    fi

    if [[ ! -d "$d" ]]; then
        die "${label} not found: $d"
    fi
}

# ------------------------------------------------------------
# Samplesheet validation
# ------------------------------------------------------------

validate_asm_type() {
    local asm="$1"

    if [[ "$asm" != "long" && "$asm" != "hybrid" ]]; then
        die "asm_type must be 'long' or 'hybrid' (got: $asm)"
    fi
}

validate_genome_size() {
    local g="$1"

    if ! [[ "$g" =~ ^[0-9]+$ ]]; then
        die "expected_genome_size must be an integer (got: $g)"
    fi
}

# ------------------------------------------------------------
# File existence check
# ------------------------------------------------------------

check_file_exists() {
    local f="$1"

    if [[ ! -f "$f" ]]; then
        log_warn "File not found: $f"
        return 1
    fi

    return 0
}

# ------------------------------------------------------------
# Samplesheet reader helper
# ------------------------------------------------------------

read_samplesheet() {
    local samples="$1"
    require_file "$samples" "samplesheet"
    awk 'NR>1' "$samples"
}

# ------------------------------------------------------------
# Path utilities
# ------------------------------------------------------------

abs_path() {
    local p="$1"
    if [[ -d "$p" ]]; then
        (cd "$p" && pwd)
    else
        (cd "$(dirname "$p")" && echo "$(pwd)/$(basename "$p")")
    fi
}

# ------------------------------------------------------------
# Empty files check (low input files)
# ------------------------------------------------------------

check_file_not_empty() {
    local file="$1"

    if [[ ! -s "$file" ]]; then
        log_warn "File is empty or missing: $file"
        return 1
    fi

    return 0
}


# ------------------------------------------------------------
# check FASTQ files (low number of reads)
# ------------------------------------------------------------

check_min_reads() {
    local fastq="$1"
    local min_reads="$2"

    local count
    count=$(grep -c "^+$" "$fastq" 2>/dev/null || echo 0)

    if [[ "$count" -lt "$min_reads" ]]; then
        log_warn "Low read count: $count (< $min_reads)"
        return 1
    fi

    return 0
}


# ------------------------------------------------------------
# Safe run, run the pipeline everytime
# ------------------------------------------------------------

run_safe() {
    local step_name="$1"
    shift

    if ! "$@"; then
        log_warn "Step failed: ${step_name}"
        return 1
    fi

    return 0
}


# ------------------------------------------------------------
# Skip failed or empty file
# ------------------------------------------------------------

skip_sample() {
    local reason="$1"

    log_warn "Skipping sample: $reason"
    return 1
}

# ==========================================
# Resume / skip utilities
# ==========================================

should_skip_sample() {
    local output_file="$1"

    if [[ "${FORCE:-false}" == "true" ]]; then
        return 1
    fi

    if [[ "${RESUME:-false}" == "true" && -s "${output_file}" ]]; then
        log_info "Skipping ${SAMPLE_ID}: output already exists -> ${output_file}"
        return 0
    fi

    return 1
}

should_skip_multiple_outputs() {
    # skip only if ALL outputs exist
    local files=("$@")

    if [[ "${FORCE:-false}" == "true" ]]; then
        return 1
    fi

    if [[ "${RESUME:-false}" != "true" ]]; then
        return 1
    fi

    for f in "${files[@]}"; do
        [[ -s "$f" ]] || return 1
    done

    log_info "Skipping ${SAMPLE_ID}: all outputs already exist"
    return 0
}

should_skip_global() {
    local output_file="$1"

    if [[ "${FORCE:-false}" == "true" ]]; then
        return 1
    fi

    if [[ "${RESUME:-false}" == "true" && -s "${output_file}" ]]; then
        log_info "Skipping step: output already exists -> ${output_file}"
        return 0
    fi

    return 1
}

# ------------------------------------------------------------
# CPU detection helper
# ------------------------------------------------------------

detect_cpus() {
    if command -v nproc >/dev/null 2>&1; then
        nproc
    else
        sysctl -n hw.ncpu 2>/dev/null || echo 1
    fi
}