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
# Log Status
# ------------------------------------------------------------

log_status() {
    local sample="$1"
    local step="$2"
    local status="$3"
    local reason="${4:--}"

    # checking
    if [[ -z "${STATUS_FILE:-}" ]]; then
        echo "[ERROR] STATUS_FILE is not defined" >&2
        return 1
    fi

    mkdir -p "$(dirname "$STATUS_FILE")"

    # Timestamp
    local ts
    ts=$(date '+%Y-%m-%d %H:%M:%S')

    # write "status_file"
    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$sample" "$step" "$status" "$reason" "$ts" >> "$STATUS_FILE"
}

# ------------------------------------------------------------
# Status files
# ------------------------------------------------------------

STATUS_STARTED="STARTED"
STATUS_COMPLETED="COMPLETED"
STATUS_FAILED="FAILED"
STATUS_SKIPPED="SKIPPED"

init_status_file() {
    mkdir -p "$(dirname "$STATUS_FILE")"

    local expected=$'sample\tstep\tstatus\treason\ttimestamp'

    if [[ ! -f "$STATUS_FILE" ]]; then
        echo "$expected" > "$STATUS_FILE"
    else
        local header
        header=$(head -1 "$STATUS_FILE")

        if [[ "$header" != "$expected" ]]; then
            echo "[ERROR] STATUS_FILE header mismatch" >&2
            echo "[ERROR] Expected: $expected" >&2
            echo "[ERROR] Found   : $header" >&2
            exit 1
        fi
    fi
}

status_start() {
    log_status "$1" "$2" "$STATUS_STARTED"
}

status_ok() {
    log_status "$1" "$2" "$STATUS_COMPLETED"
}

status_fail() {
    log_status "$1" "$2" "$STATUS_FAILED" "$3"
}

status_skip() {
    log_status "$1" "$2" "$STATUS_SKIPPED" "$3"
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

check_min_reads_nanostat() {
    local nanostat="$1"
    local min_reads="$2"

    if [[ ! -f "$nanostat" ]]; then
        log_warn "Nanostat file missing: $nanostat"
        return 1
    fi

    local reads

    reads=$(grep "Number of reads:" "$nanostat" 2>/dev/null | awk -F: '{print $2}')
    
    if [[ -z "$reads" ]]; then
        log_warn "Could not extract read count from $nanostat"
        return 1
    fi

    reads=$(echo "$reads" | tr -d ' ,' | sed 's/\..*$//')

    if ! [[ "$reads" =~ ^[0-9]+$ ]]; then
        log_warn "Invalid read count: $reads"
        return 1
    fi

    if (( reads < min_reads )); then
        log_warn "Low reads ($reads < $min_reads)"
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

    set +e
    "$@"
    local status=$?
    set -e

    if [[ $status -ne 0 ]]; then
        log_warn "Step failed: ${step_name}"
        return 1
    fi

    return 0
}


# ------------------------------------------------------------
# Skip failed or empty file
# ------------------------------------------------------------

skip_sample() {
    local sample="$1"
    local step="$2"
    local reason="$3"

    log_warn "Skipping ${sample}: $reason"
    status_skip "$sample" "$step" "$reason"

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