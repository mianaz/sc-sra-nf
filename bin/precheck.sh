#!/usr/bin/env bash
set -euo pipefail
DEBUG=${1:-0}
USE_CONTAINER=${USE_CONTAINER:-false}
CONTAINER_ENGINE=${CONTAINER_ENGINE:-docker}

# Enhanced logging function
log_info() {
    echo "INFO: $*" >&2
}

log_warn() {
    echo "WARNING: $*" >&2
}

log_error() {
    echo "ERROR: $*" >&2
}

log_debug() {
    if [[ $DEBUG -eq 1 ]]; then
        echo "DEBUG: $*" >&2
    fi
}

# Enhanced version checking
check_bin() {
    local name=$1
    local required=${2:-true}

    if ! command -v "$name" >/dev/null 2>&1; then
        if [[ "$required" == "true" ]]; then
            log_error "Required software '$name' not found in PATH"
            return 1
        else
            log_warn "Optional software '$name' not found in PATH"
            return 0
        fi
    fi

    local location=$(command -v "$name")
    log_debug "Found $name: $location"

    # Try to get version info
    local version=""
    case "$name" in
        prefetch|fasterq-dump|fastq-dump)
            version=$("$name" --version 2>/dev/null | head -1 || echo "version unknown")
            ;;
        cellranger)
            version=$("$name" --version 2>&1 | head -1 || echo "version unknown")
            ;;
        python3)
            version=$("$name" --version 2>/dev/null || echo "version unknown")
            ;;
        docker)
            version=$("$name" --version 2>/dev/null | head -1 || echo "version unknown")
            ;;
        singularity)
            version=$("$name" --version 2>/dev/null || echo "version unknown")
            ;;
        *)
            version=$("$name" --version 2>/dev/null | head -1 || echo "version unknown")
            ;;
    esac

    log_debug "$name version: $version"
    return 0
}

# Platform detection
detect_platform() {
    local platform="unknown"
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        platform="linux"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        platform="macos"
    elif [[ "$OSTYPE" == "cygwin" ]] || [[ "$OSTYPE" == "msys" ]]; then
        platform="windows"
    fi
    echo "$platform"
}

# Memory detection
check_memory() {
    local available_gb=0

    if command -v free >/dev/null 2>&1; then
        # Linux
        available_gb=$(free -g | awk '/^Mem:/{print $2}')
    elif command -v vm_stat >/dev/null 2>&1; then
        # macOS
        local pages=$(vm_stat | grep "free\|inactive" | awk '{sum += $3} END {print sum}')
        available_gb=$((pages * 4096 / 1024 / 1024 / 1024))
    fi

    log_debug "Available memory: ~${available_gb}GB"

    if [[ $available_gb -lt 16 ]]; then
        log_warn "Low memory detected (~${available_gb}GB). Cell Ranger requires at least 16GB"
    fi
}

# Storage check
check_storage() {
    local current_dir=$(pwd)
    local available_gb=0

    if command -v df >/dev/null 2>&1; then
        # macOS vs Linux df compatibility
        if [[ "$OSTYPE" == "darwin"* ]]; then
            available_gb=$(df -g "$current_dir" | tail -1 | awk '{print $4}')
        else
            available_gb=$(df -BG "$current_dir" | tail -1 | awk '{print $4}' | sed 's/G//')
        fi

        log_debug "Available storage: ~${available_gb}GB"

        if [[ $available_gb -lt 100 ]]; then
            log_warn "Low storage space (~${available_gb}GB). Single-cell data can be large"
        fi
    fi
}

# Container engine specific checks
check_container_engine() {
    local engine=$1

    case "$engine" in
        docker)
            if ! docker info >/dev/null 2>&1; then
                log_error "Docker daemon not running or not accessible"
                return 1
            fi
            log_debug "Docker daemon is running"

            # Check if we can pull images
            if [[ $DEBUG -eq 1 ]]; then
                log_debug "Testing Docker image access..."
                if docker images etycksen/cellranger >/dev/null 2>&1; then
                    log_debug "Cell Ranger container image available locally"
                fi
            fi
            ;;
        singularity)
            # Basic singularity check
            if ! singularity --version >/dev/null 2>&1; then
                log_error "Singularity not properly installed"
                return 1
            fi
            log_debug "Singularity is available"
            ;;
    esac

    return 0
}

# Main precheck logic
main() {
    local platform=$(detect_platform)
    log_info "Running precheck on $platform platform"

    # Check system resources
    check_memory
    check_storage

    # Core SRA tools
    log_info "Checking SRA toolkit..."
    if check_bin prefetch false; then
        local version=$(prefetch --version 2>&1 | head -1)
        log_info "  prefetch: $version"
    fi

    # FASTQ conversion tools (at least one required)
    local fastq_tool_found=false
    if check_bin fasterq-dump false; then
        fastq_tool_found=true
        local version=$(fasterq-dump --version 2>&1 | head -1)
        log_info "  fasterq-dump: $version"
    elif check_bin fastq-dump false; then
        fastq_tool_found=true
        local version=$(fastq-dump --version 2>&1 | head -1)
        log_info "  fastq-dump: $version"
        log_warn "Using older fastq-dump. Consider upgrading to fasterq-dump for better performance"
    fi

    if [[ "$fastq_tool_found" == "false" ]]; then
        log_error "No FASTQ conversion tool found. Install SRA toolkit (fasterq-dump or fastq-dump)"
        exit 1
    fi

    # Python
    log_info "Checking Python..."
    if ! check_bin python3 true; then
        exit 1
    fi
    local py_version=$(python3 --version 2>&1)
    log_info "  $py_version"

    # Check quantifiers based on container mode
    QUANTIFIER="${QUANTIFIER:-auto}"

    log_info ""
    log_info "=== Quantifier Configuration ==="
    log_info "Quantifier mode: ${QUANTIFIER}"

    case "${USE_CONTAINER}" in
        true)
            # Force container mode - check container engine
            log_info "Container mode: FORCE (use_container=true)"
            log_info "Checking ${CONTAINER_ENGINE}..."
            if ! check_bin "${CONTAINER_ENGINE}" true; then
                exit 1
            fi
            if ! check_container_engine "${CONTAINER_ENGINE}"; then
                exit 1
            fi

            # Show which quantifiers will use containers
            case "$QUANTIFIER" in
                cellranger|auto)
                    log_info ""
                    log_info "Quantifier: CellRanger (forced container)"
                    log_info "  Container: cumulusprod/cellranger:9.0.1"
                    log_info "  Engine: ${CONTAINER_ENGINE}"
                    ;;
                starsolo)
                    log_info ""
                    log_info "Quantifier: STARsolo (forced container)"
                    log_info "  Container: quay.io/biocontainers/star:2.7.11a--h0033a41_0"
                    log_info "  Engine: ${CONTAINER_ENGINE}"
                    ;;
                alevin)
                    log_info ""
                    log_info "Quantifier: Salmon/Alevin (forced container)"
                    log_info "  Container: quay.io/biocontainers/salmon:1.10.0--h7e5ed60_0"
                    log_info "  Engine: ${CONTAINER_ENGINE}"
                    ;;
                kallisto)
                    log_info ""
                    log_info "Quantifier: Kallisto/BUStools (forced container)"
                    log_info "  Container: quay.io/biocontainers/kallisto:0.48.0--h15996b6_2"
                    log_info "  Engine: ${CONTAINER_ENGINE}"
                    ;;
            esac
            ;;

        false)
            # Never use containers - all software must be installed
            log_info "Container mode: DISABLED (use_container=false)"
            log_info "Checking native installations..."

            case "$QUANTIFIER" in
                cellranger|auto)
                    log_info ""
                    log_info "Quantifier: CellRanger (native)"
                    if ! check_bin cellranger true; then
                        log_error "Cell Ranger not found and containers disabled"
                        exit 1
                    fi
                    local cr_version=$(cellranger --version 2>&1 | head -1)
                    log_info "  Version: $cr_version"
                    log_info "  Path: $(which cellranger)"
                    ;;
                starsolo)
                    log_info ""
                    log_info "Quantifier: STARsolo (native)"
                    if ! check_bin STAR true; then
                        log_error "STAR not found and containers disabled"
                        exit 1
                    fi
                    local star_version=$(STAR --version 2>&1 | head -1)
                    log_info "  Version: $star_version"
                    log_info "  Path: $(which STAR)"
                    ;;
                alevin)
                    log_info ""
                    log_info "Quantifier: Salmon/Alevin (native)"
                    if ! check_bin salmon true; then
                        log_error "Salmon not found and containers disabled"
                        exit 1
                    fi
                    local salmon_version=$(salmon --version 2>&1 | head -1)
                    log_info "  Version: $salmon_version"
                    log_info "  Path: $(which salmon)"
                    ;;
                kallisto)
                    log_info ""
                    log_info "Quantifier: Kallisto/BUStools (native)"
                    if ! check_bin kallisto true || ! check_bin bustools true; then
                        log_error "Kallisto/BUStools not found and containers disabled"
                        exit 1
                    fi
                    local kallisto_version=$(kallisto version 2>&1 | head -1)
                    local bustools_version=$(bustools version 2>&1 | head -1)
                    log_info "  Kallisto version: $kallisto_version"
                    log_info "  BUStools version: $bustools_version"
                    log_info "  Paths: $(which kallisto), $(which bustools)"
                    ;;
            esac
            ;;

        auto)
            # Smart mode - check what's available, use containers as fallback
            log_info "Container mode: AUTO (use_container=auto)"
            log_info "Detecting available software..."

            # Check container engine availability first
            local container_available=false
            if check_bin "${CONTAINER_ENGINE}" false; then
                if check_container_engine "${CONTAINER_ENGINE}" 2>/dev/null; then
                    container_available=true
                    log_info "Container engine (${CONTAINER_ENGINE}) available as fallback"
                fi
            fi

            # Check software based on quantifier
            case "$QUANTIFIER" in
                cellranger|auto)
                    log_info ""
                    if check_bin cellranger false; then
                        local version=$(cellranger --version 2>&1 | head -1)
                        log_info "Quantifier: CellRanger (native - auto-detected)"
                        log_info "  ✓ Found in PATH"
                        log_info "  Version: $version"
                        log_info "  Path: $(which cellranger)"
                        echo "CELLRANGER_NATIVE=true" >> precheck_results.txt
                    elif [[ "$container_available" == "true" ]]; then
                        log_info "Quantifier: CellRanger (container - fallback)"
                        log_info "  ⚠ Not found in PATH, using container"
                        log_info "  Container: cumulusprod/cellranger:9.0.1"
                        log_info "  Engine: ${CONTAINER_ENGINE}"
                        echo "CELLRANGER_CONTAINER=true" >> precheck_results.txt
                    else
                        log_error "CellRanger not found and no container engine available"
                        exit 1
                    fi
                    ;;
                starsolo)
                    log_info ""
                    if check_bin STAR false; then
                        local version=$(STAR --version 2>&1 | head -1)
                        log_info "Quantifier: STARsolo (native - auto-detected)"
                        log_info "  ✓ Found in PATH"
                        log_info "  Version: $version"
                        log_info "  Path: $(which STAR)"
                        echo "STAR_NATIVE=true" >> precheck_results.txt
                    elif [[ "$container_available" == "true" ]]; then
                        log_info "Quantifier: STARsolo (container - fallback)"
                        log_info "  ⚠ Not found in PATH, using container"
                        log_info "  Container: quay.io/biocontainers/star:2.7.11a--h0033a41_0"
                        log_info "  Engine: ${CONTAINER_ENGINE}"
                        echo "STAR_CONTAINER=true" >> precheck_results.txt
                    else
                        log_error "STAR not found and no container engine available"
                        exit 1
                    fi
                    ;;
                alevin)
                    log_info ""
                    if check_bin salmon false; then
                        local version=$(salmon --version 2>&1 | head -1)
                        log_info "Quantifier: Salmon/Alevin (native - auto-detected)"
                        log_info "  ✓ Found in PATH"
                        log_info "  Version: $version"
                        log_info "  Path: $(which salmon)"
                        echo "SALMON_NATIVE=true" >> precheck_results.txt
                    elif [[ "$container_available" == "true" ]]; then
                        log_info "Quantifier: Salmon/Alevin (container - fallback)"
                        log_info "  ⚠ Not found in PATH, using container"
                        log_info "  Container: quay.io/biocontainers/salmon:1.10.0--h7e5ed60_0"
                        log_info "  Engine: ${CONTAINER_ENGINE}"
                        echo "SALMON_CONTAINER=true" >> precheck_results.txt
                    else
                        log_error "Salmon not found and no container engine available"
                        exit 1
                    fi
                    ;;
                kallisto)
                    log_info ""
                    if check_bin kallisto false && check_bin bustools false; then
                        local kallisto_version=$(kallisto version 2>&1 | head -1)
                        local bustools_version=$(bustools version 2>&1 | head -1)
                        log_info "Quantifier: Kallisto/BUStools (native - auto-detected)"
                        log_info "  ✓ Found in PATH"
                        log_info "  Kallisto version: $kallisto_version"
                        log_info "  BUStools version: $bustools_version"
                        log_info "  Paths: $(which kallisto), $(which bustools)"
                        echo "KALLISTO_NATIVE=true" >> precheck_results.txt
                    elif [[ "$container_available" == "true" ]]; then
                        log_info "Quantifier: Kallisto/BUStools (container - fallback)"
                        log_info "  ⚠ Not found in PATH, using container"
                        log_info "  Container: quay.io/biocontainers/kallisto:0.48.0--h15996b6_2"
                        log_info "  Engine: ${CONTAINER_ENGINE}"
                        echo "KALLISTO_CONTAINER=true" >> precheck_results.txt
                    else
                        log_error "Kallisto/BUStools not found and no container engine available"
                        exit 1
                    fi
                    ;;
            esac
            ;;

        *)
            log_error "Invalid use_container value: ${USE_CONTAINER}. Must be 'true', 'false', or 'auto'"
            exit 1
            ;;
    esac

    # Platform-specific warnings
    case "$platform" in
        macos)
            log_info "macOS detected - Docker recommended for Cell Ranger"
            if [[ "${USE_CONTAINER}" != "true" ]]; then
                log_warn "Consider using --profile macos for containerized Cell Ranger"
            fi
            ;;
        linux)
            log_info "Linux detected - native Cell Ranger recommended"
            ;;
    esac

    log_info "Precheck completed successfully"
}

main "$@"
