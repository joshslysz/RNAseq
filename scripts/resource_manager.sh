#!/bin/bash

# Resource Management Module
# Author: Joshua Slysz, PhD - Dalhousie University
# Memory and system resource monitoring for RNA-seq pipeline
# Provides dynamic resource allocation and optimization

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Logging functions
log() { echo -e "${BLUE}[RESOURCE]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }

# Get total system memory in GB
get_total_memory() {
    if [[ "$(uname)" == "Darwin" ]]; then
        local mem_bytes
        mem_bytes=$(sysctl -n hw.memsize 2>/dev/null)
        if [[ -z "$mem_bytes" || "$mem_bytes" -eq 0 ]]; then
            # Alternative method for older macOS versions
            mem_bytes=$(system_profiler SPHardwareDataType | grep "Memory:" | awk '{print $2}' | sed 's/GB//' 2>/dev/null)
            if [[ -n "$mem_bytes" ]]; then
                echo "$mem_bytes"
            else
                # Final fallback - estimate from vm_stat
                echo "8"  # Conservative estimate
            fi
        else
            echo $((mem_bytes / 1024 / 1024 / 1024))
        fi
    else
        local mem_kb
        mem_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
        echo $((mem_kb / 1024 / 1024))
    fi
}

# Get available memory in GB
get_available_memory() {
    if [[ "$(uname)" == "Darwin" ]]; then
        # Parse vm_stat output for available memory
        local page_size free_pages inactive_pages
        
        # Get page size - handle different vm_stat output formats
        page_size=$(vm_stat | head -1 | grep -o '[0-9]*' | head -1)
        if [[ -z "$page_size" ]]; then
            # Fallback: standard macOS page size is 4096 bytes
            page_size=4096
        fi
        
        # Get free and inactive pages
        free_pages=$(vm_stat | grep "Pages free" | grep -o '[0-9]*' | head -1)
        inactive_pages=$(vm_stat | grep "Pages inactive" | grep -o '[0-9]*' | head -1)
        
        # Handle case where values might be empty
        free_pages=${free_pages:-0}
        inactive_pages=${inactive_pages:-0}
        
        # Calculate available memory in GB
        echo $(( (free_pages + inactive_pages) * page_size / 1024 / 1024 / 1024 ))
    else
        local mem_kb
        mem_kb=$(grep MemAvailable /proc/meminfo | awk '{print $2}')
        echo $((mem_kb / 1024 / 1024))
    fi
}

# Get memory usage in GB
get_memory_usage() {
    if [[ "$(uname)" == "Darwin" ]]; then
        local total_gb available_gb
        total_gb=$(get_total_memory)
        available_gb=$(get_available_memory)
        echo $((total_gb - available_gb))
    else
        local total_kb used_kb
        total_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
        used_kb=$(( total_kb - $(grep MemAvailable /proc/meminfo | awk '{print $2}') ))
        echo $((used_kb / 1024 / 1024))
    fi
}

# Get memory usage percentage
get_memory_percentage() {
    if [[ "$(uname)" == "Darwin" ]]; then
        local total_gb used_gb
        total_gb=$(get_total_memory)
        used_gb=$(get_memory_usage)
        echo $(( used_gb * 100 / total_gb ))
    else
        local total_kb used_kb
        total_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
        used_kb=$(( total_kb - $(grep MemAvailable /proc/meminfo | awk '{print $2}') ))
        echo $(( used_kb * 100 / total_kb ))
    fi
}

# Get CPU count
get_cpu_count() {
    if [[ "$(uname)" == "Darwin" ]]; then
        sysctl -n hw.ncpu
    else
        nproc
    fi
}

# Get disk usage for a directory in GB
get_disk_usage() {
    local dir=$1
    du -s "$dir" 2>/dev/null | awk '{print int($1/1024/1024)}'
}

# Get available disk space in GB
get_available_disk() {
    local dir=$1
    df "$dir" | tail -1 | awk '{print int($4/1024/1024)}'
}

# Check if enough memory is available
check_memory_requirement() {
    local required_gb=$1
    local available_gb
    available_gb=$(get_available_memory)
    
    if [[ $available_gb -lt $required_gb ]]; then
        return 1
    else
        return 0
    fi
}

# Check if enough disk space is available
check_disk_requirement() {
    local dir=$1
    local required_gb=$2
    local available_gb
    available_gb=$(get_available_disk "$dir")
    
    if [[ $available_gb -lt $required_gb ]]; then
        return 1
    else
        return 0
    fi
}

# Select optimal RNA-seq aligner based on available memory
select_aligner() {
    local available_mem
    available_mem=$(get_available_memory)
    
    # Check for RNA-seq aligners first
    if command -v hisat2 >/dev/null 2>&1 && [[ $available_mem -ge 4 ]]; then
        echo "hisat2"
    elif command -v STAR >/dev/null 2>&1 && [[ $available_mem -ge 30 ]]; then
        echo "star"
    elif command -v hisat2 >/dev/null 2>&1; then
        echo "hisat2"  # HISAT2 works even with limited memory
    elif [[ $available_mem -ge 2 ]]; then
        echo "minimap2"  # Fallback for very limited systems
    else
        echo "hisat2"  # Default recommendation
    fi
}

# Get optimal thread count based on system resources
get_optimal_threads() {
    local cpu_count mem_gb operation
    cpu_count=$(get_cpu_count)
    mem_gb=$(get_available_memory)
    operation=${1:-"default"}
    
    case "$operation" in
        "alignment")
            # Alignment is memory-intensive, limit threads based on memory
            if [[ $mem_gb -ge 16 ]]; then
                echo $cpu_count
            elif [[ $mem_gb -ge 8 ]]; then
                echo $((cpu_count / 2))
            else
                echo 1
            fi
            ;;
        "compression"|"fastqc")
            # These operations benefit from more threads
            echo $cpu_count
            ;;
        *)
            # Conservative default
            echo $((cpu_count / 2))
            ;;
    esac
}

# Monitor memory usage during a process
monitor_memory() {
    local pid=$1
    local log_file=$2
    local interval=${3:-5}
    
    echo "timestamp,memory_gb,memory_percent" > "$log_file"
    
    while kill -0 "$pid" 2>/dev/null; do
        local timestamp memory_gb memory_percent
        timestamp=$(date '+%Y-%m-%d %H:%M:%S')
        memory_gb=$(get_memory_usage)
        memory_percent=$(get_memory_percentage)
        
        echo "$timestamp,$memory_gb,$memory_percent" >> "$log_file"
        
        # Alert if memory usage is too high
        if [[ $memory_percent -gt 90 ]]; then
            warning "High memory usage: ${memory_percent}%"
        fi
        
        sleep "$interval"
    done
}

# Clean up temporary files
cleanup_temp_files() {
    local temp_dir=$1
    local keep_files=${2:-false}
    
    if [[ "$keep_files" == "false" && -d "$temp_dir" ]]; then
        log "Cleaning up temporary files in $temp_dir"
        rm -rf "$temp_dir"/*
        success "Temporary files cleaned"
    fi
}

# Create memory-optimized configuration
create_memory_config() {
    local mem_gb=$1
    local config_file=$2
    
    cat > "$config_file" << EOF
# Memory configuration for ${mem_gb}GB system
# Generated on $(date)

TOTAL_MEMORY_GB=$mem_gb
AVAILABLE_MEMORY_GB=$(get_available_memory)

# Tool selection
ALIGNER=$(select_aligner)
THREADS_ALIGNMENT=$(get_optimal_threads "alignment")
THREADS_COMPRESSION=$(get_optimal_threads "compression")
THREADS_QC=$(get_optimal_threads "fastqc")

# Memory limits (80% of available)
MAX_MEMORY_GB=$((mem_gb * 80 / 100))

# Tool-specific settings
if [[ $mem_gb -ge 8 ]]; then
    BOWTIE2_MODE="standard"
    CHUNK_SIZE=0  # No chunking needed
    PARALLEL_SAMPLES=1
elif [[ $mem_gb -ge 4 ]]; then
    BOWTIE2_MODE="lowmem"
    CHUNK_SIZE=10000000  # 10M reads per chunk
    PARALLEL_SAMPLES=1
else
    BOWTIE2_MODE="none"
    USE_MINIMAP2=true
    CHUNK_SIZE=5000000   # 5M reads per chunk
    PARALLEL_SAMPLES=1
fi

# Disk space requirements (GB)
TEMP_SPACE_REQUIRED=10
OUTPUT_SPACE_REQUIRED=5
EOF

    success "Memory configuration saved: $config_file"
}

# System resource summary
print_resource_summary() {
    local total_mem available_mem cpu_count disk_available
    total_mem=$(get_total_memory)
    available_mem=$(get_available_memory)
    cpu_count=$(get_cpu_count)
    disk_available=$(get_available_disk ".")
    
    echo ""
    echo "System Resource Summary"
    echo "======================="
    echo "Total Memory:     ${total_mem}GB"
    echo "Available Memory: ${available_mem}GB"
    echo "CPU Cores:        $cpu_count"
    echo "Available Disk:   ${disk_available}GB"
    echo ""
    echo "Recommended Settings"
    echo "===================="
    echo "Aligner:          $(select_aligner)"
    echo "Alignment Threads: $(get_optimal_threads "alignment")"
    echo "Compression Threads: $(get_optimal_threads "compression")"
    echo ""
}

# Check system requirements
check_system_requirements() {
    local min_memory=${1:-2}
    local min_disk=${2:-10}
    local requirements_met=true
    
    log "Checking system requirements..."
    
    # Check memory
    local available_mem
    available_mem=$(get_available_memory)
    if [[ $available_mem -lt $min_memory ]]; then
        error "Insufficient memory: ${available_mem}GB available, ${min_memory}GB required"
        requirements_met=false
    else
        success "Memory requirement met: ${available_mem}GB available"
    fi
    
    # Check disk space
    local available_disk
    available_disk=$(get_available_disk ".")
    if [[ $available_disk -lt $min_disk ]]; then
        error "Insufficient disk space: ${available_disk}GB available, ${min_disk}GB required"
        requirements_met=false
    else
        success "Disk space requirement met: ${available_disk}GB available"
    fi
    
    if [[ "$requirements_met" == "true" ]]; then
        success "All system requirements met!"
        return 0
    else
        error "System requirements not met!"
        return 1
    fi
}

# Get memory usage for a specific process
get_process_memory() {
    local pid=$1
    if [[ -f "/proc/$pid/status" ]]; then
        grep VmRSS "/proc/$pid/status" | awk '{print int($2/1024/1024)}'
    else
        echo 0
    fi
}

# Wait for memory to become available
wait_for_memory() {
    local required_gb=$1
    local timeout=${2:-300}  # 5 minutes default
    local elapsed=0
    
    log "Waiting for ${required_gb}GB memory to become available..."
    
    while [[ $elapsed -lt $timeout ]]; do
        if check_memory_requirement "$required_gb"; then
            success "Required memory is now available"
            return 0
        fi
        
        sleep 10
        elapsed=$((elapsed + 10))
        
        if [[ $((elapsed % 60)) -eq 0 ]]; then
            log "Still waiting... ${elapsed}s elapsed ($(get_available_memory)GB available)"
        fi
    done
    
    error "Timeout waiting for memory to become available"
    return 1
}

# Main function for CLI usage
main() {
    case "${1:-summary}" in
        "summary")
            print_resource_summary
            ;;
        "memory")
            echo "$(get_available_memory)"
            ;;
        "threads")
            echo "$(get_optimal_threads "${2:-default}")"
            ;;
        "aligner")
            echo "$(select_aligner)"
            ;;
        "check")
            check_system_requirements "${2:-2}" "${3:-10}"
            ;;
        "config")
            create_memory_config "$(get_total_memory)" "${2:-memory_config.conf}"
            ;;
        *)
            echo "Usage: $0 [summary|memory|threads|aligner|check|config]"
            echo ""
            echo "Commands:"
            echo "  summary     Show system resource summary"
            echo "  memory      Show available memory in GB"
            echo "  threads     Show optimal thread count"
            echo "  aligner     Show recommended aligner"
            echo "  check       Check system requirements"
            echo "  config      Create memory configuration file"
            ;;
    esac
}

# Run main function if script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi