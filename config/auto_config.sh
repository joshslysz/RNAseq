#!/bin/bash

# Automatic Configuration Selection
# Detects system resources and selects optimal configuration
# Part of memory-optimized RNA-seq pipeline

set -euo pipefail

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source resource manager
source "$SCRIPT_DIR/../scripts/resource_manager.sh"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Logging functions
log() { echo -e "${BLUE}[CONFIG]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }

# Default parameters
OUTPUT_FILE=""
FORCE_CONFIG=""
VERBOSE=false

# Usage function
usage() {
    cat << EOF
Automatic Configuration Selection

Usage: $0 [options]

Options:
  -o, --output FILE     Output configuration file [default: auto_generated.conf]
  -f, --force CONFIG    Force specific configuration (low|medium|high)
  -v, --verbose         Verbose output
  -h, --help            Show this help message

Examples:
  # Auto-detect and generate configuration
  $0 -o pipeline.conf

  # Force high memory configuration
  $0 -o pipeline.conf -f high

  # Verbose auto-detection
  $0 -o pipeline.conf -v
EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -f|--force)
            FORCE_CONFIG="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Set default output file if not specified
if [[ -z "$OUTPUT_FILE" ]]; then
    OUTPUT_FILE="auto_generated.conf"
fi

# System resource detection
detect_system_resources() {
    local total_mem available_mem cpu_count disk_space
    
    total_mem=$(get_total_memory)
    available_mem=$(get_available_memory)
    cpu_count=$(get_cpu_count)
    disk_space=$(get_available_disk ".")
    
    if [[ "$VERBOSE" == "true" ]]; then
        log "System resource detection:"
        log "  Total memory: ${total_mem}GB"
        log "  Available memory: ${available_mem}GB"
        log "  CPU cores: $cpu_count"
        log "  Available disk space: ${disk_space}GB"
    fi
    
    # Return values for configuration selection
    echo "$total_mem,$available_mem,$cpu_count,$disk_space"
}

# Select optimal configuration
select_configuration() {
    local resources=$1
    IFS=',' read -r total_mem available_mem cpu_count disk_space <<< "$resources"
    
    local config_type=""
    
    # Force configuration if specified
    if [[ -n "$FORCE_CONFIG" ]]; then
        case "$FORCE_CONFIG" in
            "low"|"medium"|"high")
                config_type="$FORCE_CONFIG"
                log "Forced configuration: $config_type"
                ;;
            *)
                error "Invalid forced configuration: $FORCE_CONFIG"
                error "Valid options: low, medium, high"
                exit 1
                ;;
        esac
    else
        # Auto-detect based on available memory
        if [[ $available_mem -ge 8 ]]; then
            config_type="high"
            log "Auto-selected: high memory configuration (${available_mem}GB available)"
        elif [[ $available_mem -ge 4 ]]; then
            config_type="medium"
            log "Auto-selected: medium memory configuration (${available_mem}GB available)"
        else
            config_type="low"
            log "Auto-selected: low memory configuration (${available_mem}GB available)"
        fi
    fi
    
    # Validate minimum requirements
    case "$config_type" in
        "low")
            if [[ $available_mem -lt 2 ]]; then
                error "Insufficient memory for low configuration: ${available_mem}GB (minimum 2GB)"
                exit 1
            fi
            ;;
        "medium")
            if [[ $available_mem -lt 4 ]]; then
                warning "Borderline memory for medium configuration: ${available_mem}GB"
                warning "Consider using low memory configuration"
            fi
            ;;
        "high")
            if [[ $available_mem -lt 8 ]]; then
                warning "Borderline memory for high configuration: ${available_mem}GB"
                warning "Consider using medium memory configuration"
            fi
            ;;
    esac
    
    echo "$config_type"
}

# Generate custom configuration
generate_configuration() {
    local config_type=$1
    local resources=$2
    IFS=',' read -r total_mem available_mem cpu_count disk_space <<< "$resources"
    
    # Load base configuration
    local base_config="$SCRIPT_DIR/${config_type}_memory.conf"
    
    if [[ ! -f "$base_config" ]]; then
        error "Base configuration file not found: $base_config"
        exit 1
    fi
    
    log "Generating configuration based on: $base_config"
    
    # Create header
    cat > "$OUTPUT_FILE" << EOF
# Auto-generated RNA-seq Pipeline Configuration
# Generated on: $(date)
# System: ${total_mem}GB total memory, ${available_mem}GB available, ${cpu_count} CPU cores
# Configuration type: $config_type

EOF
    
    # Copy base configuration
    cat "$base_config" >> "$OUTPUT_FILE"
    
    # Add system-specific overrides
    cat >> "$OUTPUT_FILE" << EOF

# System-specific overrides (auto-generated)
DETECTED_TOTAL_MEMORY_GB=$total_mem
DETECTED_AVAILABLE_MEMORY_GB=$available_mem
DETECTED_CPU_COUNT=$cpu_count
DETECTED_DISK_SPACE_GB=$disk_space
GENERATION_DATE="$(date)"
CONFIGURATION_TYPE="$config_type"
AUTO_GENERATED=true
EOF
    
    # Optimize thread count based on CPU cores and memory
    local optimal_threads
    case "$config_type" in
        "low")
            optimal_threads=1
            ;;
        "medium")
            optimal_threads=$((cpu_count > 2 ? 2 : cpu_count))
            ;;
        "high")
            optimal_threads=$((cpu_count > 4 ? 4 : cpu_count))
            ;;
    esac
    
    # Add optimized settings
    cat >> "$OUTPUT_FILE" << EOF

# Optimized settings based on detected hardware
OPTIMAL_THREADS=$optimal_threads
OPTIMAL_MEMORY_LIMIT=$((available_mem * 80 / 100))G
OPTIMAL_SAMTOOLS_THREADS=$((optimal_threads > 1 ? optimal_threads / 2 : 1))
EOF
    
    success "Configuration generated: $OUTPUT_FILE"
}

# Validate configuration
validate_configuration() {
    local config_file=$1
    
    if [[ ! -f "$config_file" ]]; then
        error "Configuration file not found: $config_file"
        exit 1
    fi
    
    # Check for required variables
    local required_vars=("ALIGNER" "QUANTIFIER" "MIN_MEMORY_GB" "RECOMMENDED_THREADS")
    local missing_vars=()
    
    for var in "${required_vars[@]}"; do
        if ! grep -q "^$var=" "$config_file"; then
            missing_vars+=("$var")
        fi
    done
    
    if [[ ${#missing_vars[@]} -gt 0 ]]; then
        error "Missing required configuration variables: ${missing_vars[*]}"
        exit 1
    fi
    
    success "Configuration validation passed"
}

# Display configuration summary
display_summary() {
    local config_file=$1
    
    if [[ "$VERBOSE" == "true" ]]; then
        log "Configuration summary:"
        echo ""
        
        # Extract key settings
        local aligner quantifier memory_gb threads
        aligner=$(grep "^ALIGNER=" "$config_file" | cut -d= -f2 || echo "Unknown")
        quantifier=$(grep "^QUANTIFIER=" "$config_file" | cut -d= -f2 || echo "Unknown")
        memory_gb=$(grep "^MIN_MEMORY_GB=" "$config_file" | cut -d= -f2 || echo "Unknown")
        threads=$(grep "^RECOMMENDED_THREADS=" "$config_file" | cut -d= -f2 || echo "Unknown")
        
        echo "  Aligner: $aligner"
        echo "  Quantifier: $quantifier"
        echo "  Memory requirement: ${memory_gb}GB"
        echo "  Recommended threads: $threads"
        echo ""
        echo "Configuration file: $config_file"
        echo ""
    fi
}

# Main function
main() {
    log "Starting automatic configuration selection"
    
    # Detect system resources
    local resources
    resources=$(detect_system_resources)
    
    # Select configuration
    local config_type
    config_type=$(select_configuration "$resources")
    
    # Generate configuration
    generate_configuration "$config_type" "$resources"
    
    # Validate configuration
    validate_configuration "$OUTPUT_FILE"
    
    # Display summary
    display_summary "$OUTPUT_FILE"
    
    success "Configuration selection completed successfully"
    
    # Usage instructions
    echo ""
    echo "To use this configuration:"
    echo "  ./rna_preprocessing_local.sh --genome-config your_genome_config.txt --data-dir your_data/"
    echo ""
    echo "Or source the configuration in your script:"
    echo "  source $OUTPUT_FILE"
}

# Run main function
main