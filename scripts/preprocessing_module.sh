#!/bin/bash

# RNA-seq Preprocessing Module
# Author: Joshua Slysz, PhD - Dalhousie University
# Handles adapter trimming, quality filtering, and rRNA removal for bulk RNA-seq

set -euo pipefail

# Source resource manager
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/resource_manager.sh"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Logging functions
log() { echo -e "${BLUE}[PREPROCESS]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }

# Default parameters
READ1=""
READ2=""
OUTPUT_DIR=""
SAMPLE_ID=""
THREADS=0
MEMORY_GB=0
TEMP_DIR=""
SINGLE_END=false
ADAPTER_TRIMMING=true
RRNA_REMOVAL=false
MIN_LENGTH=36
QUALITY_THRESHOLD=20
VERBOSE=false

# Tool selection
TRIMMER="auto"  # auto, trimmomatic, cutadapt, fastp
RRNA_TOOL="auto"  # auto, sortmerna, bowtie2

# Parse arguments
parse_preprocessing_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --read1)
                READ1="$2"
                shift 2
                ;;
            --read2)
                READ2="$2"
                shift 2
                ;;
            --output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --sample-id)
                SAMPLE_ID="$2"
                shift 2
                ;;
            --threads)
                THREADS="$2"
                shift 2
                ;;
            --memory)
                MEMORY_GB="$2"
                shift 2
                ;;
            --temp-dir)
                TEMP_DIR="$2"
                shift 2
                ;;
            --single-end)
                SINGLE_END=true
                shift
                ;;
            --no-adapter-trimming)
                ADAPTER_TRIMMING=false
                shift
                ;;
            --rrna-removal)
                RRNA_REMOVAL=true
                shift
                ;;
            --min-length)
                MIN_LENGTH="$2"
                shift 2
                ;;
            --quality-threshold)
                QUALITY_THRESHOLD="$2"
                shift 2
                ;;
            --trimmer)
                TRIMMER="$2"
                shift 2
                ;;
            --verbose)
                VERBOSE=true
                shift
                ;;
            *)
                error "Unknown preprocessing option: $1"
                exit 1
                ;;
        esac
    done
}

# Validate parameters
validate_preprocessing_params() {
    local errors=0
    
    if [[ -z "$READ1" || ! -f "$READ1" ]]; then
        error "Read 1 file not found: $READ1"
        errors=$((errors + 1))
    fi
    
    if [[ "$SINGLE_END" == "false" && (-z "$READ2" || ! -f "$READ2") ]]; then
        error "Read 2 file not found: $READ2"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        error "Output directory must be specified"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$SAMPLE_ID" ]]; then
        error "Sample ID must be specified"
        errors=$((errors + 1))
    fi
    
    if [[ $errors -gt 0 ]]; then
        exit 1
    fi
    
    # Create output and temp directories
    mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"
}

# Auto-detect best trimming tool
auto_detect_trimmer() {
    if [[ "$TRIMMER" != "auto" ]]; then
        return 0
    fi
    
    # Check available tools in order of preference
    if command -v fastp >/dev/null 2>&1; then
        TRIMMER="fastp"
        log "Selected trimmer: fastp (fast and comprehensive)"
    elif command -v trimmomatic >/dev/null 2>&1; then
        TRIMMER="trimmomatic"
        log "Selected trimmer: trimmomatic (robust and well-tested)"
    elif command -v cutadapt >/dev/null 2>&1; then
        TRIMMER="cutadapt"
        log "Selected trimmer: cutadapt (simple and effective)"
    else
        warning "No adapter trimming tool found"
        warning "Install one of: fastp, trimmomatic, cutadapt"
        ADAPTER_TRIMMING=false
        TRIMMER="none"
    fi
}

# Set preprocessing parameters
set_preprocessing_params() {
    if [[ $THREADS -eq 0 ]]; then
        THREADS=$(get_optimal_threads "preprocessing")
    fi
    
    if [[ $MEMORY_GB -eq 0 ]]; then
        MEMORY_GB=$(get_available_memory)
    fi
    
    if [[ -z "$TEMP_DIR" ]]; then
        TEMP_DIR="/tmp/preprocessing_$$"
        mkdir -p "$TEMP_DIR"
    fi
    
    auto_detect_trimmer
    
    log "Preprocessing parameters:"
    log "  Sample: $SAMPLE_ID"
    log "  Adapter trimming: $ADAPTER_TRIMMING ($TRIMMER)"
    log "  rRNA removal: $RRNA_REMOVAL"
    log "  Threads: $THREADS"
    log "  Quality threshold: $QUALITY_THRESHOLD"
    log "  Min length: $MIN_LENGTH"
}

# Run fastp trimming
run_fastp_trimming() {
    log "Running fastp adapter trimming and quality filtering"
    
    local output1="${OUTPUT_DIR}/${SAMPLE_ID}_trimmed_1.fq.gz"
    local output2="${OUTPUT_DIR}/${SAMPLE_ID}_trimmed_2.fq.gz"
    local report="${OUTPUT_DIR}/${SAMPLE_ID}_fastp.html"
    local json="${OUTPUT_DIR}/${SAMPLE_ID}_fastp.json"
    
    local fastp_opts="--thread $THREADS --qualified_quality_phred $QUALITY_THRESHOLD"
    fastp_opts="$fastp_opts --length_required $MIN_LENGTH"
    fastp_opts="$fastp_opts --detect_adapter_for_pe"
    fastp_opts="$fastp_opts --correction --cut_front --cut_tail"
    fastp_opts="$fastp_opts --html $report --json $json"
    
    if [[ "$SINGLE_END" == "true" ]]; then
        output1="${OUTPUT_DIR}/${SAMPLE_ID}_trimmed.fq.gz"
        fastp_opts="$fastp_opts --in1 $READ1 --out1 $output1"
    else
        fastp_opts="$fastp_opts --in1 $READ1 --in2 $READ2 --out1 $output1 --out2 $output2"
    fi
    
    log "Executing: fastp $fastp_opts"
    
    if eval "fastp $fastp_opts"; then
        success "fastp trimming completed"
        # Update read paths for downstream processing
        READ1="$output1"
        if [[ "$SINGLE_END" == "false" ]]; then
            READ2="$output2"
        fi
    else
        error "fastp trimming failed"
        exit 1
    fi
}

# Run trimmomatic trimming
run_trimmomatic_trimming() {
    log "Running trimmomatic adapter trimming and quality filtering"
    
    local output1="${OUTPUT_DIR}/${SAMPLE_ID}_trimmed_1.fq.gz"
    local output2="${OUTPUT_DIR}/${SAMPLE_ID}_trimmed_2.fq.gz"
    local unpaired1="${TEMP_DIR}/${SAMPLE_ID}_unpaired_1.fq.gz"
    local unpaired2="${TEMP_DIR}/${SAMPLE_ID}_unpaired_2.fq.gz"
    
    # Create adapter file if not exists
    local adapter_file="${TEMP_DIR}/adapters.fa"
    cat > "$adapter_file" << 'EOF'
>TruSeq3_IndexedAdapter
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>TruSeq3_UniversalAdapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
EOF
    
    local trimmomatic_opts="ILLUMINACLIP:${adapter_file}:2:30:10"
    trimmomatic_opts="$trimmomatic_opts LEADING:3 TRAILING:3"
    trimmomatic_opts="$trimmomatic_opts SLIDINGWINDOW:4:${QUALITY_THRESHOLD}"
    trimmomatic_opts="$trimmomatic_opts MINLEN:${MIN_LENGTH}"
    
    if [[ "$SINGLE_END" == "true" ]]; then
        output1="${OUTPUT_DIR}/${SAMPLE_ID}_trimmed.fq.gz"
        if trimmomatic SE -threads "$THREADS" "$READ1" "$output1" $trimmomatic_opts; then
            success "trimmomatic trimming completed"
            READ1="$output1"
        else
            error "trimmomatic trimming failed"
            exit 1
        fi
    else
        if trimmomatic PE -threads "$THREADS" \
           "$READ1" "$READ2" \
           "$output1" "$unpaired1" \
           "$output2" "$unpaired2" \
           $trimmomatic_opts; then
            success "trimmomatic trimming completed"
            READ1="$output1"
            READ2="$output2"
        else
            error "trimmomatic trimming failed"
            exit 1
        fi
    fi
}

# Run cutadapt trimming
run_cutadapt_trimming() {
    log "Running cutadapt adapter trimming"
    
    local output1="${OUTPUT_DIR}/${SAMPLE_ID}_trimmed_1.fq.gz"
    local output2="${OUTPUT_DIR}/${SAMPLE_ID}_trimmed_2.fq.gz"
    
    local cutadapt_opts="-j $THREADS -q $QUALITY_THRESHOLD -m $MIN_LENGTH"
    cutadapt_opts="$cutadapt_opts -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    cutadapt_opts="$cutadapt_opts -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
    
    if [[ "$SINGLE_END" == "true" ]]; then
        output1="${OUTPUT_DIR}/${SAMPLE_ID}_trimmed.fq.gz"
        if cutadapt $cutadapt_opts -o "$output1" "$READ1"; then
            success "cutadapt trimming completed"
            READ1="$output1"
        else
            error "cutadapt trimming failed"
            exit 1
        fi
    else
        if cutadapt $cutadapt_opts -o "$output1" -p "$output2" "$READ1" "$READ2"; then
            success "cutadapt trimming completed"
            READ1="$output1"
            READ2="$output2"
        else
            error "cutadapt trimming failed"
            exit 1
        fi
    fi
}

# Run adapter trimming based on selected tool
run_adapter_trimming() {
    if [[ "$ADAPTER_TRIMMING" == "false" ]]; then
        log "Skipping adapter trimming"
        return 0
    fi
    
    case "$TRIMMER" in
        "fastp")
            run_fastp_trimming
            ;;
        "trimmomatic")
            run_trimmomatic_trimming
            ;;
        "cutadapt")
            run_cutadapt_trimming
            ;;
        "none")
            log "No trimming tool available, skipping adapter trimming"
            ;;
        *)
            error "Unknown trimmer: $TRIMMER"
            exit 1
            ;;
    esac
}

# Run rRNA removal (placeholder - implement if needed)
run_rrna_removal() {
    if [[ "$RRNA_REMOVAL" == "false" ]]; then
        log "Skipping rRNA removal"
        return 0
    fi
    
    warning "rRNA removal not yet implemented"
    warning "Consider using SortMeRNA or SILVA databases"
    
    # TODO: Implement SortMeRNA or bowtie2-based rRNA removal
    # This would involve:
    # 1. Download rRNA databases (SILVA, GreenGenes, etc.)
    # 2. Index rRNA sequences
    # 3. Align reads to rRNA and remove matches
    # 4. Keep non-rRNA reads for downstream analysis
}

# Generate preprocessing statistics
generate_preprocessing_stats() {
    log "Generating preprocessing statistics"
    
    local stats_file="${OUTPUT_DIR}/${SAMPLE_ID}_preprocessing_stats.txt"
    
    cat > "$stats_file" << EOF
PREPROCESSING STATISTICS
=======================
Sample: $SAMPLE_ID
Date: $(date)

Input files:
$(if [[ "$SINGLE_END" == "true" ]]; then
    echo "- Read 1: $(basename "$READ1")"
else
    echo "- Read 1: $(basename "$READ1")"
    echo "- Read 2: $(basename "$READ2")"
fi)

Processing steps:
- Adapter trimming: $ADAPTER_TRIMMING ($TRIMMER)
- rRNA removal: $RRNA_REMOVAL
- Quality threshold: $QUALITY_THRESHOLD
- Minimum length: $MIN_LENGTH

Output files:
$(if [[ "$SINGLE_END" == "true" ]]; then
    echo "- Trimmed reads: ${SAMPLE_ID}_trimmed.fq.gz"
else
    echo "- Trimmed reads 1: ${SAMPLE_ID}_trimmed_1.fq.gz"
    echo "- Trimmed reads 2: ${SAMPLE_ID}_trimmed_2.fq.gz"
fi)

EOF
    
    # Add read counts if possible
    if command -v seqkit >/dev/null 2>&1; then
        echo "Read counts:" >> "$stats_file"
        if [[ "$SINGLE_END" == "true" ]]; then
            local count
            count=$(seqkit stats -T "$READ1" | tail -n 1 | cut -f4)
            echo "- Total reads: $count" >> "$stats_file"
        else
            local count1 count2
            count1=$(seqkit stats -T "$READ1" | tail -n 1 | cut -f4)
            count2=$(seqkit stats -T "$READ2" | tail -n 1 | cut -f4)
            echo "- Read 1: $count1" >> "$stats_file"
            echo "- Read 2: $count2" >> "$stats_file"
        fi
    fi
    
    success "Preprocessing statistics saved to: $stats_file"
}

# Main preprocessing function
run_preprocessing() {
    parse_preprocessing_args "$@"
    set_preprocessing_params
    validate_preprocessing_params
    
    log "Starting preprocessing for sample: $SAMPLE_ID"
    
    # Run processing steps
    run_adapter_trimming
    run_rrna_removal
    
    # Generate statistics
    generate_preprocessing_stats
    
    # Output final read paths for downstream processing
    echo "PROCESSED_READ1=$READ1"
    if [[ "$SINGLE_END" == "false" ]]; then
        echo "PROCESSED_READ2=$READ2"
    fi
    
    success "Preprocessing completed for sample: $SAMPLE_ID"
    return 0
}

# Usage function
preprocessing_usage() {
    cat << EOF
RNA-seq Preprocessing Module Usage:

$0 --read1 FILE --sample-id ID --output-dir DIR [options]

Required arguments:
  --read1 FILE          Input FASTQ file (read 1)
  --sample-id ID        Sample identifier
  --output-dir DIR      Output directory for processed files

Optional arguments:
  --read2 FILE          Input FASTQ file (read 2, for paired-end)
  --single-end          Single-end sequencing mode
  --threads N           Number of threads [default: auto-detect]
  --memory N            Available memory in GB [default: auto-detect]
  --temp-dir DIR        Temporary directory
  --no-adapter-trimming Skip adapter trimming
  --rrna-removal        Enable rRNA removal (experimental)
  --min-length N        Minimum read length after trimming [default: 36]
  --quality-threshold N Quality threshold for trimming [default: 20]
  --trimmer TOOL        Trimming tool (fastp|trimmomatic|cutadapt|auto) [default: auto]
  --verbose             Verbose output

Examples:
  # Basic preprocessing with auto-detected tools
  $0 --read1 sample_1.fq.gz --read2 sample_2.fq.gz \\
     --sample-id sample --output-dir preprocessed/

  # Single-end with specific trimmer
  $0 --read1 sample.fq.gz --single-end \\
     --sample-id sample --output-dir preprocessed/ --trimmer fastp

EOF
}

# Run preprocessing if script is called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    if [[ $# -eq 0 ]]; then
        preprocessing_usage
        exit 1
    fi
    
    run_preprocessing "$@"
fi