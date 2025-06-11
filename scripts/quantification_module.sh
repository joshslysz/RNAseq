#!/bin/bash

# Quantification Module
# Author: Joshua Slysz, PhD - Dalhousie University
# Memory-efficient gene expression quantification using featureCounts and Salmon
# Optimized for local systems with limited resources

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
log() { echo -e "${BLUE}[QUANT]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }

# Default parameters
QUANTIFIER=""
BAM_FILE=""
GENOME_GTF=""
OUTPUT_FILE=""
SAMPLE_ID=""
THREADS=0
MEMORY_GB=0
TEMP_DIR=""
SINGLE_END=false
VERBOSE=false

# featureCounts specific options
FEATURE_TYPE="exon"
ATTRIBUTE_TYPE="gene_id"
STRAND_SPECIFIC=0  # 0=unstranded, 1=stranded, 2=reverse-stranded

# Salmon specific options
SALMON_INDEX=""
READ1=""
READ2=""

# Statistics
STATS_FILE=""

# Parse quantification module arguments
parse_quantification_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --quantifier)
                QUANTIFIER="$2"
                shift 2
                ;;
            --bam-file)
                BAM_FILE="$2"
                shift 2
                ;;
            --genome-gtf)
                GENOME_GTF="$2"
                shift 2
                ;;
            --output)
                OUTPUT_FILE="$2"
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
            --feature-type)
                FEATURE_TYPE="$2"
                shift 2
                ;;
            --attribute-type)
                ATTRIBUTE_TYPE="$2"
                shift 2
                ;;
            --strand-specific)
                STRAND_SPECIFIC="$2"
                shift 2
                ;;
            --salmon-index)
                SALMON_INDEX="$2"
                shift 2
                ;;
            --read1)
                READ1="$2"
                shift 2
                ;;
            --read2)
                READ2="$2"
                shift 2
                ;;
            --verbose)
                VERBOSE=true
                shift
                ;;
            *)
                error "Unknown quantification option: $1"
                exit 1
                ;;
        esac
    done
}

# Validate quantification parameters
validate_quantification_params() {
    local errors=0
    
    if [[ -z "$QUANTIFIER" ]]; then
        error "Quantifier must be specified with --quantifier"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$SAMPLE_ID" ]]; then
        error "Sample ID must be specified with --sample-id"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$OUTPUT_FILE" ]]; then
        error "Output file must be specified with --output"
        errors=$((errors + 1))
    fi
    
    # Validate quantifier-specific parameters
    if [[ "$QUANTIFIER" == "featurecounts" ]]; then
        if [[ -z "$BAM_FILE" ]]; then
            error "BAM file must be specified with --bam-file for featureCounts"
            errors=$((errors + 1))
        elif [[ ! -f "$BAM_FILE" ]]; then
            error "BAM file not found: $BAM_FILE"
            errors=$((errors + 1))
        fi
        
        if [[ -z "$GENOME_GTF" ]]; then
            error "GTF file must be specified with --genome-gtf for featureCounts"
            errors=$((errors + 1))
        elif [[ ! -f "$GENOME_GTF" ]]; then
            error "GTF file not found: $GENOME_GTF"
            errors=$((errors + 1))
        fi
    elif [[ "$QUANTIFIER" == "salmon" ]]; then
        if [[ -z "$SALMON_INDEX" ]]; then
            error "Salmon index must be specified with --salmon-index"
            errors=$((errors + 1))
        elif [[ ! -d "$SALMON_INDEX" ]]; then
            error "Salmon index directory not found: $SALMON_INDEX"
            errors=$((errors + 1))
        fi
        
        if [[ -z "$READ1" ]]; then
            error "Read 1 file must be specified with --read1 for Salmon"
            errors=$((errors + 1))
        elif [[ ! -f "$READ1" ]]; then
            error "Read 1 file not found: $READ1"
            errors=$((errors + 1))
        fi
        
        if [[ "$SINGLE_END" == "false" ]]; then
            if [[ -z "$READ2" ]]; then
                error "Read 2 file must be specified with --read2 for paired-end Salmon"
                errors=$((errors + 1))
            elif [[ ! -f "$READ2" ]]; then
                error "Read 2 file not found: $READ2"
                errors=$((errors + 1))
            fi
        fi
    else
        error "Unknown quantifier: $QUANTIFIER"
        error "Supported quantifiers: featurecounts, salmon"
        errors=$((errors + 1))
    fi
    
    if [[ $errors -gt 0 ]]; then
        error "Quantification parameter validation failed"
        exit 1
    fi
}

# Set quantification parameters
set_quantification_params() {
    # Auto-detect system resources if not specified
    if [[ $THREADS -eq 0 ]]; then
        THREADS=$(get_optimal_threads "quantification")
    fi
    
    if [[ $MEMORY_GB -eq 0 ]]; then
        MEMORY_GB=$(get_available_memory)
    fi
    
    # Create temporary directory if not specified
    if [[ -z "$TEMP_DIR" ]]; then
        TEMP_DIR="/tmp/quantification_$$"
        mkdir -p "$TEMP_DIR"
    fi
    
    # Set statistics file
    STATS_FILE="${TEMP_DIR}/${SAMPLE_ID}_quantification_stats.txt"
    
    log "Quantification parameters:"
    log "  Quantifier: $QUANTIFIER"
    log "  Threads: $THREADS"
    log "  Memory: ${MEMORY_GB}GB"
    log "  Sample: $SAMPLE_ID"
    log "  Single-end: $SINGLE_END"
}

# Run featureCounts quantification
run_featurecounts_quantification() {
    log "Running featureCounts quantification for sample: $SAMPLE_ID"
    
    # Check memory requirements (featureCounts is memory efficient)
    if ! check_memory_requirement 1; then
        error "Insufficient memory for featureCounts"
        exit 1
    fi
    
    # Set featureCounts parameters
    local fc_opts="-T $THREADS"
    local log_file="${TEMP_DIR}/${SAMPLE_ID}_featurecounts.log"
    
    # Feature and attribute type
    fc_opts="$fc_opts -t $FEATURE_TYPE -g $ATTRIBUTE_TYPE"
    
    # Strand specificity
    fc_opts="$fc_opts -s $STRAND_SPECIFIC"
    
    # Paired-end options
    if [[ "$SINGLE_END" == "false" ]]; then
        fc_opts="$fc_opts -p"  # Paired-end
        fc_opts="$fc_opts -B"  # Count both reads in a pair
    fi
    
    # Additional options for better quantification
    fc_opts="$fc_opts -C"  # Do not count chimeric fragments
    
    # Memory optimization: process fragments in chunks if needed
    if [[ $MEMORY_GB -lt 4 ]]; then
        fc_opts="$fc_opts --ignoreDup"  # Ignore duplicate reads to save memory
    fi
    
    # Output file and annotation
    fc_opts="$fc_opts -a $GENOME_GTF -o $OUTPUT_FILE"
    
    log "featureCounts options: $fc_opts"
    
    # Progress tracking
    echo "Starting featureCounts quantification at $(date)" > "$log_file"
    
    # Run featureCounts
    local fc_cmd="featureCounts $fc_opts $BAM_FILE"
    
    log "Executing: $fc_cmd"
    
    if eval "$fc_cmd" > "$log_file" 2>&1; then
        success "featureCounts quantification completed successfully"
    else
        error "featureCounts quantification failed"
        if [[ -f "$log_file" ]]; then
            error "Check log file: $log_file"
            tail -10 "$log_file"
        fi
        exit 1
    fi
    
    # Extract quantification statistics
    extract_featurecounts_stats "$log_file"
    
    # Create simplified count file (gene_id and count only)
    create_simplified_count_file
}

# Run Salmon quantification
run_salmon_quantification() {
    log "Running Salmon quantification for sample: $SAMPLE_ID"
    
    # Check memory requirements (Salmon is very memory efficient)
    if ! check_memory_requirement 1; then
        error "Insufficient memory for Salmon"
        exit 1
    fi
    
    local salmon_dir="${TEMP_DIR}/${SAMPLE_ID}_salmon"
    local log_file="${TEMP_DIR}/${SAMPLE_ID}_salmon.log"
    
    # Set Salmon parameters
    local salmon_opts="--threads $THREADS --validateMappings"
    
    # Library type auto-detection
    salmon_opts="$salmon_opts --libType A"
    
    # Bias correction
    salmon_opts="$salmon_opts --seqBias --gcBias"
    
    # Memory optimization for low-memory systems
    if [[ $MEMORY_GB -lt 4 ]]; then
        salmon_opts="$salmon_opts --reduceGCMemory"
    fi
    
    log "Salmon options: $salmon_opts"
    
    # Progress tracking
    echo "Starting Salmon quantification at $(date)" > "$log_file"
    
    # Run Salmon
    local salmon_cmd
    if [[ "$SINGLE_END" == "true" ]]; then
        salmon_cmd="salmon quant -i $SALMON_INDEX -r $READ1 -o $salmon_dir $salmon_opts"
    else
        salmon_cmd="salmon quant -i $SALMON_INDEX -1 $READ1 -2 $READ2 -o $salmon_dir $salmon_opts"
    fi
    
    log "Executing: $salmon_cmd"
    
    if eval "$salmon_cmd" >> "$log_file" 2>&1; then
        success "Salmon quantification completed successfully"
    else
        error "Salmon quantification failed"
        if [[ -f "$log_file" ]]; then
            error "Check log file: $log_file"
            tail -10 "$log_file"
        fi
        exit 1
    fi
    
    # Convert Salmon output to standard format
    convert_salmon_output "$salmon_dir"
    
    # Extract quantification statistics
    extract_salmon_stats "$salmon_dir" "$log_file"
}

# Extract featureCounts statistics
extract_featurecounts_stats() {
    local log_file=$1
    
    if [[ ! -f "$log_file" ]]; then
        warning "featureCounts log file not found: $log_file"
        return
    fi
    
    # Parse featureCounts summary
    local total_reads assigned_reads assignment_rate
    total_reads=$(grep "Total reads" "$log_file" | awk '{print $NF}' || echo "Unknown")
    assigned_reads=$(grep "Successfully assigned reads" "$log_file" | awk '{print $NF}' || echo "Unknown")
    
    if [[ "$total_reads" != "Unknown" && "$assigned_reads" != "Unknown" ]]; then
        assignment_rate=$(echo "scale=2; $assigned_reads * 100 / $total_reads" | bc -l 2>/dev/null || echo "N/A")
    else
        assignment_rate="Unknown"
    fi
    
    # Count total genes quantified
    local total_genes
    if [[ -f "$OUTPUT_FILE" ]]; then
        total_genes=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
    else
        total_genes="Unknown"
    fi
    
    # Write statistics
    cat > "$STATS_FILE" << EOF
Sample: $SAMPLE_ID
Quantifier: featureCounts
Total_reads: $total_reads
Assigned_reads: $assigned_reads
Assignment_rate_percent: $assignment_rate
Total_genes: $total_genes
Feature_type: $FEATURE_TYPE
Attribute_type: $ATTRIBUTE_TYPE
Strand_specific: $STRAND_SPECIFIC
Output_file: $OUTPUT_FILE
Log_file: $log_file
EOF
    
    log "Assignment rate: ${assignment_rate}%"
    log "Total genes quantified: $total_genes"
}

# Extract Salmon statistics
extract_salmon_stats() {
    local salmon_dir=$1
    local log_file=$2
    
    local meta_file="${salmon_dir}/aux_info/meta_info.json"
    local mapping_rate="Unknown"
    local total_reads="Unknown"
    
    # Extract mapping rate from Salmon output
    if [[ -f "$meta_file" ]]; then
        # Parse JSON for mapping rate (this is a simple extraction)
        mapping_rate=$(grep "percent_mapped" "$meta_file" | sed 's/.*: *\([0-9.]*\).*/\1/' || echo "Unknown")
    fi
    
    # Count total transcripts quantified
    local total_transcripts
    if [[ -f "${salmon_dir}/quant.sf" ]]; then
        total_transcripts=$(tail -n +2 "${salmon_dir}/quant.sf" | wc -l)
    else
        total_transcripts="Unknown"
    fi
    
    # Write statistics
    cat > "$STATS_FILE" << EOF
Sample: $SAMPLE_ID
Quantifier: Salmon
Mapping_rate_percent: $mapping_rate
Total_transcripts: $total_transcripts
Salmon_dir: $salmon_dir
Output_file: $OUTPUT_FILE
Log_file: $log_file
EOF
    
    log "Mapping rate: ${mapping_rate}%"
    log "Total transcripts quantified: $total_transcripts"
}

# Create simplified count file from featureCounts output
create_simplified_count_file() {
    if [[ ! -f "$OUTPUT_FILE" ]]; then
        warning "featureCounts output file not found: $OUTPUT_FILE"
        return
    fi
    
    # Create simplified file with just gene_id and count
    local simple_file="${OUTPUT_FILE%.txt}_simple.txt"
    
    # Extract gene_id (column 1) and count (last column)
    awk 'NR>1 {print $1 "\t" $NF}' "$OUTPUT_FILE" > "$simple_file"
    
    # Add header
    sed -i '1i gene_id\tcount' "$simple_file"
    
    log "Simplified count file created: $simple_file"
}

# Convert Salmon output to standard format
convert_salmon_output() {
    local salmon_dir=$1
    local quant_file="${salmon_dir}/quant.sf"
    
    if [[ ! -f "$quant_file" ]]; then
        error "Salmon quantification file not found: $quant_file"
        exit 1
    fi
    
    # Convert Salmon TPM to count-like format
    # Extract transcript_id and NumReads (estimated counts)
    awk 'NR>1 {print $1 "\t" $5}' "$quant_file" > "$OUTPUT_FILE"
    
    # Add header
    sed -i '1i transcript_id\tcount' "$OUTPUT_FILE"
    
    log "Salmon output converted to standard format: $OUTPUT_FILE"
}

# Check quantification quality
check_quantification_quality() {
    if [[ ! -f "$STATS_FILE" ]]; then
        warning "Cannot check quantification quality - stats file not found"
        return
    fi
    
    # Check assignment/mapping rate
    local rate_line rate_value
    
    if [[ "$QUANTIFIER" == "featurecounts" ]]; then
        rate_line=$(grep "Assignment_rate_percent:" "$STATS_FILE" || echo "")
        rate_value=$(echo "$rate_line" | cut -d: -f2 | tr -d ' ')
    else
        rate_line=$(grep "Mapping_rate_percent:" "$STATS_FILE" || echo "")
        rate_value=$(echo "$rate_line" | cut -d: -f2 | tr -d ' ')
    fi
    
    if [[ -n "$rate_value" && "$rate_value" != "Unknown" ]]; then
        local rate_int
        rate_int=$(echo "$rate_value" | cut -d. -f1)
        
        if [[ $rate_int -lt 30 ]]; then
            error "Low ${QUANTIFIER} rate detected: ${rate_value}%"
            error "Consider checking:"
            error "  - GTF file compatibility with genome"
            error "  - Strand specificity settings"
            error "  - Read quality and alignment"
        elif [[ $rate_int -lt 50 ]]; then
            warning "Moderate ${QUANTIFIER} rate: ${rate_value}%"
            warning "This may be acceptable depending on your experimental conditions"
        else
            success "Good ${QUANTIFIER} rate: ${rate_value}%"
        fi
    fi
    
    # Check output file
    if [[ -f "$OUTPUT_FILE" ]]; then
        local line_count
        line_count=$(wc -l < "$OUTPUT_FILE")
        if [[ $line_count -gt 1 ]]; then
            success "Output file contains $((line_count - 1)) quantified features"
        else
            error "Output file appears to be empty or malformed"
        fi
    fi
}

# Cleanup quantification temporary files
cleanup_quantification() {
    if [[ -n "$TEMP_DIR" && -d "$TEMP_DIR" ]]; then
        # Keep statistics file
        if [[ -f "$STATS_FILE" ]]; then
            local basename_file
            basename_file=$(basename "$STATS_FILE")
            cp "$STATS_FILE" "$(dirname "$OUTPUT_FILE")/${basename_file}"
        fi
        
        # Remove temporary directory
        rm -rf "$TEMP_DIR"
        log "Temporary quantification files cleaned"
    fi
}

# Main quantification function
run_quantification() {
    # Parse arguments
    parse_quantification_args "$@"
    
    # Set parameters
    set_quantification_params
    
    # Validate parameters
    validate_quantification_params
    
    log "Starting quantification process"
    log "Quantifier: $QUANTIFIER"
    log "Output: $OUTPUT_FILE"
    
    # Run quantification based on selected method
    case "$QUANTIFIER" in
        "featurecounts")
            run_featurecounts_quantification
            ;;
        "salmon")
            run_salmon_quantification
            ;;
        *)
            error "Unknown quantifier: $QUANTIFIER"
            exit 1
            ;;
    esac
    
    # Quality check
    check_quantification_quality
    
    # Cleanup
    cleanup_quantification
    
    success "Quantification module completed successfully"
    return 0
}

# Usage function for quantification module
quantification_usage() {
    cat << EOF
Quantification Module Usage:

$0 --quantifier TOOL --output FILE --sample-id ID [options]

Required arguments:
  --quantifier TOOL     Quantification method (featurecounts|salmon)
  --output FILE         Output count file path
  --sample-id ID        Sample identifier

featureCounts specific arguments:
  --bam-file FILE       Input BAM file
  --genome-gtf FILE     Gene annotation GTF file
  --feature-type TYPE   Feature type to count [default: exon]
  --attribute-type TYPE Attribute type for grouping [default: gene_id]
  --strand-specific N   Strand specificity (0=unstranded, 1=stranded, 2=reverse) [default: 0]

Salmon specific arguments:
  --salmon-index DIR    Salmon index directory
  --read1 FILE          Input FASTQ file (read 1)
  --read2 FILE          Input FASTQ file (read 2, for paired-end)

Optional arguments:
  --single-end          Single-end sequencing mode
  --threads N           Number of threads [default: auto-detect]
  --memory N            Available memory in GB [default: auto-detect]
  --temp-dir DIR        Temporary directory [default: auto-create]
  --verbose             Verbose output

Examples:
  # featureCounts quantification
  $0 --quantifier featurecounts --bam-file sample.bam --genome-gtf genes.gtf \\
     --output sample.counts.txt --sample-id sample

  # Salmon quantification
  $0 --quantifier salmon --salmon-index salmon_index/ --read1 sample_1.fq.gz \\
     --read2 sample_2.fq.gz --output sample.counts.txt --sample-id sample

  # Single-end with custom parameters
  $0 --quantifier featurecounts --bam-file sample.bam --genome-gtf genes.gtf \\
     --output sample.counts.txt --sample-id sample --single-end --strand-specific 1
EOF
}

# Run quantification if script is called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    if [[ $# -eq 0 ]]; then
        quantification_usage
        exit 1
    fi
    
    run_quantification "$@"
fi