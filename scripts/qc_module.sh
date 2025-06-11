#!/bin/bash

# Quality Control Module
# Author: Joshua Slysz, PhD - Dalhousie University
# Comprehensive QC analysis with FastQC, MultiQC, and RSeQC
# Memory-optimized for local systems

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
log() { echo -e "${BLUE}[QC]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }

# Default parameters
QC_TYPE=""
FASTQ_FILES=()
BAM_FILE=""
OUTPUT_DIR=""
SAMPLE_ID=""
THREADS=0
MEMORY_GB=0
TEMP_DIR=""
SINGLE_END=false
VERBOSE=false

# Reference files for RSeQC
GENOME_GTF=""
BED_FILE=""

# QC specific options
FASTQC_KMERS=7
RUN_RSEQC=true
RUN_FASTQC=true
RUN_MULTIQC=true

# Statistics
STATS_FILE=""

# Parse QC module arguments
parse_qc_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --qc-type)
                QC_TYPE="$2"
                shift 2
                ;;
            --fastq-files)
                IFS=',' read -ra FASTQ_FILES <<< "$2"
                shift 2
                ;;
            --bam-file)
                BAM_FILE="$2"
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
            --genome-gtf)
                GENOME_GTF="$2"
                shift 2
                ;;
            --bed-file)
                BED_FILE="$2"
                shift 2
                ;;
            --fastqc-kmers)
                FASTQC_KMERS="$2"
                shift 2
                ;;
            --skip-rseqc)
                RUN_RSEQC=false
                shift
                ;;
            --skip-fastqc)
                RUN_FASTQC=false
                shift
                ;;
            --skip-multiqc)
                RUN_MULTIQC=false
                shift
                ;;
            --verbose)
                VERBOSE=true
                shift
                ;;
            *)
                error "Unknown QC option: $1"
                exit 1
                ;;
        esac
    done
}

# Validate QC parameters
validate_qc_params() {
    local errors=0
    
    if [[ -z "$QC_TYPE" ]]; then
        error "QC type must be specified with --qc-type (fastq|bam|all)"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        error "Output directory must be specified with --output-dir"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$SAMPLE_ID" ]]; then
        error "Sample ID must be specified with --sample-id"
        errors=$((errors + 1))
    fi
    
    # Validate QC type specific parameters
    if [[ "$QC_TYPE" == "fastq" || "$QC_TYPE" == "all" ]]; then
        if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
            error "FASTQ files must be specified with --fastq-files for FASTQ QC"
            errors=$((errors + 1))
        else
            for fastq in "${FASTQ_FILES[@]}"; do
                if [[ ! -f "$fastq" ]]; then
                    error "FASTQ file not found: $fastq"
                    errors=$((errors + 1))
                fi
            done
        fi
    fi
    
    if [[ "$QC_TYPE" == "bam" || "$QC_TYPE" == "all" ]]; then
        if [[ -z "$BAM_FILE" ]]; then
            error "BAM file must be specified with --bam-file for BAM QC"
            errors=$((errors + 1))
        elif [[ ! -f "$BAM_FILE" ]]; then
            error "BAM file not found: $BAM_FILE"
            errors=$((errors + 1))
        fi
        
        if [[ "$RUN_RSEQC" == "true" ]]; then
            if [[ -z "$GENOME_GTF" ]]; then
                error "GTF file must be specified with --genome-gtf for RSeQC"
                errors=$((errors + 1))
            elif [[ ! -f "$GENOME_GTF" ]]; then
                error "GTF file not found: $GENOME_GTF"
                errors=$((errors + 1))
            fi
        fi
    fi
    
    if [[ $errors -gt 0 ]]; then
        error "QC parameter validation failed"
        exit 1
    fi
}

# Set QC parameters
set_qc_params() {
    # Auto-detect system resources if not specified
    if [[ $THREADS -eq 0 ]]; then
        THREADS=$(get_optimal_threads "fastqc")
    fi
    
    if [[ $MEMORY_GB -eq 0 ]]; then
        MEMORY_GB=$(get_available_memory)
    fi
    
    # Create output directories
    mkdir -p "$OUTPUT_DIR"/{fastqc,rseqc,multiqc}
    OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
    
    # Create temporary directory if not specified
    if [[ -z "$TEMP_DIR" ]]; then
        TEMP_DIR="/tmp/qc_$$"
        mkdir -p "$TEMP_DIR"
    fi
    
    # Set statistics file
    STATS_FILE="${TEMP_DIR}/${SAMPLE_ID}_qc_stats.txt"
    
    # Disable resource-intensive QC for very low memory systems
    if [[ $MEMORY_GB -lt 2 ]]; then
        warning "Low memory detected (${MEMORY_GB}GB), disabling RSeQC"
        RUN_RSEQC=false
    fi
    
    log "QC parameters:"
    log "  QC type: $QC_TYPE"
    log "  Threads: $THREADS"
    log "  Memory: ${MEMORY_GB}GB"
    log "  Sample: $SAMPLE_ID"
    log "  Output: $OUTPUT_DIR"
    log "  Run FastQC: $RUN_FASTQC"
    log "  Run RSeQC: $RUN_RSEQC"
    log "  Run MultiQC: $RUN_MULTIQC"
}

# Run FastQC analysis
run_fastqc_analysis() {
    if [[ "$RUN_FASTQC" == "false" ]]; then
        log "Skipping FastQC analysis"
        return 0
    fi
    
    log "Running FastQC analysis for sample: $SAMPLE_ID"
    
    # Check memory requirements
    if ! check_memory_requirement 1; then
        warning "Insufficient memory for FastQC, skipping"
        return 0
    fi
    
    local fastqc_dir="$OUTPUT_DIR/fastqc"
    local log_file="${TEMP_DIR}/${SAMPLE_ID}_fastqc.log"
    
    # FastQC options
    local fastqc_opts="-o $fastqc_dir -t $THREADS"
    fastqc_opts="$fastqc_opts --kmers $FASTQC_KMERS"
    
    # Memory optimization for low-memory systems
    if [[ $MEMORY_GB -lt 4 ]]; then
        fastqc_opts="$fastqc_opts --nogroup"  # Disable position grouping to save memory
    fi
    
    log "FastQC options: $fastqc_opts"
    
    # Progress tracking
    echo "Starting FastQC analysis at $(date)" > "$log_file"
    
    # Run FastQC on all FASTQ files
    local failed_files=()
    for fastq in "${FASTQ_FILES[@]}"; do
        log "Processing FASTQ file: $(basename "$fastq")"
        
        if fastqc $fastqc_opts "$fastq" >> "$log_file" 2>&1; then
            success "FastQC completed for: $(basename "$fastq")"
        else
            error "FastQC failed for: $(basename "$fastq")"
            failed_files+=("$fastq")
        fi
    done
    
    # Report results
    if [[ ${#failed_files[@]} -eq 0 ]]; then
        success "FastQC analysis completed successfully for all files"
    else
        warning "FastQC failed for ${#failed_files[@]} files: ${failed_files[*]}"
    fi
    
    # Extract FastQC statistics
    extract_fastqc_stats "$fastqc_dir"
}

# Run RSeQC analysis
run_rseqc_analysis() {
    if [[ "$RUN_RSEQC" == "false" ]]; then
        log "Skipping RSeQC analysis"
        return 0
    fi
    
    log "Running RSeQC analysis for sample: $SAMPLE_ID"
    
    # Check memory requirements
    if ! check_memory_requirement 2; then
        warning "Insufficient memory for RSeQC, skipping"
        return 0
    fi
    
    local rseqc_dir="$OUTPUT_DIR/rseqc"
    local log_file="${TEMP_DIR}/${SAMPLE_ID}_rseqc.log"
    
    # Create BED file from GTF if not provided
    if [[ -z "$BED_FILE" ]]; then
        log "Converting GTF to BED format for RSeQC..."
        BED_FILE="${TEMP_DIR}/genes.bed"
        gtf_to_bed "$GENOME_GTF" "$BED_FILE"
    fi
    
    # Progress tracking
    echo "Starting RSeQC analysis at $(date)" > "$log_file"
    
    # RSeQC modules to run (memory-optimized selection)
    local rseqc_modules=()
    
    if [[ $MEMORY_GB -ge 4 ]]; then
        rseqc_modules+=("bam_stat" "read_distribution" "infer_experiment")
    fi
    
    if [[ $MEMORY_GB -ge 6 ]]; then
        rseqc_modules+=("junction_annotation" "read_duplication")
    fi
    
    if [[ $MEMORY_GB -ge 8 ]]; then
        rseqc_modules+=("inner_distance" "read_GC")
    fi
    
    log "Running RSeQC modules: ${rseqc_modules[*]}"
    
    # Run selected RSeQC modules
    local failed_modules=()
    for module in "${rseqc_modules[@]}"; do
        log "Running RSeQC module: $module"
        
        case "$module" in
            "bam_stat")
                if bam_stat.py -i "$BAM_FILE" > "$rseqc_dir/${SAMPLE_ID}.bam_stat.txt" 2>>"$log_file"; then
                    success "RSeQC bam_stat completed"
                else
                    failed_modules+=("bam_stat")
                fi
                ;;
            "read_distribution")
                if read_distribution.py -i "$BAM_FILE" -r "$BED_FILE" > "$rseqc_dir/${SAMPLE_ID}.read_distribution.txt" 2>>"$log_file"; then
                    success "RSeQC read_distribution completed"
                else
                    failed_modules+=("read_distribution")
                fi
                ;;
            "infer_experiment")
                if infer_experiment.py -i "$BAM_FILE" -r "$BED_FILE" > "$rseqc_dir/${SAMPLE_ID}.infer_experiment.txt" 2>>"$log_file"; then
                    success "RSeQC infer_experiment completed"
                else
                    failed_modules+=("infer_experiment")
                fi
                ;;
            "junction_annotation")
                if junction_annotation.py -i "$BAM_FILE" -r "$BED_FILE" -o "$rseqc_dir/${SAMPLE_ID}.junction" >> "$log_file" 2>&1; then
                    success "RSeQC junction_annotation completed"
                else
                    failed_modules+=("junction_annotation")
                fi
                ;;
            "read_duplication")
                if read_duplication.py -i "$BAM_FILE" -o "$rseqc_dir/${SAMPLE_ID}.read_dup" >> "$log_file" 2>&1; then
                    success "RSeQC read_duplication completed"
                else
                    failed_modules+=("read_duplication")
                fi
                ;;
            "inner_distance")
                if [[ "$SINGLE_END" == "false" ]]; then
                    if inner_distance.py -i "$BAM_FILE" -r "$BED_FILE" -o "$rseqc_dir/${SAMPLE_ID}.inner_distance" >> "$log_file" 2>&1; then
                        success "RSeQC inner_distance completed"
                    else
                        failed_modules+=("inner_distance")
                    fi
                else
                    log "Skipping inner_distance (single-end data)"
                fi
                ;;
            "read_GC")
                if read_GC.py -i "$BAM_FILE" -o "$rseqc_dir/${SAMPLE_ID}.read_GC" >> "$log_file" 2>&1; then
                    success "RSeQC read_GC completed"
                else
                    failed_modules+=("read_GC")
                fi
                ;;
        esac
    done
    
    # Report results
    if [[ ${#failed_modules[@]} -eq 0 ]]; then
        success "RSeQC analysis completed successfully for all modules"
    else
        warning "RSeQC failed for modules: ${failed_modules[*]}"
    fi
    
    # Extract RSeQC statistics
    extract_rseqc_stats "$rseqc_dir"
}

# Run MultiQC analysis
run_multiqc_analysis() {
    if [[ "$RUN_MULTIQC" == "false" ]]; then
        log "Skipping MultiQC analysis"
        return 0
    fi
    
    log "Running MultiQC analysis"
    
    # Check memory requirements
    if ! check_memory_requirement 1; then
        warning "Insufficient memory for MultiQC, skipping"
        return 0
    fi
    
    local multiqc_dir="$OUTPUT_DIR/multiqc"
    local log_file="${TEMP_DIR}/multiqc.log"
    
    # MultiQC options
    local multiqc_opts="-o $multiqc_dir -f"  # Force overwrite
    multiqc_opts="$multiqc_opts --filename ${SAMPLE_ID}_multiqc_report"
    
    # Include all QC output directories
    local input_dirs=("$OUTPUT_DIR")
    
    log "MultiQC options: $multiqc_opts"
    
    # Progress tracking
    echo "Starting MultiQC analysis at $(date)" > "$log_file"
    
    # Run MultiQC
    if multiqc $multiqc_opts "${input_dirs[@]}" >> "$log_file" 2>&1; then
        success "MultiQC analysis completed successfully"
    else
        error "MultiQC analysis failed"
        if [[ -f "$log_file" ]]; then
            error "Check log file: $log_file"
            tail -10 "$log_file"
        fi
        return 1
    fi
    
    # Check if report was generated
    local report_file="$multiqc_dir/${SAMPLE_ID}_multiqc_report.html"
    if [[ -f "$report_file" ]]; then
        success "MultiQC report generated: $report_file"
    else
        warning "MultiQC report not found at expected location"
    fi
}

# Convert GTF to BED format
gtf_to_bed() {
    local gtf_file=$1
    local bed_file=$2
    
    log "Converting GTF to BED format..."
    
    # Simple GTF to BED conversion (extract exons)
    awk '$3=="exon" {print $1 "\t" ($4-1) "\t" $5 "\t" $10 "\t" "0" "\t" $7}' "$gtf_file" | \
    sed 's/[";]//g' > "$bed_file"
    
    if [[ -f "$bed_file" && -s "$bed_file" ]]; then
        success "GTF converted to BED format: $bed_file"
    else
        error "Failed to convert GTF to BED format"
        exit 1
    fi
}

# Extract FastQC statistics
extract_fastqc_stats() {
    local fastqc_dir=$1
    
    # Count FastQC reports
    local fastqc_reports
    fastqc_reports=$(find "$fastqc_dir" -name "*_fastqc.html" | wc -l)
    
    # Initialize statistics
    cat > "$STATS_FILE" << EOF
Sample: $SAMPLE_ID
QC_type: FastQC
FastQC_reports: $fastqc_reports
FastQC_dir: $fastqc_dir
EOF
    
    # Extract basic statistics from FastQC data files
    local total_sequences=0
    local poor_quality=0
    
    for data_file in "$fastqc_dir"/*_fastqc/fastqc_data.txt; do
        if [[ -f "$data_file" ]]; then
            # Extract total sequences
            local seq_count
            seq_count=$(grep "Total Sequences" "$data_file" | awk '{print $3}' || echo "0")
            total_sequences=$((total_sequences + seq_count))
            
            # Check for poor quality warning
            if grep -q "FAIL" "$data_file"; then
                poor_quality=$((poor_quality + 1))
            fi
        fi
    done
    
    # Update statistics
    cat >> "$STATS_FILE" << EOF
Total_sequences: $total_sequences
Poor_quality_files: $poor_quality
EOF
    
    log "FastQC statistics: $fastqc_reports reports, $total_sequences total sequences"
}

# Extract RSeQC statistics
extract_rseqc_stats() {
    local rseqc_dir=$1
    
    # Count RSeQC output files
    local rseqc_files
    rseqc_files=$(find "$rseqc_dir" -name "${SAMPLE_ID}*" | wc -l)
    
    # Update statistics
    cat >> "$STATS_FILE" << EOF
RSeQC_files: $rseqc_files
RSeQC_dir: $rseqc_dir
EOF
    
    # Extract specific statistics if available
    local bam_stat_file="$rseqc_dir/${SAMPLE_ID}.bam_stat.txt"
    if [[ -f "$bam_stat_file" ]]; then
        local total_reads
        total_reads=$(grep "Total records:" "$bam_stat_file" | awk '{print $3}' || echo "Unknown")
        
        cat >> "$STATS_FILE" << EOF
RSeQC_total_reads: $total_reads
EOF
    fi
    
    log "RSeQC statistics: $rseqc_files output files generated"
}

# Check QC quality
check_qc_quality() {
    if [[ ! -f "$STATS_FILE" ]]; then
        warning "Cannot check QC quality - stats file not found"
        return
    fi
    
    # Check FastQC results
    if grep -q "FastQC_reports:" "$STATS_FILE"; then
        local fastqc_reports poor_quality
        fastqc_reports=$(grep "FastQC_reports:" "$STATS_FILE" | cut -d: -f2 | tr -d ' ')
        poor_quality=$(grep "Poor_quality_files:" "$STATS_FILE" | cut -d: -f2 | tr -d ' ' || echo "0")
        
        if [[ $poor_quality -gt 0 ]]; then
            warning "FastQC detected quality issues in $poor_quality files"
            warning "Check FastQC reports for details"
        else
            success "All FastQC reports passed quality checks"
        fi
    fi
    
    # Check if MultiQC report exists
    local multiqc_report="$OUTPUT_DIR/multiqc/${SAMPLE_ID}_multiqc_report.html"
    if [[ -f "$multiqc_report" ]]; then
        success "MultiQC report available: $multiqc_report"
    fi
}

# Cleanup QC temporary files
cleanup_qc() {
    if [[ -n "$TEMP_DIR" && -d "$TEMP_DIR" ]]; then
        # Keep statistics file
        if [[ -f "$STATS_FILE" ]]; then
            local basename_file
            basename_file=$(basename "$STATS_FILE")
            cp "$STATS_FILE" "$OUTPUT_DIR/${basename_file}"
        fi
        
        # Remove temporary directory
        rm -rf "$TEMP_DIR"
        log "Temporary QC files cleaned"
    fi
}

# Main QC function
run_qc() {
    # Parse arguments
    parse_qc_args "$@"
    
    # Set parameters
    set_qc_params
    
    # Validate parameters
    validate_qc_params
    
    log "Starting QC analysis"
    log "QC type: $QC_TYPE"
    log "Output directory: $OUTPUT_DIR"
    
    # Run QC based on type
    case "$QC_TYPE" in
        "fastq")
            run_fastqc_analysis
            ;;
        "bam")
            run_rseqc_analysis
            ;;
        "all")
            run_fastqc_analysis
            run_rseqc_analysis
            ;;
        *)
            error "Unknown QC type: $QC_TYPE"
            exit 1
            ;;
    esac
    
    # Run MultiQC
    run_multiqc_analysis
    
    # Quality check
    check_qc_quality
    
    # Cleanup
    cleanup_qc
    
    success "QC module completed successfully"
    return 0
}

# Usage function for QC module
qc_usage() {
    cat << EOF
Quality Control Module Usage:

$0 --qc-type TYPE --output-dir DIR --sample-id ID [options]

Required arguments:
  --qc-type TYPE        QC analysis type (fastq|bam|all)
  --output-dir DIR      Output directory for QC results
  --sample-id ID        Sample identifier

FASTQ QC arguments:
  --fastq-files FILES   Comma-separated list of FASTQ files

BAM QC arguments:
  --bam-file FILE       Input BAM file
  --genome-gtf FILE     Gene annotation GTF file (for RSeQC)
  --bed-file FILE       BED file (optional, will convert from GTF if not provided)

Optional arguments:
  --single-end          Single-end sequencing mode
  --threads N           Number of threads [default: auto-detect]
  --memory N            Available memory in GB [default: auto-detect]
  --temp-dir DIR        Temporary directory [default: auto-create]
  --fastqc-kmers N      K-mer size for FastQC [default: 7]
  --skip-rseqc          Skip RSeQC analysis
  --skip-fastqc         Skip FastQC analysis
  --skip-multiqc        Skip MultiQC analysis
  --verbose             Verbose output

Examples:
  # FASTQ QC only
  $0 --qc-type fastq --fastq-files sample_1.fq.gz,sample_2.fq.gz \\
     --output-dir qc_results/ --sample-id sample

  # BAM QC only
  $0 --qc-type bam --bam-file sample.bam --genome-gtf genes.gtf \\
     --output-dir qc_results/ --sample-id sample

  # Complete QC analysis
  $0 --qc-type all --fastq-files sample_1.fq.gz,sample_2.fq.gz \\
     --bam-file sample.bam --genome-gtf genes.gtf \\
     --output-dir qc_results/ --sample-id sample

  # Low memory system (skip RSeQC)
  $0 --qc-type fastq --fastq-files sample_1.fq.gz,sample_2.fq.gz \\
     --output-dir qc_results/ --sample-id sample --skip-rseqc
EOF
}

# Run QC if script is called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    if [[ $# -eq 0 ]]; then
        qc_usage
        exit 1
    fi
    
    run_qc "$@"
fi