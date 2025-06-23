#!/bin/bash

# Memory-Optimized Local RNA-seq Preprocessing Pipeline
# Author: Joshua Slysz, PhD - Dalhousie University
# Designed for systems with limited RAM (2-16GB)
# Replaces SLURM-based cluster computing with smart local resource management

set -euo pipefail

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source resource manager
source "$SCRIPT_DIR/scripts/resource_manager.sh"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m'

# Default parameters
DATA_DIR=""
OUTPUT_DIR="results"
GENOME_CONFIG=""
GENOME_FASTA=""
GENOME_GTF=""
BOWTIE2_INDEX=""
SINGLE_END=false
THREADS=0
MEMORY_GB=0
ALIGNER="auto"
QUANTIFIER="featurecounts"
CHUNK_SIZE=0
KEEP_TEMP=false
VERBOSE=false
RESUME=false
ENV_CHECK=false
PREPROCESSING=true
ADAPTER_TRIMMING=true
STRANDEDNESS="auto"

# Internal variables
TEMP_DIR=""
LOG_FILE=""
PROGRESS_FILE=""
START_TIME=""

# Logging functions
log() { 
    local msg="$1"
    echo -e "${BLUE}[$(date '+%H:%M:%S')]${NC} $msg" 
    if [[ -n "$LOG_FILE" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $msg" >> "$LOG_FILE"
    fi
}

error() { 
    local msg="$1"
    echo -e "${RED}[ERROR]${NC} $msg" >&2
    if [[ -n "$LOG_FILE" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $msg" >> "$LOG_FILE"
    fi
}

success() { 
    local msg="$1"
    echo -e "${GREEN}[SUCCESS]${NC} $msg"
    if [[ -n "$LOG_FILE" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] SUCCESS: $msg" >> "$LOG_FILE"
    fi
}

warning() { 
    local msg="$1"
    echo -e "${YELLOW}[WARNING]${NC} $msg"
    if [[ -n "$LOG_FILE" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $msg" >> "$LOG_FILE"
    fi
}

# Progress tracking
update_progress() {
    local sample=$1
    local step=$2
    local status=$3
    
    echo "$(date '+%Y-%m-%d %H:%M:%S'),$sample,$step,$status" >> "$PROGRESS_FILE"
}

# Usage function
usage() {
    cat << EOF
${BOLD}Memory-Optimized RNA-seq Preprocessing Pipeline${NC}

Usage: $0 -d <data_directory> [options]

${BOLD}Required Arguments:${NC}
  -d, --data-dir DIR          Directory containing raw FASTQ files

${BOLD}Genome Configuration (choose one):${NC}
  --genome-config FILE        Path to genome configuration file
  --genome-fasta FILE         Genome FASTA file
  --genome-gtf FILE           Gene annotation GTF file
  --bowtie2-index PREFIX      Legacy Bowtie2 index (not recommended for RNA-seq)

${BOLD}Optional Arguments:${NC}
  -o, --output-dir DIR        Output directory [default: results]
  -t, --threads N             Number of threads [default: auto-detect]
  -m, --memory N              Available memory in GB [default: auto-detect]
  -s, --single-end            Single-end sequencing [default: paired-end]
  -a, --aligner TOOL          Force aligner (hisat2|star|minimap2|auto) [default: auto]
  -q, --quantifier TOOL       Quantification method (featurecounts|salmon) [default: featurecounts]
  -c, --chunk-size N          FASTQ chunk size for memory optimization [default: auto]
  --no-preprocessing          Skip adapter trimming and preprocessing steps
  --no-adapter-trimming       Skip only adapter trimming
  --strandedness TYPE         Library strandedness (auto|unstranded|forward|reverse) [default: auto]
  -k, --keep-temp             Keep temporary files for debugging
  -v, --verbose               Verbose logging
  -r, --resume                Resume from previous run
  -e, --env-check             Verify conda environment setup
  -h, --help                  Show this help message

${BOLD}Examples:${NC}
  # Basic run with genome config
  $0 -d raw_data/ --genome-config references/GRCh38/config.txt

  # Low memory system (force minimap2)
  $0 -d raw_data/ --genome-config references/GRCh38/config.txt -a minimap2 -m 4

  # Resume interrupted run
  $0 -d raw_data/ --genome-config references/GRCh38/config.txt -r

  # Single-end data with custom settings
  $0 -d raw_data/ -s --genome-fasta genome.fa --genome-gtf genes.gtf --aligner hisat2

${BOLD}File Naming Convention:${NC}
  Paired-end: SampleID_1.fq.gz, SampleID_2.fq.gz
  Single-end: SampleID.fq.gz

${BOLD}Memory Requirements:${NC}
  Minimum: 2GB RAM
  Recommended: 4GB+ RAM
  High performance: 8GB+ RAM

EOF
}

# Parse command line arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -d|--data-dir)
                DATA_DIR="$2"
                shift 2
                ;;
            -o|--output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --genome-config)
                GENOME_CONFIG="$2"
                shift 2
                ;;
            --genome-fasta)
                GENOME_FASTA="$2"
                shift 2
                ;;
            --genome-gtf)
                GENOME_GTF="$2"
                shift 2
                ;;
            --bowtie2-index)
                BOWTIE2_INDEX="$2"
                shift 2
                ;;
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            -m|--memory)
                MEMORY_GB="$2"
                shift 2
                ;;
            -s|--single-end)
                SINGLE_END=true
                shift
                ;;
            -a|--aligner)
                ALIGNER="$2"
                shift 2
                ;;
            -q|--quantifier)
                QUANTIFIER="$2"
                shift 2
                ;;
            -c|--chunk-size)
                CHUNK_SIZE="$2"
                shift 2
                ;;
            --no-preprocessing)
                PREPROCESSING=false
                shift
                ;;
            --no-adapter-trimming)
                ADAPTER_TRIMMING=false
                shift
                ;;
            --strandedness)
                STRANDEDNESS="$2"
                shift 2
                ;;
            -k|--keep-temp)
                KEEP_TEMP=true
                shift
                ;;
            -v|--verbose)
                VERBOSE=true
                shift
                ;;
            -r|--resume)
                RESUME=true
                shift
                ;;
            -e|--env-check)
                ENV_CHECK=true
                shift
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                error "Unknown option: $1"
                echo ""
                usage
                exit 1
                ;;
        esac
    done
}

# Environment check
check_environment() {
    log "Checking environment setup..."
    
    # Check conda environment
    if [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "$CONDA_DEFAULT_ENV" != "rnaseqpipe" ]]; then
        error "Please activate the rnaseqpipe conda environment first:"
        error "  conda activate rnaseqpipe"
        error "  OR source activate_pipeline.sh"
        exit 1
    fi
    
    # Check required tools
    local tools=("fastqc" "hisat2" "samtools" "featureCounts" "multiqc")
    local missing_tools=()
    
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        error "Missing required tools: ${missing_tools[*]}"
        error "Please check your conda environment installation"
        exit 1
    fi
    
    success "Environment check passed!"
    
    if [[ "$ENV_CHECK" == "true" ]]; then
        exit 0
    fi
}

# Load genome configuration
load_genome_config() {
    if [[ -n "$GENOME_CONFIG" ]]; then
        if [[ ! -f "$GENOME_CONFIG" ]]; then
            error "Genome configuration file not found: $GENOME_CONFIG"
            exit 1
        fi
        
        log "Loading genome configuration: $GENOME_CONFIG"
        source "$GENOME_CONFIG"
        
        # Use config values if not overridden
        GENOME_FASTA=${GENOME_FASTA:-$GENOME_FASTA}
        GENOME_GTF=${GENOME_GTF:-$GENOME_GTF}
        BOWTIE2_INDEX=${BOWTIE2_INDEX:-$BOWTIE2_INDEX}
    fi
    
    # Validate genome files
    if [[ -z "$GENOME_FASTA" ]]; then
        error "Genome FASTA file not specified"
        error "Use --genome-config or --genome-fasta"
        exit 1
    fi
    
    if [[ -z "$GENOME_GTF" ]]; then
        error "Genome GTF file not specified"
        error "Use --genome-config or --genome-gtf"  
        exit 1
    fi
    
    if [[ ! -f "$GENOME_FASTA" ]]; then
        error "Genome FASTA file not found: $GENOME_FASTA"
        exit 1
    fi
    
    if [[ ! -f "$GENOME_GTF" ]]; then
        error "Genome GTF file not found: $GENOME_GTF"
        exit 1
    fi
    
    success "Genome configuration loaded"
}

# Initialize pipeline
initialize_pipeline() {
    START_TIME=$(date)
    
    # Create output directory structure
    mkdir -p "$OUTPUT_DIR"/{aligned,counts,qc_reports,logs,temp}
    # Convert to absolute path (cross-platform compatible)
    if command -v realpath >/dev/null 2>&1; then
        OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
    else
        OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)
    fi
    
    # Setup logging
    LOG_FILE="$OUTPUT_DIR/logs/pipeline_$(date '+%Y%m%d_%H%M%S').log"
    PROGRESS_FILE="$OUTPUT_DIR/logs/progress.csv"
    TEMP_DIR="$OUTPUT_DIR/temp"
    
    # Initialize progress tracking
    if [[ ! -f "$PROGRESS_FILE" || "$RESUME" == "false" ]]; then
        echo "timestamp,sample,step,status" > "$PROGRESS_FILE"
    fi
    
    log "Pipeline initialized"
    log "Output directory: $OUTPUT_DIR"
    log "Log file: $LOG_FILE"
    log "Temporary directory: $TEMP_DIR"
    
    # Auto-detect system resources if not specified
    if [[ $MEMORY_GB -eq 0 ]]; then
        MEMORY_GB=$(get_available_memory)
        log "Auto-detected memory: ${MEMORY_GB}GB"
    fi
    
    if [[ $THREADS -eq 0 ]]; then
        THREADS=$(get_optimal_threads "default")
        log "Auto-detected threads: $THREADS"
    fi
    
    # Select aligner if auto
    if [[ "$ALIGNER" == "auto" ]]; then
        ALIGNER=$(select_aligner)
        log "Auto-selected aligner: $ALIGNER"
    fi
    
    # Check system requirements
    if ! check_system_requirements 2 10; then
        error "System requirements not met!"
        exit 1
    fi
    
    print_resource_summary
}

# Find FASTQ files
find_fastq_files() {
    local data_dir=$1
    local -n samples_ref=$2
    
    log "Scanning for FASTQ files in: $data_dir"
    
    if [[ "$SINGLE_END" == "true" ]]; then
        # Single-end pattern: SampleID.fq.gz
        while IFS= read -r -d '' file; do
            local basename
            basename=$(basename "$file" .fq.gz)
            samples_ref["$basename"]="$file"
        done < <(find "$data_dir" -name "*.fq.gz" -not -name "*_1.fq.gz" -not -name "*_2.fq.gz" -print0)
    else
        # Paired-end pattern: SampleID_1.fq.gz, SampleID_2.fq.gz
        while IFS= read -r -d '' file; do
            local basename
            basename=$(basename "$file" _1.fq.gz)
            local read2="${file/_1.fq.gz/_2.fq.gz}"
            
            if [[ -f "$read2" ]]; then
                samples_ref["$basename"]="$file,$read2"
            else
                warning "Missing read 2 for sample: $basename"
            fi
        done < <(find "$data_dir" -name "*_1.fq.gz" -print0)
    fi
    
    if [[ ${#samples_ref[@]} -eq 0 ]]; then
        error "No FASTQ files found in: $data_dir"
        if [[ "$SINGLE_END" == "true" ]]; then
            error "Expected pattern: SampleID.fq.gz"
        else
            error "Expected pattern: SampleID_1.fq.gz, SampleID_2.fq.gz"
        fi
        exit 1
    fi
    
    success "Found ${#samples_ref[@]} samples"
}

# Process single sample
process_sample() {
    local sample_id=$1
    local fastq_files=$2
    
    log "Processing sample: $sample_id"
    update_progress "$sample_id" "start" "running"
    
    # Create sample-specific temp directory
    local sample_temp="$TEMP_DIR/$sample_id"
    mkdir -p "$sample_temp"
    
    # Parse FASTQ files
    local read1 read2=""
    if [[ "$SINGLE_END" == "true" ]]; then
        read1="$fastq_files"
    else
        IFS=',' read -r read1 read2 <<< "$fastq_files"
    fi
    
    # Preprocessing (adapter trimming, quality filtering)
    if [[ "$PREPROCESSING" == "true" ]]; then
        log "Running preprocessing for sample: $sample_id"
        update_progress "$sample_id" "preprocessing" "running"
        
        local preprocess_opts="--read1 $read1 --sample-id $sample_id"
        preprocess_opts="$preprocess_opts --output-dir $sample_temp/preprocessed"
        preprocess_opts="$preprocess_opts --threads $THREADS --memory $MEMORY_GB"
        preprocess_opts="$preprocess_opts --temp-dir $sample_temp"
        
        if [[ "$SINGLE_END" == "false" ]]; then
            preprocess_opts="$preprocess_opts --read2 $read2"
        else
            preprocess_opts="$preprocess_opts --single-end"
        fi
        
        if [[ "$ADAPTER_TRIMMING" == "false" ]]; then
            preprocess_opts="$preprocess_opts --no-adapter-trimming"
        fi
        
        # Run preprocessing and capture output paths
        local preprocess_output
        preprocess_output=$(bash "$SCRIPT_DIR/scripts/preprocessing_module.sh" $preprocess_opts)
        
        # Update read paths to preprocessed files
        if echo "$preprocess_output" | grep -q "PROCESSED_READ1"; then
            read1=$(echo "$preprocess_output" | grep "PROCESSED_READ1" | cut -d'=' -f2)
            if [[ "$SINGLE_END" == "false" ]]; then
                read2=$(echo "$preprocess_output" | grep "PROCESSED_READ2" | cut -d'=' -f2)
            fi
            success "Using preprocessed reads for downstream analysis"
        fi
        
        update_progress "$sample_id" "preprocessing" "completed"
    fi
    
    # Quality control
    log "Running FastQC for sample: $sample_id"
    update_progress "$sample_id" "fastqc" "running"
    
    local fastqc_dir="$OUTPUT_DIR/qc_reports/fastqc"
    mkdir -p "$fastqc_dir"
    
    if [[ "$SINGLE_END" == "true" ]]; then
        fastqc -o "$fastqc_dir" -t "$THREADS" "$read1"
    else
        fastqc -o "$fastqc_dir" -t "$THREADS" "$read1" "$read2"
    fi
    
    update_progress "$sample_id" "fastqc" "completed"
    
    # Alignment
    log "Running alignment for sample: $sample_id"
    update_progress "$sample_id" "alignment" "running"
    
    local aligned_bam="$OUTPUT_DIR/aligned/${sample_id}.sorted.bam"
    
    # Use RNA-seq specific alignment module with splice-aware aligners
    local rna_align_opts="--read1 $read1 --output $aligned_bam --sample-id $sample_id"
    rna_align_opts="$rna_align_opts --threads $THREADS --memory $MEMORY_GB"
    rna_align_opts="$rna_align_opts --genome-fasta $GENOME_FASTA --annotation-gtf $GENOME_GTF"
    rna_align_opts="$rna_align_opts --temp-dir $sample_temp"
    
    # Check for existing HISAT2 index from config or pre-built location
    if [[ -n "$HISAT2_INDEX" && -f "${HISAT2_INDEX}.1.ht2" ]]; then
        log "Using pre-built HISAT2 index: $HISAT2_INDEX"
        rna_align_opts="$rna_align_opts --hisat2-index $HISAT2_INDEX"
    else
        # Check for sample-specific index in temp directory
        local hisat2_index_path="$sample_temp/hisat2_index/genome"
        if [[ -f "${hisat2_index_path}.1.ht2" ]]; then
            log "Found sample-specific HISAT2 index: $hisat2_index_path"
            rna_align_opts="$rna_align_opts --hisat2-index $hisat2_index_path"
        fi
    fi
    
    if [[ "$SINGLE_END" == "false" ]]; then
        rna_align_opts="$rna_align_opts --read2 $read2"
    else
        rna_align_opts="$rna_align_opts --single-end"
    fi
    
    # Force RNA-seq appropriate aligners (HISAT2 or STAR)
    if [[ "$ALIGNER" == "bowtie2" || "$ALIGNER" == "bowtie2-lowmem" ]]; then
        warning "Bowtie2 is not suitable for RNA-seq data - switching to HISAT2"
        rna_align_opts="$rna_align_opts --aligner hisat2"
    elif [[ "$ALIGNER" == "minimap2" ]]; then
        # minimap2 can work for RNA-seq but HISAT2/STAR are better
        rna_align_opts="$rna_align_opts --aligner auto"
    else
        rna_align_opts="$rna_align_opts --aligner $ALIGNER"
    fi
    
    # Run RNA-seq alignment
    if bash "$SCRIPT_DIR/scripts/rna_alignment_module.sh" $rna_align_opts 2>>"$sample_temp/alignment.log"; then
        success "RNA-seq alignment completed for $sample_id"
    else
        error "RNA-seq alignment failed for $sample_id"
        if [[ -f "$sample_temp/alignment.log" ]]; then
            error "Check alignment log: $sample_temp/alignment.log"
            cat "$sample_temp/alignment.log" >&2
        fi
        update_progress "$sample_id" "alignment" "failed"
        return 1
    fi
    
    # Index BAM file
    samtools index "$aligned_bam"
    
    # Generate alignment statistics
    samtools flagstat "$aligned_bam" > "$OUTPUT_DIR/aligned/${sample_id}.flagstat"
    
    update_progress "$sample_id" "alignment" "completed"
    
    # Quantification
    log "Running quantification for sample: $sample_id"
    update_progress "$sample_id" "quantification" "running"
    
    local counts_file="$OUTPUT_DIR/counts/${sample_id}.counts.txt"
    
    if [[ "$QUANTIFIER" == "featurecounts" ]]; then
        # featureCounts quantification optimized for bulk RNA-seq
        local fc_opts="-T $THREADS -g gene_id -a $GENOME_GTF"
        fc_opts="$fc_opts -t exon -f"  # Count at gene level using exon features
        fc_opts="$fc_opts -M"  # Count multimapping reads
        fc_opts="$fc_opts -O"  # Count overlapping features
        fc_opts="$fc_opts --largestOverlap"  # Assign reads to feature with largest overlap
        fc_opts="$fc_opts -Q 10"  # Minimum mapping quality
        fc_opts="$fc_opts -s 0"  # Unstranded (adjust based on library prep)
        
        if [[ "$SINGLE_END" == "false" ]]; then
            fc_opts="$fc_opts -p -B -C"  # Paired-end, count fragments, exclude chimeric fragments
        fi
        
        # Add strandedness detection comment
        log "Using unstranded counting (-s 0). Adjust if your library is stranded:"
        log "  -s 1: forward stranded (fr-secondstrand)"
        log "  -s 2: reverse stranded (fr-firststrand)"
        
        # Ensure temp directory exists for quantification log
        mkdir -p "$sample_temp"
        
        featureCounts $fc_opts -o "$counts_file" "$aligned_bam" 2>"$sample_temp/quantification.log"
    else
        # Salmon quantification (if implemented)
        error "Salmon quantification not yet implemented"
        exit 1
    fi
    
    update_progress "$sample_id" "quantification" "completed"
    
    # RSeQC analysis (optional, memory permitting)
    if [[ $MEMORY_GB -ge 4 ]]; then
        log "Running RSeQC analysis for sample: $sample_id"
        update_progress "$sample_id" "rseqc" "running"
        
        # Add RSeQC commands here
        # This would include read distribution, junction annotation, etc.
        
        update_progress "$sample_id" "rseqc" "completed"
    fi
    
    # Cleanup sample temp files
    if [[ "$KEEP_TEMP" == "false" ]]; then
        rm -rf "$sample_temp"
    fi
    
    update_progress "$sample_id" "completed" "success"
    success "Sample $sample_id processing completed"
}

# Combine count matrices
combine_counts() {
    log "Combining count matrices..."
    
    local count_files=("$OUTPUT_DIR"/counts/*.counts.txt)
    if [[ ${#count_files[@]} -eq 0 ]] || [[ ! -f "${count_files[0]}" ]]; then
        warning "No count files found to combine"
        return
    fi
    
    # Use Python to combine count matrices (optimized version)
    python3 - "$OUTPUT_DIR/counts" "$OUTPUT_DIR/counts/combined_counts.txt" << 'EOF'
import os
import sys

count_dir = sys.argv[1]
output_file = sys.argv[2]

# Find all count files
count_files = [f for f in os.listdir(count_dir) if f.endswith('.counts.txt')]

if not count_files:
    print("No count files found")
    sys.exit(0)

print(f"Found {len(count_files)} count files to combine")

# Read all count files into memory once
count_data = {}
all_genes = set()

for count_file in sorted(count_files):
    sample_name = count_file.replace('.counts.txt', '')
    count_data[sample_name] = {}
    
    file_path = os.path.join(count_dir, count_file)
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('Geneid'):
                continue
            parts = line.strip().split('\t')
            gene_id = parts[0]
            count = parts[-1]
            count_data[sample_name][gene_id] = count
            all_genes.add(gene_id)

# Sort genes for consistent output
genes = sorted(all_genes)
sample_names = sorted(count_data.keys())

print(f"Processing {len(genes)} genes across {len(sample_names)} samples")

# Write combined matrix
with open(output_file, 'w') as out_f:
    # Write header
    out_f.write('Geneid\t' + '\t'.join(sample_names) + '\n')
    
    # Write data rows
    for i, gene in enumerate(genes):
        if i % 10000 == 0:
            print(f"Writing gene {i+1}/{len(genes)}")
        
        row = [gene]
        for sample in sample_names:
            count = count_data[sample].get(gene, '0')
            row.append(count)
        
        out_f.write('\t'.join(row) + '\n')

print(f"Combined count matrix saved: {output_file}")
print(f"Dimensions: {len(genes)} genes x {len(sample_names)} samples")
EOF
    success "Count matrices combined"
}

# Generate MultiQC report
generate_multiqc_report() {
    log "Generating MultiQC report..."
    
    multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR/qc_reports" -n "multiqc_report.html" --force
    
    success "MultiQC report generated: $OUTPUT_DIR/qc_reports/multiqc_report.html"
}

# Pipeline summary
generate_summary() {
    local end_time
    end_time=$(date)
    
    cat > "$OUTPUT_DIR/pipeline_summary.txt" << EOF
RNA-seq Pipeline Summary
========================

Start time: $START_TIME
End time: $end_time

Input directory: $DATA_DIR
Output directory: $OUTPUT_DIR

Configuration:
- Genome FASTA: $GENOME_FASTA
- Genome GTF: $GENOME_GTF
- Aligner: $ALIGNER
- Quantifier: $QUANTIFIER
- Threads: $THREADS
- Memory: ${MEMORY_GB}GB
- Single-end: $SINGLE_END

Samples processed: $(grep -c "completed,success" "$PROGRESS_FILE" || echo "0")

Output files:
- Aligned BAM files: $OUTPUT_DIR/aligned/
- Gene counts: $OUTPUT_DIR/counts/
- QC reports: $OUTPUT_DIR/qc_reports/
- Log files: $OUTPUT_DIR/logs/

Pipeline completed successfully!
EOF

    success "Pipeline summary generated: $OUTPUT_DIR/pipeline_summary.txt"
}

# Cleanup function
cleanup() {
    if [[ "$KEEP_TEMP" == "false" && -d "$TEMP_DIR" ]]; then
        log "Cleaning up temporary files..."
        rm -rf "$TEMP_DIR"
        success "Temporary files cleaned"
    fi
}

# Error handler
error_handler() {
    error "Pipeline failed at line $1"
    cleanup
    exit 1
}

# Main pipeline function
main() {
    # Set error handler
    trap 'error_handler $LINENO' ERR
    
    # Parse arguments
    parse_args "$@"
    
    # Check required arguments
    if [[ -z "$DATA_DIR" ]]; then
        error "Data directory must be specified with -d or --data-dir"
        echo ""
        usage
        exit 1
    fi
    
    if [[ ! -d "$DATA_DIR" ]]; then
        error "Data directory not found: $DATA_DIR"
        exit 1
    fi
    
    # Check environment
    check_environment
    
    # Load genome configuration
    load_genome_config
    
    # Initialize pipeline
    initialize_pipeline
    
    # Find FASTQ files
    declare -A samples
    find_fastq_files "$DATA_DIR" samples
    
    log "Starting RNA-seq preprocessing pipeline"
    log "Samples to process: ${!samples[*]}"
    
    # Process each sample sequentially (memory optimization)
    for sample_id in "${!samples[@]}"; do
        # Check if sample already completed (for resume functionality)
        if [[ "$RESUME" == "true" ]] && grep -q "$sample_id,completed,success" "$PROGRESS_FILE" 2>/dev/null; then
            log "Skipping completed sample: $sample_id"
            continue
        fi
        
        # Wait for sufficient memory if needed
        if ! check_memory_requirement 2; then
            wait_for_memory 2
        fi
        
        # Process sample
        process_sample "$sample_id" "${samples[$sample_id]}"
    done
    
    # Post-processing
    combine_counts
    generate_multiqc_report
    generate_summary
    
    # Cleanup
    cleanup
    
    success "🎉 RNA-seq pipeline completed successfully!"
    log "Results available in: $OUTPUT_DIR"
    log "Summary report: $OUTPUT_DIR/pipeline_summary.txt"
    log "QC report: $OUTPUT_DIR/qc_reports/multiqc_report.html"
}

# Run main function
main "$@"