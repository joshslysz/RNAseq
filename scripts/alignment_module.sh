#!/bin/bash

# Alignment Module
# Author: Joshua Slysz, PhD - Dalhousie University
# Memory-optimized alignment with Bowtie2 and minimap2 support
# Dynamically selects optimal aligner based on available resources

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
log() { echo -e "${BLUE}[ALIGN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }

# Default parameters
ALIGNER=""
GENOME_FASTA=""
BOWTIE2_INDEX=""
THREADS=0
MEMORY_GB=0
READ1=""
READ2=""
OUTPUT_BAM=""
SAMPLE_ID=""
TEMP_DIR=""
SINGLE_END=false
VERBOSE=false

# Alignment statistics
STATS_FILE=""

# Parse alignment module arguments
parse_alignment_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --aligner)
                ALIGNER="$2"
                shift 2
                ;;
            --genome-fasta)
                GENOME_FASTA="$2"
                shift 2
                ;;
            --bowtie2-index)
                BOWTIE2_INDEX="$2"
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
            --read1)
                READ1="$2"
                shift 2
                ;;
            --read2)
                READ2="$2"
                shift 2
                ;;
            --output)
                OUTPUT_BAM="$2"
                shift 2
                ;;
            --sample-id)
                SAMPLE_ID="$2"
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
            --verbose)
                VERBOSE=true
                shift
                ;;
            *)
                error "Unknown alignment option: $1"
                exit 1
                ;;
        esac
    done
}

# Validate alignment parameters
validate_alignment_params() {
    local errors=0
    
    if [[ -z "$READ1" ]]; then
        error "Read 1 file must be specified with --read1"
        errors=$((errors + 1))
    elif [[ ! -f "$READ1" ]]; then
        error "Read 1 file not found: $READ1"
        errors=$((errors + 1))
    fi
    
    if [[ "$SINGLE_END" == "false" ]]; then
        if [[ -z "$READ2" ]]; then
            error "Read 2 file must be specified with --read2 for paired-end data"
            errors=$((errors + 1))
        elif [[ ! -f "$READ2" ]]; then
            error "Read 2 file not found: $READ2"
            errors=$((errors + 1))
        fi
    fi
    
    if [[ -z "$OUTPUT_BAM" ]]; then
        error "Output BAM file must be specified with --output"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$SAMPLE_ID" ]]; then
        error "Sample ID must be specified with --sample-id"
        errors=$((errors + 1))
    fi
    
    # Check aligner-specific requirements
    if [[ "$ALIGNER" == "bowtie2" || "$ALIGNER" == "bowtie2-lowmem" ]]; then
        if [[ -z "$BOWTIE2_INDEX" ]]; then
            error "Bowtie2 index must be specified with --bowtie2-index"
            errors=$((errors + 1))
        elif [[ ! -f "${BOWTIE2_INDEX}.1.bt2" ]]; then
            error "Bowtie2 index not found: ${BOWTIE2_INDEX}.1.bt2"
            errors=$((errors + 1))
        fi
    elif [[ "$ALIGNER" == "minimap2" ]]; then
        if [[ -z "$GENOME_FASTA" ]]; then
            error "Genome FASTA must be specified with --genome-fasta for minimap2"
            errors=$((errors + 1))
        elif [[ ! -f "$GENOME_FASTA" ]]; then
            error "Genome FASTA not found: $GENOME_FASTA"
            errors=$((errors + 1))
        fi
    fi
    
    if [[ $errors -gt 0 ]]; then
        error "Alignment parameter validation failed"
        exit 1
    fi
}

# Auto-detect optimal aligner
auto_detect_aligner() {
    if [[ -n "$ALIGNER" && "$ALIGNER" != "auto" ]]; then
        return 0  # Aligner already specified
    fi
    
    local available_mem
    available_mem=$(get_available_memory)
    
    log "Auto-detecting optimal aligner (Available memory: ${available_mem}GB)"
    
    if [[ $available_mem -ge 8 ]]; then
        ALIGNER="bowtie2"
        log "Selected: Bowtie2 (standard mode) - sufficient memory available"
    elif [[ $available_mem -ge 4 ]]; then
        ALIGNER="bowtie2-lowmem"
        log "Selected: Bowtie2 (low memory mode) - moderate memory available"
    else
        ALIGNER="minimap2"
        log "Selected: minimap2 - low memory system detected"
    fi
}

# Set alignment parameters based on system resources
set_alignment_params() {
    # Auto-detect system resources if not specified
    if [[ $THREADS -eq 0 ]]; then
        THREADS=$(get_optimal_threads "alignment")
    fi
    
    if [[ $MEMORY_GB -eq 0 ]]; then
        MEMORY_GB=$(get_available_memory)
    fi
    
    # Create temporary directory if not specified
    if [[ -z "$TEMP_DIR" ]]; then
        TEMP_DIR="/tmp/alignment_$$"
        mkdir -p "$TEMP_DIR"
    fi
    
    # Set statistics file
    STATS_FILE="${TEMP_DIR}/${SAMPLE_ID}_alignment_stats.txt"
    
    log "Alignment parameters:"
    log "  Aligner: $ALIGNER"
    log "  Threads: $THREADS"
    log "  Memory: ${MEMORY_GB}GB"
    log "  Sample: $SAMPLE_ID"
    log "  Single-end: $SINGLE_END"
}

# Run Bowtie2 alignment
run_bowtie2_alignment() {
    log "Running Bowtie2 alignment for sample: $SAMPLE_ID"
    
    # Check memory requirements
    local required_mem=4
    if [[ "$ALIGNER" == "bowtie2-lowmem" ]]; then
        required_mem=3
    fi
    
    if ! check_memory_requirement $required_mem; then
        warning "Insufficient memory for Bowtie2, switching to minimap2"
        ALIGNER="minimap2"
        run_minimap2_alignment
        return
    fi
    
    # Set Bowtie2 parameters
    local bt2_opts=""
    local log_file="${TEMP_DIR}/${SAMPLE_ID}_bowtie2.log"
    
    # Memory optimization options
    if [[ "$ALIGNER" == "bowtie2-lowmem" ]]; then
        bt2_opts="--mm"  # Memory-mapped I/O
        log "Using Bowtie2 low-memory mode"
    fi
    
    # Additional Bowtie2 options for RNA-seq
    bt2_opts="$bt2_opts --sensitive --no-discordant --no-mixed"
    
    # Progress tracking
    echo "Starting Bowtie2 alignment at $(date)" > "$log_file"
    
    # Run alignment
    local alignment_cmd
    if [[ "$SINGLE_END" == "true" ]]; then
        alignment_cmd="bowtie2 -x $BOWTIE2_INDEX -U $READ1 -p $THREADS $bt2_opts"
    else
        alignment_cmd="bowtie2 -x $BOWTIE2_INDEX -1 $READ1 -2 $READ2 -p $THREADS $bt2_opts"
    fi
    
    log "Executing: $alignment_cmd"
    
    # Run alignment and sort
    if eval "$alignment_cmd" 2>>"$log_file" | \
       samtools sort -@ "$THREADS" -m 1G -o "$OUTPUT_BAM" -; then
        success "Bowtie2 alignment completed successfully"
    else
        error "Bowtie2 alignment failed"
        if [[ -f "$log_file" ]]; then
            error "Check log file: $log_file"
            tail -10 "$log_file"
        fi
        exit 1
    fi
    
    # Extract alignment statistics
    extract_bowtie2_stats "$log_file"
}

# Run minimap2 alignment
run_minimap2_alignment() {
    log "Running minimap2 alignment for sample: $SAMPLE_ID"
    
    # Check memory requirements (minimap2 is very memory efficient)
    if ! check_memory_requirement 1; then
        error "Insufficient memory even for minimap2"
        exit 1
    fi
    
    local log_file="${TEMP_DIR}/${SAMPLE_ID}_minimap2.log"
    
    # minimap2 options for RNA-seq
    local mm2_opts="-ax sr"  # Short read mode
    
    # Progress tracking
    echo "Starting minimap2 alignment at $(date)" > "$log_file"
    
    # Run alignment
    local alignment_cmd
    if [[ "$SINGLE_END" == "true" ]]; then
        alignment_cmd="minimap2 $mm2_opts -t $THREADS $GENOME_FASTA $READ1"
    else
        alignment_cmd="minimap2 $mm2_opts -t $THREADS $GENOME_FASTA $READ1 $READ2"
    fi
    
    log "Executing: $alignment_cmd"
    
    # Run alignment and sort
    if eval "$alignment_cmd" 2>>"$log_file" | \
       samtools sort -@ "$THREADS" -m 1G -o "$OUTPUT_BAM" -; then
        success "minimap2 alignment completed successfully"
    else
        error "minimap2 alignment failed"
        if [[ -f "$log_file" ]]; then
            error "Check log file: $log_file"
            tail -10 "$log_file"
        fi
        exit 1
    fi
    
    # Extract alignment statistics
    extract_minimap2_stats "$log_file"
}

# Extract Bowtie2 alignment statistics
extract_bowtie2_stats() {
    local log_file=$1
    
    if [[ ! -f "$log_file" ]]; then
        warning "Bowtie2 log file not found: $log_file"
        return
    fi
    
    # Parse Bowtie2 output
    local total_reads overall_alignment_rate
    total_reads=$(grep "reads; of these:" "$log_file" | awk '{print $1}')
    overall_alignment_rate=$(grep "overall alignment rate" "$log_file" | awk '{print $1}')
    
    # Write statistics
    cat > "$STATS_FILE" << EOF
Sample: $SAMPLE_ID
Aligner: Bowtie2
Total_reads: ${total_reads:-"Unknown"}
Overall_alignment_rate: ${overall_alignment_rate:-"Unknown"}
Log_file: $log_file
EOF
    
    log "Alignment rate: ${overall_alignment_rate:-"Unknown"}"
}

# Extract minimap2 alignment statistics
extract_minimap2_stats() {
    local log_file=$1
    
    if [[ ! -f "$log_file" ]]; then
        warning "minimap2 log file not found: $log_file"
        return
    fi
    
    # minimap2 doesn't provide detailed stats in stderr like Bowtie2
    # We'll use samtools flagstat for basic statistics
    
    cat > "$STATS_FILE" << EOF
Sample: $SAMPLE_ID
Aligner: minimap2
Log_file: $log_file
Note: Use samtools flagstat for detailed alignment statistics
EOF
    
    log "minimap2 alignment completed (use samtools flagstat for detailed stats)"
}

# Generate BAM file index and statistics
post_alignment_processing() {
    log "Post-alignment processing for sample: $SAMPLE_ID"
    
    # Index BAM file
    log "Indexing BAM file..."
    if samtools index "$OUTPUT_BAM"; then
        success "BAM file indexed successfully"
    else
        error "Failed to index BAM file"
        exit 1
    fi
    
    # Generate flagstat statistics
    local flagstat_file="${OUTPUT_BAM}.flagstat"
    log "Generating alignment statistics..."
    if samtools flagstat "$OUTPUT_BAM" > "$flagstat_file"; then
        success "Alignment statistics generated: $flagstat_file"
    else
        error "Failed to generate alignment statistics"
        exit 1
    fi
    
    # Extract key statistics for summary
    local total_reads mapped_reads mapping_rate
    total_reads=$(grep "in total" "$flagstat_file" | awk '{print $1}')
    mapped_reads=$(grep "mapped (" "$flagstat_file" | head -1 | awk '{print $1}')
    
    if [[ -n "$total_reads" && "$total_reads" -gt 0 ]]; then
        mapping_rate=$(echo "scale=2; $mapped_reads * 100 / $total_reads" | bc -l 2>/dev/null || echo "N/A")
        log "Mapping summary: $mapped_reads/$total_reads reads mapped (${mapping_rate}%)"
    fi
    
    # Update statistics file
    cat >> "$STATS_FILE" << EOF
Total_reads_flagstat: ${total_reads:-"Unknown"}
Mapped_reads: ${mapped_reads:-"Unknown"}
Mapping_rate_percent: ${mapping_rate:-"Unknown"}
BAM_file: $OUTPUT_BAM
Index_file: ${OUTPUT_BAM}.bai
Flagstat_file: $flagstat_file
EOF
}

# Check alignment quality
check_alignment_quality() {
    local flagstat_file="${OUTPUT_BAM}.flagstat"
    
    if [[ ! -f "$flagstat_file" ]]; then
        warning "Cannot check alignment quality - flagstat file not found"
        return
    fi
    
    # Extract mapping rate
    local total_reads mapped_reads
    total_reads=$(grep "in total" "$flagstat_file" | awk '{print $1}')
    mapped_reads=$(grep "mapped (" "$flagstat_file" | head -1 | awk '{print $1}')
    
    if [[ -n "$total_reads" && "$total_reads" -gt 0 ]]; then
        local mapping_rate
        mapping_rate=$(echo "scale=0; $mapped_reads * 100 / $total_reads" | bc -l 2>/dev/null || echo "0")
        
        if [[ "${mapping_rate%.*}" -lt 50 ]]; then
            error "Low mapping rate detected: ${mapping_rate}%"
            error "Consider checking:"
            error "  - Read quality"
            error "  - Reference genome compatibility"
            error "  - Adapter contamination"
        elif [[ "${mapping_rate%.*}" -lt 70 ]]; then
            warning "Moderate mapping rate: ${mapping_rate}%"
            warning "This may be acceptable depending on your experimental conditions"
        else
            success "Good mapping rate: ${mapping_rate}%"
        fi
    fi
}

# Cleanup alignment temporary files
cleanup_alignment() {
    if [[ -n "$TEMP_DIR" && -d "$TEMP_DIR" ]]; then
        # Keep statistics file and important logs
        local keep_files=("$STATS_FILE")
        
        for file in "${keep_files[@]}"; do
            if [[ -f "$file" ]]; then
                local basename_file
                basename_file=$(basename "$file")
                cp "$file" "$(dirname "$OUTPUT_BAM")/${basename_file}"
            fi
        done
        
        # Remove temporary directory
        rm -rf "$TEMP_DIR"
        log "Temporary alignment files cleaned"
    fi
}

# Main alignment function
run_alignment() {
    # Parse arguments
    parse_alignment_args "$@"
    
    # Auto-detect aligner if needed
    auto_detect_aligner
    
    # Set parameters
    set_alignment_params
    
    # Validate parameters
    validate_alignment_params
    
    log "Starting alignment process"
    log "Input: $READ1 $([ "$SINGLE_END" == "false" ] && echo "$READ2")"
    log "Output: $OUTPUT_BAM"
    
    # Monitor memory usage during alignment
    local memory_log="${TEMP_DIR}/${SAMPLE_ID}_memory.log"
    
    # Run alignment based on selected aligner
    case "$ALIGNER" in
        "bowtie2"|"bowtie2-lowmem")
            run_bowtie2_alignment
            ;;
        "minimap2")
            run_minimap2_alignment
            ;;
        *)
            error "Unknown aligner: $ALIGNER"
            exit 1
            ;;
    esac
    
    # Post-processing
    post_alignment_processing
    
    # Quality check
    check_alignment_quality
    
    # Cleanup
    cleanup_alignment
    
    success "Alignment module completed successfully"
    return 0
}

# Usage function for alignment module
alignment_usage() {
    cat << EOF
Alignment Module Usage:

$0 --read1 FILE --output BAM --sample-id ID [options]

Required arguments:
  --read1 FILE          Input FASTQ file (read 1)
  --output BAM          Output BAM file path
  --sample-id ID        Sample identifier

Optional arguments:
  --read2 FILE          Input FASTQ file (read 2, for paired-end)
  --single-end          Single-end sequencing mode
  --aligner TOOL        Aligner to use (bowtie2|minimap2|auto) [default: auto]
  --bowtie2-index PREFIX  Bowtie2 index prefix
  --genome-fasta FILE   Genome FASTA file (for minimap2)
  --threads N           Number of threads [default: auto-detect]
  --memory N            Available memory in GB [default: auto-detect]
  --temp-dir DIR        Temporary directory [default: auto-create]
  --verbose             Verbose output

Examples:
  # Auto-detect aligner (recommended)
  $0 --read1 sample_1.fq.gz --read2 sample_2.fq.gz --output sample.bam --sample-id sample

  # Force Bowtie2 with specific index
  $0 --read1 sample_1.fq.gz --read2 sample_2.fq.gz --output sample.bam \\
     --sample-id sample --aligner bowtie2 --bowtie2-index /path/to/index

  # Force minimap2 for low memory systems
  $0 --read1 sample_1.fq.gz --read2 sample_2.fq.gz --output sample.bam \\
     --sample-id sample --aligner minimap2 --genome-fasta genome.fa
EOF
}

# Run alignment if script is called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    if [[ $# -eq 0 ]]; then
        alignment_usage
        exit 1
    fi
    
    run_alignment "$@"
fi