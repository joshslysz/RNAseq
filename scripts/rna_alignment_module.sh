#!/bin/bash

# RNA-seq Specific Alignment Module
# Author: Joshua Slysz, PhD - Dalhousie University
# Splice-aware alignment with HISAT2, STAR, and troubleshooting support

set -euo pipefail

# Source resource manager
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/resource_manager.sh"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
NC='\033[0m'

# Logging functions
log() { echo -e "${BLUE}[RNA-ALIGN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }

# Default parameters
ALIGNER=""
GENOME_FASTA=""
ANNOTATION_GTF=""
HISAT2_INDEX=""
STAR_INDEX=""
THREADS=0
MEMORY_GB=0
READ1=""
READ2=""
OUTPUT_BAM=""
SAMPLE_ID=""
TEMP_DIR=""
SINGLE_END=false
STRANDNESS=""
TROUBLESHOOTING_MODE=false
VERBOSE=false

# Parse RNA alignment arguments
parse_rna_alignment_args() {
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
            --annotation-gtf)
                ANNOTATION_GTF="$2"
                shift 2
                ;;
            --hisat2-index)
                HISAT2_INDEX="$2"
                shift 2
                ;;
            --star-index)
                STAR_INDEX="$2"
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
            --strandness)
                STRANDNESS="$2"
                shift 2
                ;;
            --troubleshooting)
                TROUBLESHOOTING_MODE=true
                shift
                ;;
            --verbose)
                VERBOSE=true
                shift
                ;;
            *)
                error "Unknown RNA alignment option: $1"
                exit 1
                ;;
        esac
    done
}

# Auto-detect optimal RNA-seq aligner
auto_detect_rna_aligner() {
    if [[ -n "$ALIGNER" && "$ALIGNER" != "auto" ]]; then
        return 0  # Aligner already specified
    fi
    
    local available_mem
    available_mem=$(get_available_memory)
    
    log "Auto-detecting optimal RNA-seq aligner (Available memory: ${available_mem}GB)"
    
    # Check which aligners are available
    local hisat2_available=false
    local star_available=false
    
    if command -v hisat2 >/dev/null 2>&1; then
        hisat2_available=true
    fi
    
    if command -v STAR >/dev/null 2>&1; then
        star_available=true
    fi
    
    # Select based on availability and memory - prioritize memory-efficient aligners
    if [[ $available_mem -ge 32 && "$star_available" == "true" ]]; then
        ALIGNER="star"
        log "Selected: STAR - high memory system with STAR available"
    elif [[ $available_mem -le 8 ]] && command -v minimap2 >/dev/null 2>&1; then
        ALIGNER="minimap2"
        log "Selected: minimap2 - memory-efficient aligner for limited RAM systems"
    elif [[ "$hisat2_available" == "true" && $available_mem -ge 16 ]]; then
        ALIGNER="hisat2"
        log "Selected: HISAT2 - splice-aware aligner for RNA-seq"
    elif command -v minimap2 >/dev/null 2>&1; then
        ALIGNER="minimap2"
        warning "Selected: minimap2 - HISAT2 requires more memory than available"
    elif [[ "$star_available" == "true" ]]; then
        ALIGNER="star"
        warning "Selected: STAR with limited memory - may be slow"
    else
        error "No suitable RNA-seq aligner found"
        error "Please install minimap2 (recommended for limited memory) or HISAT2:"
        error "  conda install minimap2"
        error "  conda install hisat2"
        exit 1
    fi
}

# Validate RNA alignment parameters
validate_rna_alignment_params() {
    local errors=0
    
    # Basic file checks
    if [[ -z "$READ1" || ! -f "$READ1" ]]; then
        error "Read 1 file not found: $READ1"
        errors=$((errors + 1))
    fi
    
    if [[ "$SINGLE_END" == "false" && (-z "$READ2" || ! -f "$READ2") ]]; then
        error "Read 2 file not found: $READ2"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$OUTPUT_BAM" ]]; then
        error "Output BAM file must be specified"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$SAMPLE_ID" ]]; then
        error "Sample ID must be specified"
        errors=$((errors + 1))
    fi
    
    # Aligner-specific checks
    case "$ALIGNER" in
        "hisat2")
            if [[ -z "$HISAT2_INDEX" ]]; then
                if [[ -n "$GENOME_FASTA" ]]; then
                    log "HISAT2 index not provided, will build from FASTA"
                else
                    error "Either HISAT2 index or genome FASTA must be provided"
                    errors=$((errors + 1))
                fi
            elif [[ ! -f "${HISAT2_INDEX}.1.ht2" ]]; then
                error "HISAT2 index not found: ${HISAT2_INDEX}.1.ht2"
                errors=$((errors + 1))
            fi
            ;;
        "star")
            if [[ -z "$STAR_INDEX" ]]; then
                if [[ -n "$GENOME_FASTA" ]]; then
                    log "STAR index not provided, will build from FASTA"
                else
                    error "Either STAR index or genome FASTA must be provided"
                    errors=$((errors + 1))
                fi
            elif [[ ! -d "$STAR_INDEX" ]]; then
                error "STAR index directory not found: $STAR_INDEX"
                errors=$((errors + 1))
            fi
            ;;
        "minimap2")
            if [[ -z "$GENOME_FASTA" ]]; then
                error "Genome FASTA must be provided for minimap2"
                errors=$((errors + 1))
            elif [[ ! -f "$GENOME_FASTA" ]]; then
                error "Genome FASTA file not found: $GENOME_FASTA"
                errors=$((errors + 1))
            fi
            ;;
    esac
    
    if [[ $errors -gt 0 ]]; then
        exit 1
    fi
}

# Build HISAT2 index if needed
build_hisat2_index() {
    if [[ -n "$HISAT2_INDEX" && -f "${HISAT2_INDEX}.1.ht2" ]]; then
        return 0  # Index already exists
    fi
    
    if [[ -z "$GENOME_FASTA" ]]; then
        error "Cannot build HISAT2 index without genome FASTA"
        exit 1
    fi
    
    log "Building HISAT2 index from genome FASTA"
    
    # Set index path if not provided
    if [[ -z "$HISAT2_INDEX" ]]; then
        HISAT2_INDEX="${TEMP_DIR}/hisat2_index/genome"
        mkdir -p "$(dirname "$HISAT2_INDEX")"
    fi
    
    # Build index with annotation if available
    local build_cmd="hisat2-build"
    local build_opts=""
    
    if [[ -n "$ANNOTATION_GTF" && -f "$ANNOTATION_GTF" ]]; then
        log "Using gene annotation for splice site detection"
        
        # Extract splice sites and exons
        local ss_file="${TEMP_DIR}/splicesites.txt"
        local exon_file="${TEMP_DIR}/exons.txt"
        
        hisat2_extract_splice_sites.py "$ANNOTATION_GTF" > "$ss_file"
        hisat2_extract_exons.py "$ANNOTATION_GTF" > "$exon_file"
        
        build_opts="--ss $ss_file --exon $exon_file"
    fi
    
    # Run index building
    log "Building HISAT2 index: $build_cmd $build_opts $GENOME_FASTA $HISAT2_INDEX"
    if eval "$build_cmd $build_opts -p $THREADS $GENOME_FASTA $HISAT2_INDEX"; then
        success "HISAT2 index built successfully"
    else
        error "Failed to build HISAT2 index"
        exit 1
    fi
}

# Build STAR index if needed
build_star_index() {
    if [[ -n "$STAR_INDEX" && -d "$STAR_INDEX" && -f "$STAR_INDEX/genomeParameters.txt" ]]; then
        return 0  # Index already exists
    fi
    
    if [[ -z "$GENOME_FASTA" ]]; then
        error "Cannot build STAR index without genome FASTA"
        exit 1
    fi
    
    log "Building STAR index from genome FASTA"
    
    # Set index path if not provided
    if [[ -z "$STAR_INDEX" ]]; then
        STAR_INDEX="${TEMP_DIR}/star_index"
    fi
    mkdir -p "$STAR_INDEX"
    
    # Check memory requirements for STAR
    local required_mem=30
    if ! check_memory_requirement $required_mem; then
        warning "STAR typically requires 30GB+ RAM for human genome"
        warning "Proceeding with limited memory mode"
    fi
    
    # Build index with annotation if available
    local build_cmd="STAR --runMode genomeGenerate --genomeDir $STAR_INDEX"
    build_cmd="$build_cmd --genomeFastaFiles $GENOME_FASTA"
    build_cmd="$build_cmd --runThreadN $THREADS"
    
    if [[ -n "$ANNOTATION_GTF" && -f "$ANNOTATION_GTF" ]]; then
        log "Using gene annotation for splice junction database"
        build_cmd="$build_cmd --sjdbGTFfile $ANNOTATION_GTF --sjdbOverhang 99"
    fi
    
    # Add memory optimization for limited memory systems
    if [[ $MEMORY_GB -lt 32 ]]; then
        build_cmd="$build_cmd --limitGenomeGenerateRAM $((MEMORY_GB * 1000000000))"
        build_cmd="$build_cmd --genomeSAindexNbases 12"  # Reduce memory usage
    fi
    
    log "Building STAR index: $build_cmd"
    if eval "$build_cmd"; then
        success "STAR index built successfully"
    else
        error "Failed to build STAR index"
        exit 1
    fi
}

# Run HISAT2 alignment
run_hisat2_rna_alignment() {
    log "Running HISAT2 RNA-seq alignment for sample: $SAMPLE_ID"
    
    # Build index if needed
    build_hisat2_index
    
    local log_file="${TEMP_DIR}/${SAMPLE_ID}_hisat2.log"
    local sam_file="${TEMP_DIR}/${SAMPLE_ID}_hisat2.sam"
    
    # Set HISAT2 parameters for RNA-seq
    local hisat2_opts=""
    
    # Strandness settings
    if [[ -n "$STRANDNESS" ]]; then
        case "$STRANDNESS" in
            "fr-firststrand"|"RF"|"reverse")
                hisat2_opts="$hisat2_opts --rna-strandness RF"
                ;;
            "fr-secondstrand"|"FR"|"forward")
                hisat2_opts="$hisat2_opts --rna-strandness FR"
                ;;
            "unstranded"|"none")
                # No strandness option
                ;;
        esac
    fi
    
    # Troubleshooting mode with more permissive settings
    if [[ "$TROUBLESHOOTING_MODE" == "true" ]]; then
        warning "Running in troubleshooting mode with permissive settings"
        hisat2_opts="$hisat2_opts --mp 6,2 --rdg 5,3 --rfg 5,3 --score-min L,0,-0.6"
    fi
    
    # Additional RNA-seq optimizations
    hisat2_opts="$hisat2_opts --dta --no-discordant --no-mixed"
    
    # Build alignment command
    local alignment_cmd="hisat2 -x $HISAT2_INDEX -p $THREADS $hisat2_opts"
    
    if [[ "$SINGLE_END" == "true" ]]; then
        alignment_cmd="$alignment_cmd -U $READ1"
    else
        alignment_cmd="$alignment_cmd -1 $READ1 -2 $READ2"
    fi
    
    alignment_cmd="$alignment_cmd -S $sam_file"
    
    log "Executing: $alignment_cmd"
    
    # Run alignment
    if eval "$alignment_cmd" 2>"$log_file"; then
        success "HISAT2 alignment completed"
    else
        error "HISAT2 alignment failed"
        if [[ -f "$log_file" ]]; then
            error "Check log file: $log_file"
            tail -10 "$log_file"
        fi
        exit 1
    fi
    
    # Convert SAM to BAM and sort
    log "Converting SAM to sorted BAM"
    if samtools view -bS "$sam_file" | samtools sort -@ "$THREADS" -m 1G -o "$OUTPUT_BAM" -; then
        success "BAM conversion and sorting completed"
        rm -f "$sam_file"  # Clean up SAM file
    else
        error "BAM conversion failed"
        exit 1
    fi
    
    # Extract alignment statistics
    extract_hisat2_stats "$log_file"
}

# Run STAR alignment
run_star_rna_alignment() {
    log "Running STAR RNA-seq alignment for sample: $SAMPLE_ID"
    
    # Build index if needed
    build_star_index
    
    local star_prefix="${TEMP_DIR}/${SAMPLE_ID}_star_"
    
    # Set STAR parameters for RNA-seq
    local star_opts=""
    star_opts="$star_opts --outSAMtype BAM SortedByCoordinate"
    star_opts="$star_opts --outFilterMismatchNmax 2"
    star_opts="$star_opts --alignSJoverhangMin 8"
    star_opts="$star_opts --alignSJDBoverhangMin 1"
    star_opts="$star_opts --outFilterMultimapNmax 20"
    star_opts="$star_opts --alignIntronMin 20"
    star_opts="$star_opts --alignIntronMax 1000000"
    
    # Strandness settings
    if [[ -n "$STRANDNESS" ]]; then
        case "$STRANDNESS" in
            "fr-firststrand"|"RF"|"reverse")
                star_opts="$star_opts --outSAMstrandField intronMotif"
                ;;
            "fr-secondstrand"|"FR"|"forward")
                star_opts="$star_opts --outSAMstrandField intronMotif"
                ;;
        esac
    fi
    
    # Troubleshooting mode
    if [[ "$TROUBLESHOOTING_MODE" == "true" ]]; then
        warning "Running in troubleshooting mode with permissive settings"
        star_opts="$star_opts --outFilterMismatchNmax 5"
        star_opts="$star_opts --outFilterMultimapNmax 50"
    fi
    
    # Memory optimization for limited systems
    if [[ $MEMORY_GB -lt 32 ]]; then
        star_opts="$star_opts --limitBAMsortRAM $((MEMORY_GB * 1000000000))"
    fi
    
    # Build alignment command
    local alignment_cmd="STAR --genomeDir $STAR_INDEX"
    alignment_cmd="$alignment_cmd --runThreadN $THREADS"
    alignment_cmd="$alignment_cmd --outFileNamePrefix $star_prefix"
    alignment_cmd="$alignment_cmd $star_opts"
    
    if [[ "$SINGLE_END" == "true" ]]; then
        alignment_cmd="$alignment_cmd --readFilesIn $READ1"
    else
        alignment_cmd="$alignment_cmd --readFilesIn $READ1 $READ2"
    fi
    
    # Handle compressed files
    if [[ "$READ1" == *.gz ]]; then
        alignment_cmd="$alignment_cmd --readFilesCommand zcat"
    fi
    
    log "Executing: $alignment_cmd"
    
    # Run alignment
    if eval "$alignment_cmd"; then
        success "STAR alignment completed"
    else
        error "STAR alignment failed"
        if [[ -f "${star_prefix}Log.final.out" ]]; then
            error "Check log file: ${star_prefix}Log.final.out"
            tail -10 "${star_prefix}Log.final.out"
        fi
        exit 1
    fi
    
    # Move output BAM to final location
    local star_bam="${star_prefix}Aligned.sortedByCoord.out.bam"
    if [[ -f "$star_bam" ]]; then
        mv "$star_bam" "$OUTPUT_BAM"
        success "STAR output moved to: $OUTPUT_BAM"
    else
        error "STAR output BAM not found: $star_bam"
        exit 1
    fi
    
    # Extract alignment statistics
    extract_star_stats "${star_prefix}Log.final.out"
}

# Run minimap2 alignment
run_minimap2_rna_alignment() {
    log "Running minimap2 RNA-seq alignment for sample: $SAMPLE_ID"
    
    local log_file="${TEMP_DIR}/${SAMPLE_ID}_minimap2.log"
    local sam_file="${TEMP_DIR}/${SAMPLE_ID}_minimap2.sam"
    
    # Set minimap2 parameters for RNA-seq
    local minimap2_opts="-ax splice -uf -k14"
    
    # Strandness settings
    if [[ -n "$STRANDNESS" ]]; then
        case "$STRANDNESS" in
            "fr-firststrand"|"RF"|"reverse")
                minimap2_opts="$minimap2_opts --for-only"
                ;;
            "fr-secondstrand"|"FR"|"forward")
                minimap2_opts="$minimap2_opts --rev-only"
                ;;
        esac
    fi
    
    # Build alignment command
    local alignment_cmd="minimap2 $minimap2_opts -t $THREADS $GENOME_FASTA"
    
    if [[ "$SINGLE_END" == "true" ]]; then
        alignment_cmd="$alignment_cmd $READ1"
    else
        alignment_cmd="$alignment_cmd $READ1 $READ2"
    fi
    
    alignment_cmd="$alignment_cmd > $sam_file"
    
    log "Executing: $alignment_cmd"
    
    # Run alignment
    if eval "$alignment_cmd" 2>"$log_file"; then
        success "minimap2 alignment completed"
    else
        error "minimap2 alignment failed"
        if [[ -f "$log_file" ]]; then
            error "Check log file: $log_file"
            tail -10 "$log_file"
        fi
        exit 1
    fi
    
    # Convert SAM to BAM and sort
    log "Converting SAM to sorted BAM"
    if samtools view -bS "$sam_file" | samtools sort -@ "$THREADS" -m 512M -o "$OUTPUT_BAM" -; then
        success "BAM conversion and sorting completed"
        rm -f "$sam_file"  # Clean up SAM file
    else
        error "BAM conversion failed"
        exit 1
    fi
    
    # Extract alignment statistics
    extract_minimap2_stats "$log_file"
}

# Extract HISAT2 statistics
extract_hisat2_stats() {
    local log_file=$1
    
    if [[ ! -f "$log_file" ]]; then
        warning "HISAT2 log file not found: $log_file"
        return
    fi
    
    local stats_file="${TEMP_DIR}/$(basename "$OUTPUT_BAM" .bam)_hisat2_stats.txt"
    
    # Parse HISAT2 output
    local total_reads overall_alignment_rate
    total_reads=$(grep "reads; of these:" "$log_file" | awk '{print $1}' || echo "Unknown")
    overall_alignment_rate=$(grep "overall alignment rate" "$log_file" | awk '{print $1}' || echo "Unknown")
    
    cat > "$stats_file" << EOF
Sample: $SAMPLE_ID
Aligner: HISAT2
Total_reads: $total_reads
Overall_alignment_rate: $overall_alignment_rate
Log_file: $log_file
BAM_file: $OUTPUT_BAM
EOF
    
    log "HISAT2 alignment rate: $overall_alignment_rate"
}

# Extract STAR statistics
extract_star_stats() {
    local log_file=$1
    
    if [[ ! -f "$log_file" ]]; then
        warning "STAR log file not found: $log_file"
        return
    fi
    
    local stats_file="${OUTPUT_DIR}/$(basename "$OUTPUT_BAM" .bam)_star_stats.txt"
    
    # Parse STAR output
    local total_reads uniquely_mapped overall_mapped
    total_reads=$(grep "Number of input reads" "$log_file" | awk '{print $NF}' || echo "Unknown")
    uniquely_mapped=$(grep "Uniquely mapped reads %" "$log_file" | awk '{print $NF}' || echo "Unknown")
    overall_mapped=$(grep "% of reads mapped to multiple loci" "$log_file" | awk '{print $NF}' || echo "Unknown")
    
    cat > "$stats_file" << EOF
Sample: $SAMPLE_ID
Aligner: STAR
Total_reads: $total_reads
Uniquely_mapped_percent: $uniquely_mapped
Multiple_mapped_percent: $overall_mapped
Log_file: $log_file
BAM_file: $OUTPUT_BAM
EOF
    
    log "STAR uniquely mapped: $uniquely_mapped"
}

# Extract minimap2 statistics
extract_minimap2_stats() {
    local log_file=$1
    
    if [[ ! -f "$log_file" ]]; then
        warning "minimap2 log file not found: $log_file"
        return
    fi
    
    local stats_file="${OUTPUT_DIR}/$(basename "$OUTPUT_BAM" .bam)_minimap2_stats.txt"
    
    # Parse minimap2 output (basic stats from stderr)
    local mapped_reads total_reads
    mapped_reads=$(grep "mapped" "$log_file" | tail -1 | awk '{print $1}' || echo "Unknown")
    total_reads=$(grep "processed" "$log_file" | tail -1 | awk '{print $1}' || echo "Unknown")
    
    cat > "$stats_file" << EOF
Sample: $SAMPLE_ID
Aligner: minimap2
Total_reads: $total_reads
Mapped_reads: $mapped_reads
Log_file: $log_file
BAM_file: $OUTPUT_BAM
EOF
    
    log "minimap2 mapped reads: $mapped_reads"
}

# Set RNA alignment parameters
set_rna_alignment_params() {
    # Auto-detect system resources
    if [[ $THREADS -eq 0 ]]; then
        THREADS=$(get_optimal_threads "alignment")
    fi
    
    if [[ $MEMORY_GB -eq 0 ]]; then
        MEMORY_GB=$(get_available_memory)
    fi
    
    # Create temporary directory
    if [[ -z "$TEMP_DIR" ]]; then
        TEMP_DIR="/tmp/rna_alignment_$$"
        mkdir -p "$TEMP_DIR"
    fi
    
    # Auto-detect aligner
    auto_detect_rna_aligner
    
    log "RNA-seq alignment parameters:"
    log "  Aligner: $ALIGNER"
    log "  Threads: $THREADS"
    log "  Memory: ${MEMORY_GB}GB"
    log "  Sample: $SAMPLE_ID"
    log "  Troubleshooting mode: $TROUBLESHOOTING_MODE"
}

# Post-alignment processing with RNA-seq quality checks
post_rna_alignment_processing() {
    log "Post-alignment processing for RNA-seq sample: $SAMPLE_ID"
    
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
    
    # RNA-seq specific quality assessment
    assess_rna_alignment_quality "$flagstat_file"
}

# Assess RNA-seq alignment quality
assess_rna_alignment_quality() {
    local flagstat_file=$1
    
    local total_reads mapped_reads properly_paired
    total_reads=$(grep "in total" "$flagstat_file" | awk '{print $1}')
    mapped_reads=$(grep "mapped (" "$flagstat_file" | head -1 | awk '{print $1}')
    properly_paired=$(grep "properly paired" "$flagstat_file" | awk '{print $1}')
    
    if [[ -n "$total_reads" && "$total_reads" -gt 0 ]]; then
        local mapping_rate paired_rate
        mapping_rate=$(echo "scale=0; $mapped_reads * 100 / $total_reads" | bc -l 2>/dev/null || echo "0")
        paired_rate=$(echo "scale=0; $properly_paired * 100 / $total_reads" | bc -l 2>/dev/null || echo "0")
        
        log "RNA-seq Quality Assessment:"
        log "  Total reads: $total_reads"
        log "  Mapped reads: $mapped_reads (${mapping_rate}%)"
        log "  Properly paired: $properly_paired (${paired_rate}%)"
        
        # RNA-seq specific thresholds
        if [[ "${mapping_rate%.*}" -lt 70 ]]; then
            warning "Low mapping rate for RNA-seq: ${mapping_rate}%"
            warning "Consider checking:"
            warning "  - Reference genome version"
            warning "  - rRNA contamination"
            warning "  - Library preparation protocol"
        elif [[ "${mapping_rate%.*}" -lt 85 ]]; then
            log "Moderate mapping rate for RNA-seq: ${mapping_rate}%"
        else
            success "Good mapping rate for RNA-seq: ${mapping_rate}%"
        fi
        
        # Suggest running diagnostics if mapping rate is still low
        if [[ "${mapping_rate%.*}" -lt 70 ]]; then
            log "Consider running alignment diagnostics:"
            log "  bash scripts/alignment_diagnostics.sh --bam $OUTPUT_BAM --sample $SAMPLE_ID --output diagnostics/"
        fi
    fi
}

# Main RNA alignment function
run_rna_alignment() {
    # Parse arguments
    parse_rna_alignment_args "$@"
    
    # Set parameters
    set_rna_alignment_params
    
    # Validate parameters
    validate_rna_alignment_params
    
    log "Starting RNA-seq alignment process"
    log "Input: $READ1 $([ "$SINGLE_END" == "false" ] && echo "$READ2")"
    log "Output: $OUTPUT_BAM"
    
    # Run alignment based on selected aligner
    case "$ALIGNER" in
        "hisat2")
            run_hisat2_rna_alignment
            ;;
        "star")
            run_star_rna_alignment
            ;;
        "minimap2")
            run_minimap2_rna_alignment
            ;;
        *)
            error "Unknown RNA-seq aligner: $ALIGNER"
            exit 1
            ;;
    esac
    
    # Post-processing
    post_rna_alignment_processing
    
    # Cleanup temporary files
    if [[ -n "$TEMP_DIR" && -d "$TEMP_DIR" ]]; then
        rm -rf "$TEMP_DIR"
        log "Temporary files cleaned"
    fi
    
    success "RNA-seq alignment module completed successfully"
    return 0
}

# Usage function
rna_alignment_usage() {
    cat << EOF
RNA-seq Alignment Module Usage:

$0 --read1 FILE --output BAM --sample-id ID [options]

Required arguments:
  --read1 FILE          Input FASTQ file (read 1)
  --output BAM          Output BAM file path
  --sample-id ID        Sample identifier

Optional arguments:
  --read2 FILE          Input FASTQ file (read 2, for paired-end)
  --single-end          Single-end sequencing mode
  --aligner TOOL        RNA-seq aligner (hisat2|star|minimap2|auto) [default: auto]
  --genome-fasta FILE   Genome FASTA file
  --annotation-gtf FILE Gene annotation GTF file
  --hisat2-index PREFIX HISAT2 index prefix (if pre-built)
  --star-index DIR      STAR index directory (if pre-built)
  --strandness TYPE     Library strandness (fr-firststrand|fr-secondstrand|unstranded)
  --threads N           Number of threads [default: auto-detect]
  --memory N            Available memory in GB [default: auto-detect]
  --troubleshooting     Use permissive alignment parameters
  --verbose             Verbose output

Examples:
  # Auto-detect aligner with genome FASTA
  $0 --read1 sample_1.fq.gz --read2 sample_2.fq.gz \\
     --output sample.bam --sample-id sample \\
     --genome-fasta genome.fa --annotation-gtf genes.gtf

  # Use pre-built HISAT2 index
  $0 --read1 sample_1.fq.gz --read2 sample_2.fq.gz \\
     --output sample.bam --sample-id sample \\
     --aligner hisat2 --hisat2-index /path/to/hisat2_index

  # Troubleshooting mode for difficult samples
  $0 --read1 sample_1.fq.gz --read2 sample_2.fq.gz \\
     --output sample.bam --sample-id sample \\
     --troubleshooting --genome-fasta genome.fa

EOF
}

# Run RNA alignment if script is called directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    if [[ $# -eq 0 ]]; then
        rna_alignment_usage
        exit 1
    fi
    
    run_rna_alignment "$@"
fi