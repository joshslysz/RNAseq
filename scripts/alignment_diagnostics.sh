#!/bin/bash

# RNA-seq Alignment Diagnostics Tool
# Author: Joshua Slysz, PhD - Dalhousie University
# Comprehensive analysis of low mapping rates and alignment issues

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m'

# Logging functions
log() { echo -e "${BLUE}[DIAG]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
info() { echo -e "${CYAN}[INFO]${NC} $1"; }

# Default parameters
BAM_FILE=""
SAMPLE_ID=""
OUTPUT_DIR=""
FASTQ1=""
FASTQ2=""
REFERENCE_FASTA=""
ANNOTATION_GTF=""
NUM_UNMAPPED_READS=1000
VERBOSE=false

usage() {
    cat << EOF
RNA-seq Alignment Diagnostics Tool

Usage: $0 --bam BAM_FILE --sample SAMPLE_ID --output OUTPUT_DIR [options]

Required arguments:
  --bam FILE          Input BAM file to analyze
  --sample ID         Sample identifier
  --output DIR        Output directory for diagnostic reports

Optional arguments:
  --fastq1 FILE       Original FASTQ file (read 1) for unmapped read analysis
  --fastq2 FILE       Original FASTQ file (read 2) for paired-end data
  --reference FILE    Reference genome FASTA file
  --annotation FILE   Gene annotation GTF file
  --num-unmapped N    Number of unmapped reads to extract [default: 1000]
  --verbose           Verbose output

Examples:
  # Basic diagnostics
  $0 --bam sample.bam --sample SRR123 --output diagnostics/

  # Comprehensive analysis with unmapped read extraction
  $0 --bam sample.bam --sample SRR123 --output diagnostics/ \\
     --fastq1 sample_1.fq.gz --fastq2 sample_2.fq.gz \\
     --reference genome.fa --annotation genes.gtf

EOF
}

# Parse command line arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --bam)
                BAM_FILE="$2"
                shift 2
                ;;
            --sample)
                SAMPLE_ID="$2"
                shift 2
                ;;
            --output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --fastq1)
                FASTQ1="$2"
                shift 2
                ;;
            --fastq2)
                FASTQ2="$2"
                shift 2
                ;;
            --reference)
                REFERENCE_FASTA="$2"
                shift 2
                ;;
            --annotation)
                ANNOTATION_GTF="$2"
                shift 2
                ;;
            --num-unmapped)
                NUM_UNMAPPED_READS="$2"
                shift 2
                ;;
            --verbose)
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
}

# Validate parameters
validate_params() {
    local errors=0
    
    if [[ -z "$BAM_FILE" ]]; then
        error "BAM file must be specified with --bam"
        errors=$((errors + 1))
    elif [[ ! -f "$BAM_FILE" ]]; then
        error "BAM file not found: $BAM_FILE"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$SAMPLE_ID" ]]; then
        error "Sample ID must be specified with --sample"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        error "Output directory must be specified with --output"
        errors=$((errors + 1))
    fi
    
    if [[ $errors -gt 0 ]]; then
        exit 1
    fi
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
}

# Generate comprehensive alignment statistics
analyze_alignment_stats() {
    log "Analyzing alignment statistics for $SAMPLE_ID"
    
    local stats_file="${OUTPUT_DIR}/${SAMPLE_ID}_detailed_stats.txt"
    local flagstat_file="${OUTPUT_DIR}/${SAMPLE_ID}_flagstat.txt"
    
    # Generate flagstat
    samtools flagstat "$BAM_FILE" > "$flagstat_file"
    
    # Parse statistics
    local total_reads mapped_reads properly_paired singletons
    total_reads=$(grep "in total" "$flagstat_file" | awk '{print $1}')
    mapped_reads=$(grep "mapped (" "$flagstat_file" | head -1 | awk '{print $1}')
    properly_paired=$(grep "properly paired" "$flagstat_file" | awk '{print $1}')
    singletons=$(grep "singletons" "$flagstat_file" | awk '{print $1}')
    
    # Calculate percentages
    local mapping_rate properly_paired_rate singleton_rate
    mapping_rate=$(echo "scale=2; $mapped_reads * 100 / $total_reads" | bc -l 2>/dev/null || echo "N/A")
    if [[ $mapped_reads -gt 0 ]]; then
        properly_paired_rate=$(echo "scale=2; $properly_paired * 100 / $mapped_reads" | bc -l 2>/dev/null || echo "N/A")
        singleton_rate=$(echo "scale=2; $singletons * 100 / $mapped_reads" | bc -l 2>/dev/null || echo "N/A")
    else
        properly_paired_rate="N/A"
        singleton_rate="N/A"
    fi
    
    # Write detailed statistics
    cat > "$stats_file" << EOF
ALIGNMENT STATISTICS SUMMARY
============================
Sample: $SAMPLE_ID
BAM file: $BAM_FILE
Analysis date: $(date)

BASIC STATISTICS:
- Total reads: $total_reads
- Mapped reads: $mapped_reads
- Mapping rate: ${mapping_rate}%
- Properly paired: $properly_paired
- Properly paired rate: ${properly_paired_rate}%
- Singletons: $singletons
- Singleton rate: ${singleton_rate}%

QUALITY ASSESSMENT:
EOF

    # Quality assessment
    if [[ "${mapping_rate%.*}" -lt 50 ]] 2>/dev/null; then
        echo "- Mapping quality: POOR (${mapping_rate}%)" >> "$stats_file"
        echo "- Recommendation: Check reference genome, read quality, and alignment parameters" >> "$stats_file"
    elif [[ "${mapping_rate%.*}" -lt 70 ]] 2>/dev/null; then
        echo "- Mapping quality: MODERATE (${mapping_rate}%)" >> "$stats_file"
        echo "- Recommendation: May be acceptable depending on experimental conditions" >> "$stats_file"
    else
        echo "- Mapping quality: GOOD (${mapping_rate}%)" >> "$stats_file"
    fi
    
    success "Alignment statistics saved to: $stats_file"
    
    # Display key findings
    info "KEY FINDINGS:"
    info "  Total reads: $total_reads"
    info "  Mapping rate: ${mapping_rate}%"
    info "  Properly paired: ${properly_paired_rate}%"
}

# Extract unmapped reads for analysis
extract_unmapped_reads() {
    if [[ -z "$FASTQ1" ]]; then
        warning "FASTQ files not provided, skipping unmapped read extraction"
        return
    fi
    
    log "Extracting unmapped reads for analysis"
    
    local unmapped_dir="${OUTPUT_DIR}/unmapped_analysis"
    mkdir -p "$unmapped_dir"
    
    # Extract unmapped read IDs
    local unmapped_ids="${unmapped_dir}/${SAMPLE_ID}_unmapped_ids.txt"
    samtools view -f 4 "$BAM_FILE" | head -n "$NUM_UNMAPPED_READS" | cut -f1 | sort -u > "$unmapped_ids"
    
    local num_ids
    num_ids=$(wc -l < "$unmapped_ids")
    info "Extracted $num_ids unique unmapped read IDs"
    
    # Extract corresponding FASTQ sequences
    local unmapped_fastq1="${unmapped_dir}/${SAMPLE_ID}_unmapped_1.fastq"
    local unmapped_fastq2="${unmapped_dir}/${SAMPLE_ID}_unmapped_2.fastq"
    
    log "Extracting unmapped sequences from FASTQ files..."
    
    # Use seqtk to extract reads by ID (if available), otherwise use grep
    if command -v seqtk >/dev/null 2>&1; then
        if [[ "$FASTQ1" == *.gz ]]; then
            zcat "$FASTQ1" | seqtk subseq - "$unmapped_ids" > "$unmapped_fastq1"
        else
            seqtk subseq "$FASTQ1" "$unmapped_ids" > "$unmapped_fastq1"
        fi
        
        if [[ -n "$FASTQ2" ]]; then
            if [[ "$FASTQ2" == *.gz ]]; then
                zcat "$FASTQ2" | seqtk subseq - "$unmapped_ids" > "$unmapped_fastq2"
            else
                seqtk subseq "$FASTQ2" "$unmapped_ids" > "$unmapped_fastq2"
            fi
        fi
    else
        # Fallback to grep-based extraction (slower but works without seqtk)
        warning "seqtk not found, using slower grep-based extraction"
        extract_reads_with_grep "$FASTQ1" "$unmapped_ids" "$unmapped_fastq1"
        if [[ -n "$FASTQ2" ]]; then
            extract_reads_with_grep "$FASTQ2" "$unmapped_ids" "$unmapped_fastq2"
        fi
    fi
    
    # Analyze unmapped sequences
    analyze_unmapped_sequences "$unmapped_fastq1" "$unmapped_fastq2"
}

# Grep-based read extraction (fallback method)
extract_reads_with_grep() {
    local fastq_file=$1
    local id_file=$2
    local output_file=$3
    
    local temp_pattern="/tmp/pattern_$$"
    sed 's/^/@/' "$id_file" > "$temp_pattern"
    
    if [[ "$fastq_file" == *.gz ]]; then
        zcat "$fastq_file" | grep -A 3 -f "$temp_pattern" | grep -v "^--$" > "$output_file"
    else
        grep -A 3 -f "$temp_pattern" "$fastq_file" | grep -v "^--$" > "$output_file"
    fi
    
    rm -f "$temp_pattern"
}

# Analyze unmapped sequences for common issues
analyze_unmapped_sequences() {
    local unmapped_1=$1
    local unmapped_2=$2
    
    log "Analyzing unmapped sequences for common issues"
    
    local analysis_file="${OUTPUT_DIR}/unmapped_analysis/${SAMPLE_ID}_unmapped_analysis.txt"
    
    cat > "$analysis_file" << EOF
UNMAPPED READ ANALYSIS
=====================
Sample: $SAMPLE_ID
Number of reads analyzed: $NUM_UNMAPPED_READS
Analysis date: $(date)

EOF
    
    # Basic sequence statistics
    local num_reads avg_length gc_content
    num_reads=$(grep -c "^@" "$unmapped_1" 2>/dev/null || echo "0")
    
    if [[ $num_reads -gt 0 ]]; then
        # Calculate average length
        avg_length=$(awk 'NR%4==2 {sum+=length($0); count++} END {if(count>0) print sum/count; else print 0}' "$unmapped_1")
        
        # Calculate GC content
        gc_content=$(awk 'NR%4==2 {
            gsub(/[^GCgc]/, "", $0); 
            gc+=length($0); 
            gsub(/[^ATGCatgc]/, "", $0); 
            total+=length($0)
        } END {
            if(total>0) print (gc/total)*100; else print 0
        }' "$unmapped_1")
        
        cat >> "$analysis_file" << EOF
SEQUENCE CHARACTERISTICS:
- Number of unmapped reads: $num_reads
- Average read length: ${avg_length%.??}
- GC content: ${gc_content%.??}%

EOF
    fi
    
    # Check for adapter contamination
    check_adapter_contamination "$unmapped_1" "$analysis_file"
    
    # Check for low complexity sequences
    check_low_complexity "$unmapped_1" "$analysis_file"
    
    # Check for ribosomal RNA contamination
    check_rrna_contamination "$unmapped_1" "$analysis_file"
    
    success "Unmapped sequence analysis saved to: $analysis_file"
}

# Check for adapter contamination
check_adapter_contamination() {
    local fastq_file=$1
    local output_file=$2
    
    # Common Illumina adapters
    local adapters=(
        "AGATCGGAAGAGC"  # TruSeq Universal Adapter
        "CTGTCTCTTATA"   # Nextera Transposase
        "GATCGGAAGAGC"   # Alternative TruSeq
    )
    
    echo "ADAPTER CONTAMINATION CHECK:" >> "$output_file"
    
    for adapter in "${adapters[@]}"; do
        local count
        count=$(awk -v adapter="$adapter" 'NR%4==2 && index($0, adapter) > 0 {count++} END {print count+0}' "$fastq_file")
        local percentage
        if [[ $num_reads -gt 0 ]]; then
            percentage=$(echo "scale=2; $count * 100 / $num_reads" | bc -l 2>/dev/null || echo "0")
        else
            percentage="0"
        fi
        echo "- $adapter: $count reads (${percentage}%)" >> "$output_file"
    done
    echo "" >> "$output_file"
}

# Check for low complexity sequences
check_low_complexity() {
    local fastq_file=$1
    local output_file=$2
    
    echo "LOW COMPLEXITY SEQUENCES:" >> "$output_file"
    
    # Count poly-A, poly-T, poly-G, poly-C sequences
    local poly_a poly_t poly_g poly_c
    poly_a=$(awk 'NR%4==2 && gsub(/A/, "&") >= length($0)*0.8 {count++} END {print count+0}' "$fastq_file")
    poly_t=$(awk 'NR%4==2 && gsub(/T/, "&") >= length($0)*0.8 {count++} END {print count+0}' "$fastq_file")
    poly_g=$(awk 'NR%4==2 && gsub(/G/, "&") >= length($0)*0.8 {count++} END {print count+0}' "$fastq_file")
    poly_c=$(awk 'NR%4==2 && gsub(/C/, "&") >= length($0)*0.8 {count++} END {print count+0}' "$fastq_file")
    
    echo "- Poly-A sequences: $poly_a" >> "$output_file"
    echo "- Poly-T sequences: $poly_t" >> "$output_file"
    echo "- Poly-G sequences: $poly_g" >> "$output_file"
    echo "- Poly-C sequences: $poly_c" >> "$output_file"
    echo "" >> "$output_file"
}

# Check for ribosomal RNA contamination
check_rrna_contamination() {
    local fastq_file=$1
    local output_file=$2
    
    echo "rRNA CONTAMINATION CHECK:" >> "$output_file"
    
    # Common human rRNA sequences (partial)
    local rrna_patterns=(
        "CCTACCTGGTTGATCCTGCCAGTA"  # 18S rRNA
        "GCTCGAATTTCTACTAAGCGAA"    # 28S rRNA
        "GTCCCTGCCCTTTGTACACA"      # 5.8S rRNA
    )
    
    for pattern in "${rrna_patterns[@]}"; do
        local count
        count=$(awk -v pattern="$pattern" 'NR%4==2 && index($0, pattern) > 0 {count++} END {print count+0}' "$fastq_file")
        local percentage
        if [[ $num_reads -gt 0 ]]; then
            percentage=$(echo "scale=2; $count * 100 / $num_reads" | bc -l 2>/dev/null || echo "0")
        else
            percentage="0"
        fi
        echo "- rRNA pattern match: $count reads (${percentage}%)" >> "$output_file"
    done
    echo "" >> "$output_file"
}

# Diagnose alignment parameters
diagnose_alignment_parameters() {
    log "Diagnosing alignment parameters and suggesting improvements"
    
    local diagnosis_file="${OUTPUT_DIR}/${SAMPLE_ID}_alignment_diagnosis.txt"
    
    cat > "$diagnosis_file" << EOF
ALIGNMENT PARAMETER DIAGNOSIS
============================
Sample: $SAMPLE_ID
Analysis date: $(date)

CURRENT ISSUES IDENTIFIED:
1. Using DNA alignment parameters for RNA-seq data
   - Current: Bowtie2 with --sensitive --no-discordant --no-mixed
   - Problem: Cannot handle splice junctions in RNA-seq data
   
2. Missing splice-aware alignment
   - RNA-seq requires splice-aware aligners
   - Introns cause reads to span exon-exon junctions
   
RECOMMENDED SOLUTIONS:

OPTION 1: Switch to HISAT2 (Recommended)
- Splice-aware aligner designed for RNA-seq
- Memory efficient
- Commands:
  hisat2-build genome.fa genome_index
  hisat2 -x genome_index -1 reads_1.fq.gz -2 reads_2.fq.gz -S output.sam

OPTION 2: Use STAR aligner
- Very fast and accurate for RNA-seq
- Requires more memory (30GB+ for human genome)
- Commands:
  STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf
  STAR --genomeDir star_index --readFilesIn reads_1.fq.gz reads_2.fq.gz --readFilesCommand zcat

OPTION 3: TopHat2 (Legacy, not recommended)
- Older splice-aware aligner
- Built on top of Bowtie2
- Slower than HISAT2/STAR

IMMEDIATE ACTIONS:
1. Install HISAT2 or STAR
2. Build splice-aware genome index
3. Re-run alignment with RNA-seq specific parameters
4. Expect mapping rates of 80-95% for good quality RNA-seq data

PARAMETER RECOMMENDATIONS:
- For HISAT2: Use default parameters (optimized for RNA-seq)
- For STAR: --outFilterMismatchNmax 2 --alignSJoverhangMin 8
- Consider using gene annotation GTF for better splice site detection

EOF
    
    success "Alignment diagnosis saved to: $diagnosis_file"
    
    # Display immediate recommendations
    warning "CRITICAL ISSUE DETECTED:"
    warning "  Current pipeline uses DNA alignment for RNA-seq data"
    warning "  This explains the 0.82% mapping rate"
    info "IMMEDIATE FIX:"
    info "  1. Install HISAT2: conda install hisat2"
    info "  2. Build HISAT2 index: hisat2-build genome.fa genome_index"
    info "  3. Re-align with HISAT2 instead of Bowtie2"
}

# Generate corrected alignment commands
generate_corrected_commands() {
    log "Generating corrected alignment commands"
    
    local commands_file="${OUTPUT_DIR}/${SAMPLE_ID}_corrected_commands.sh"
    
    cat > "$commands_file" << EOF
#!/bin/bash

# Corrected RNA-seq Alignment Commands
# Generated by alignment diagnostics tool

set -euo pipefail

# Variables (adjust paths as needed)
SAMPLE_ID="$SAMPLE_ID"
FASTQ1="$FASTQ1"
FASTQ2="$FASTQ2"
REFERENCE_FASTA="$REFERENCE_FASTA"
ANNOTATION_GTF="$ANNOTATION_GTF"
THREADS=\${THREADS:-6}
OUTPUT_DIR="\${OUTPUT_DIR:-aligned}"

mkdir -p "\$OUTPUT_DIR"

echo "Starting corrected RNA-seq alignment for \$SAMPLE_ID"

# Option 1: HISAT2 (Recommended)
echo "Building HISAT2 index..."
hisat2-build "\$REFERENCE_FASTA" "\${OUTPUT_DIR}/hisat2_index/genome"

echo "Running HISAT2 alignment..."
hisat2 -x "\${OUTPUT_DIR}/hisat2_index/genome" \\
       -1 "\$FASTQ1" \\
       -2 "\$FASTQ2" \\
       -p "\$THREADS" \\
       --rna-strandness RF \\
       -S "\${OUTPUT_DIR}/\${SAMPLE_ID}_hisat2.sam"

# Convert to BAM and sort
samtools view -bS "\${OUTPUT_DIR}/\${SAMPLE_ID}_hisat2.sam" | \\
samtools sort -@ "\$THREADS" -o "\${OUTPUT_DIR}/\${SAMPLE_ID}_hisat2.sorted.bam" -

# Index BAM file
samtools index "\${OUTPUT_DIR}/\${SAMPLE_ID}_hisat2.sorted.bam"

# Generate statistics
samtools flagstat "\${OUTPUT_DIR}/\${SAMPLE_ID}_hisat2.sorted.bam" > \\
    "\${OUTPUT_DIR}/\${SAMPLE_ID}_hisat2.flagstat"

echo "HISAT2 alignment completed!"
echo "Expected mapping rate: 80-95%"

# Clean up SAM file
rm -f "\${OUTPUT_DIR}/\${SAMPLE_ID}_hisat2.sam"

# Option 2: STAR (Alternative)
# Uncomment if you prefer STAR:
#
# echo "Building STAR index..."
# mkdir -p "\${OUTPUT_DIR}/star_index"
# STAR --runMode genomeGenerate \\
#      --genomeDir "\${OUTPUT_DIR}/star_index" \\
#      --genomeFastaFiles "\$REFERENCE_FASTA" \\
#      --sjdbGTFfile "\$ANNOTATION_GTF" \\
#      --runThreadN "\$THREADS"
#
# echo "Running STAR alignment..."
# STAR --genomeDir "\${OUTPUT_DIR}/star_index" \\
#      --readFilesIn "\$FASTQ1" "\$FASTQ2" \\
#      --readFilesCommand zcat \\
#      --runThreadN "\$THREADS" \\
#      --outFileNamePrefix "\${OUTPUT_DIR}/\${SAMPLE_ID}_star_" \\
#      --outSAMtype BAM SortedByCoordinate \\
#      --outFilterMismatchNmax 2 \\
#      --alignSJoverhangMin 8

EOF
    
    chmod +x "$commands_file"
    success "Corrected alignment commands saved to: $commands_file"
    
    info "To run the corrected alignment:"
    info "  bash $commands_file"
}

# Generate comprehensive diagnostic report
generate_diagnostic_report() {
    log "Generating comprehensive diagnostic report"
    
    local report_file="${OUTPUT_DIR}/${SAMPLE_ID}_diagnostic_report.html"
    
    cat > "$report_file" << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>RNA-seq Alignment Diagnostic Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .header { color: #2c3e50; border-bottom: 2px solid #3498db; }
        .section { margin: 20px 0; }
        .critical { background-color: #ffebee; border-left: 4px solid #f44336; padding: 10px; }
        .warning { background-color: #fff3e0; border-left: 4px solid #ff9800; padding: 10px; }
        .success { background-color: #e8f5e8; border-left: 4px solid #4caf50; padding: 10px; }
        .code { background-color: #f5f5f5; padding: 10px; border-radius: 4px; font-family: monospace; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <div class="header">
        <h1>RNA-seq Alignment Diagnostic Report</h1>
        <p>Sample: SAMPLE_ID | Date: DATE_PLACEHOLDER</p>
    </div>
    
    <div class="section critical">
        <h2>🚨 Critical Issues Detected</h2>
        <p><strong>Primary Issue:</strong> DNA alignment parameters used for RNA-seq data</p>
        <ul>
            <li>Current aligner: Bowtie2 with DNA-specific parameters</li>
            <li>RNA-seq requires splice-aware alignment</li>
            <li>Resulting mapping rate: ~0.82% (Expected: 80-95%)</li>
        </ul>
    </div>
    
    <div class="section warning">
        <h2>⚠️ Immediate Actions Required</h2>
        <ol>
            <li>Install HISAT2 or STAR for splice-aware alignment</li>
            <li>Build appropriate genome index with splice junction information</li>
            <li>Re-run alignment with RNA-seq specific parameters</li>
        </ol>
    </div>
    
    <div class="section">
        <h2>📊 Alignment Statistics</h2>
        <table>
            <tr><th>Metric</th><th>Value</th><th>Assessment</th></tr>
            <tr><td>Total Reads</td><td id="total_reads">-</td><td>-</td></tr>
            <tr><td>Mapped Reads</td><td id="mapped_reads">-</td><td id="mapping_assessment">-</td></tr>
            <tr><td>Mapping Rate</td><td id="mapping_rate">-</td><td id="rate_assessment">-</td></tr>
            <tr><td>Properly Paired</td><td id="properly_paired">-</td><td>-</td></tr>
        </table>
    </div>
    
    <div class="section success">
        <h2>✅ Recommended Solution</h2>
        <p><strong>Use HISAT2 for splice-aware alignment:</strong></p>
        <div class="code">
# Install HISAT2<br>
conda install hisat2<br><br>
# Build index<br>
hisat2-build genome.fa genome_index<br><br>
# Align reads<br>
hisat2 -x genome_index -1 reads_1.fq.gz -2 reads_2.fq.gz -S output.sam
        </div>
    </div>
    
    <div class="section">
        <h2>🔬 Technical Details</h2>
        <h3>Why Bowtie2 Failed for RNA-seq:</h3>
        <ul>
            <li>Bowtie2 is designed for DNA sequencing (genomic, ChIP-seq, etc.)</li>
            <li>RNA-seq reads span exon-exon junctions due to splicing</li>
            <li>Bowtie2 cannot align reads across introns</li>
            <li>This results in most RNA-seq reads appearing "unmapped"</li>
        </ul>
        
        <h3>How Splice-aware Aligners Work:</h3>
        <ul>
            <li>HISAT2/STAR can detect and align across splice junctions</li>
            <li>They use genome annotation to identify known splice sites</li>
            <li>They can also discover novel splice junctions</li>
            <li>This enables proper alignment of RNA-seq data</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>📋 Next Steps</h2>
        <ol>
            <li>Review the generated corrected alignment commands</li>
            <li>Install required splice-aware aligner (HISAT2 recommended)</li>
            <li>Re-run the pipeline with corrected parameters</li>
            <li>Verify improved mapping rates (expect 80-95%)</li>
            <li>Proceed with downstream analysis (quantification, etc.)</li>
        </ol>
    </div>
</body>
</html>
EOF
    
    # Replace placeholders with actual values
    sed -i "s/SAMPLE_ID/$SAMPLE_ID/g" "$report_file"
    sed -i "s/DATE_PLACEHOLDER/$(date)/g" "$report_file"
    
    success "Diagnostic report saved to: $report_file"
    info "Open in browser: file://$report_file"
}

# Main diagnostic function
main() {
    echo -e "${MAGENTA}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║                 RNA-seq Alignment Diagnostics               ║"
    echo "║              Comprehensive Analysis Tool                     ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
    
    # Parse arguments
    parse_args "$@"
    
    # Validate parameters
    validate_params
    
    # Run diagnostic analyses
    analyze_alignment_stats
    
    if [[ -n "$FASTQ1" ]]; then
        extract_unmapped_reads
    fi
    
    diagnose_alignment_parameters
    generate_corrected_commands
    generate_diagnostic_report
    
    echo ""
    success "🎉 Diagnostic analysis completed!"
    info "Results saved in: $OUTPUT_DIR"
    info ""
    info "Summary of findings:"
    info "  ❌ Main issue: Using DNA aligner for RNA-seq data"
    info "  ✅ Solution: Switch to HISAT2 or STAR"
    info "  📈 Expected improvement: 0.82% → 80-95% mapping rate"
    info ""
    info "Next steps:"
    info "  1. Review generated files in $OUTPUT_DIR"
    info "  2. Run corrected alignment commands"
    info "  3. Verify improved mapping rates"
}

# Run main function if script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    if [[ $# -eq 0 ]]; then
        usage
        exit 1
    fi
    
    main "$@"
fi