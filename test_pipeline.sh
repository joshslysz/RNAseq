#!/bin/bash

# Pipeline Testing and Validation Script
# Author: Joshua Slysz, PhD - Dalhousie University
# Comprehensive testing suite for memory-optimized RNA-seq pipeline
# Tests environment, modules, and end-to-end functionality

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

# Test parameters
TEST_TYPE="all"
VERBOSE=false
GENERATE_TEST_DATA=false
TEST_DATA_DIR="test_data"
CLEANUP=true
QUICK_TEST=false

# Test results
TESTS_PASSED=0
TESTS_FAILED=0
FAILED_TESTS=()

# Logging functions
log() { echo -e "${BLUE}[TEST]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
header() { echo -e "${BOLD}${BLUE}=== $1 ===${NC}"; }

# Test result tracking
test_passed() {
    local test_name=$1
    TESTS_PASSED=$((TESTS_PASSED + 1))
    success "✓ $test_name"
}

test_failed() {
    local test_name=$1
    TESTS_FAILED=$((TESTS_FAILED + 1))
    FAILED_TESTS+=("$test_name")
    error "✗ $test_name"
}

# Usage function
usage() {
    cat << EOF
${BOLD}Pipeline Testing and Validation Suite${NC}

Usage: $0 [options]

Options:
  -t, --test-type TYPE      Test type (env|modules|pipeline|all) [default: all]
  -g, --generate-data       Generate test data if not present
  -d, --test-data-dir DIR   Test data directory [default: test_data]
  -q, --quick               Quick test mode (minimal datasets)
  -v, --verbose             Verbose output
  -k, --keep-files          Keep test files after completion
  -h, --help                Show this help message

Test types:
  env          Environment and tool availability
  modules      Individual module functionality
  pipeline     End-to-end pipeline testing
  all          All tests (default)

Examples:
  # Run all tests
  $0

  # Test environment only
  $0 -t env

  # Generate test data and run pipeline tests
  $0 -t pipeline -g

  # Quick test with verbose output
  $0 -q -v
EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--test-type)
            TEST_TYPE="$2"
            shift 2
            ;;
        -g|--generate-data)
            GENERATE_TEST_DATA=true
            shift
            ;;
        -d|--test-data-dir)
            TEST_DATA_DIR="$2"
            shift 2
            ;;
        -q|--quick)
            QUICK_TEST=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -k|--keep-files)
            CLEANUP=false
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

# Test environment and dependencies
test_environment() {
    header "Testing Environment"
    
    # Check conda environment
    log "Checking conda environment..."
    if [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "$CONDA_DEFAULT_ENV" != "rnaseqpipe" ]]; then
        test_failed "Conda environment (rnaseqpipe not activated)"
        return 1
    else
        test_passed "Conda environment"
    fi
    
    # Check required tools
    local tools=("fastqc" "bowtie2" "minimap2" "samtools" "featureCounts" "multiqc" "python" "salmon")
    
    for tool in "${tools[@]}"; do
        log "Checking $tool..."
        if command -v "$tool" &> /dev/null; then
            local version
            case "$tool" in
                "fastqc")
                    version=$(fastqc --version 2>&1 | head -1 || echo "Unknown")
                    ;;
                "bowtie2")
                    version=$(bowtie2 --version 2>&1 | head -1 || echo "Unknown")
                    ;;
                "minimap2")
                    version=$(minimap2 --version 2>&1 || echo "Unknown")
                    ;;
                "samtools")
                    version=$(samtools --version 2>&1 | head -1 || echo "Unknown")
                    ;;
                "featureCounts")
                    version=$(featureCounts -v 2>&1 | head -1 || echo "Unknown")
                    ;;
                "multiqc")
                    version=$(multiqc --version 2>&1 || echo "Unknown")
                    ;;
                "python")
                    version=$(python --version 2>&1 || echo "Unknown")
                    ;;
                "salmon")
                    version=$(salmon --version 2>&1 | head -1 || echo "Unknown")
                    ;;
            esac
            
            if [[ "$VERBOSE" == "true" ]]; then
                log "  $tool: $version"
            fi
            test_passed "$tool availability"
        else
            test_failed "$tool availability"
        fi
    done
    
    # Check Python packages
    log "Checking Python packages..."
    local python_packages=("pandas" "numpy" "matplotlib" "seaborn" "pysam" "psutil")
    
    for package in "${python_packages[@]}"; do
        if python -c "import $package" 2>/dev/null; then
            test_passed "Python package: $package"
        else
            test_failed "Python package: $package"
        fi
    done
    
    # Check system resources
    log "Checking system resources..."
    local total_mem available_mem cpu_count
    total_mem=$(get_total_memory)
    available_mem=$(get_available_memory)
    cpu_count=$(get_cpu_count)
    
    if [[ $available_mem -ge 2 ]]; then
        test_passed "Memory requirement (${available_mem}GB available)"
    else
        test_failed "Memory requirement (${available_mem}GB available, minimum 2GB)"
    fi
    
    if [[ $cpu_count -ge 1 ]]; then
        test_passed "CPU requirement ($cpu_count cores)"
    else
        test_failed "CPU requirement ($cpu_count cores)"
    fi
}

# Generate test data
generate_test_data() {
    header "Generating Test Data"
    
    mkdir -p "$TEST_DATA_DIR"
    
    # Generate small test FASTQ files
    log "Generating test FASTQ files..."
    
    local read_count=1000
    if [[ "$QUICK_TEST" == "true" ]]; then
        read_count=100
    fi
    
    # Generate test reads using Python
    python3 << EOF
import random
import gzip

def generate_fastq(filename, num_reads, read_length=100):
    bases = ['A', 'T', 'G', 'C']
    with gzip.open(filename, 'wt') as f:
        for i in range(num_reads):
            header = f"@test_read_{i+1}"
            sequence = ''.join(random.choices(bases, k=read_length))
            plus = "+"
            quality = ''.join(['I'] * read_length)  # High quality scores
            f.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")

# Generate paired-end reads
generate_fastq("$TEST_DATA_DIR/test_sample_1.fq.gz", $read_count)
generate_fastq("$TEST_DATA_DIR/test_sample_2.fq.gz", $read_count)

# Generate single-end reads
generate_fastq("$TEST_DATA_DIR/test_single.fq.gz", $read_count)

print(f"Generated test FASTQ files with {$read_count} reads each")
EOF
    
    # Generate minimal reference files
    log "Generating minimal reference genome..."
    
    # Create tiny reference genome
    cat > "$TEST_DATA_DIR/test_genome.fa" << 'EOF'
>chr1
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
EOF
    
    # Create minimal GTF file
    cat > "$TEST_DATA_DIR/test_genes.gtf" << 'EOF'
chr1	test	gene	1	120	.	+	.	gene_id "TESTGENE001"; gene_name "TEST1";
chr1	test	transcript	1	120	.	+	.	gene_id "TESTGENE001"; transcript_id "TESTTX001";
chr1	test	exon	1	60	.	+	.	gene_id "TESTGENE001"; transcript_id "TESTTX001"; exon_number "1";
chr1	test	exon	61	120	.	+	.	gene_id "TESTGENE001"; transcript_id "TESTTX001"; exon_number "2";
chr1	test	gene	150	240	.	-	.	gene_id "TESTGENE002"; gene_name "TEST2";
chr1	test	transcript	150	240	.	-	.	gene_id "TESTGENE002"; transcript_id "TESTTX002";
chr1	test	exon	150	240	.	-	.	gene_id "TESTGENE002"; transcript_id "TESTTX002"; exon_number "1";
EOF
    
    success "Test data generated in: $TEST_DATA_DIR"
}

# Test individual modules
test_modules() {
    header "Testing Individual Modules"
    
    # Ensure test data exists
    if [[ ! -d "$TEST_DATA_DIR" ]] || [[ ! -f "$TEST_DATA_DIR/test_sample_1.fq.gz" ]]; then
        if [[ "$GENERATE_TEST_DATA" == "true" ]]; then
            generate_test_data
        else
            error "Test data not found. Use -g to generate test data."
            return 1
        fi
    fi
    
    # Test resource manager
    log "Testing resource manager..."
    if "$SCRIPT_DIR/scripts/resource_manager.sh" summary > /dev/null; then
        test_passed "Resource manager"
    else
        test_failed "Resource manager"
    fi
    
    # Test QC module (FastQC only for quick test)
    log "Testing QC module..."
    local qc_test_dir="$TEST_DATA_DIR/qc_test"
    mkdir -p "$qc_test_dir"
    
    if "$SCRIPT_DIR/scripts/qc_module.sh" \
        --qc-type fastq \
        --fastq-files "$TEST_DATA_DIR/test_sample_1.fq.gz,$TEST_DATA_DIR/test_sample_2.fq.gz" \
        --output-dir "$qc_test_dir" \
        --sample-id test_sample \
        --skip-rseqc \
        > /dev/null 2>&1; then
        test_passed "QC module (FastQC)"
    else
        test_failed "QC module (FastQC)"
    fi
    
    # Build test reference indices (minimap2 for speed)
    log "Building test reference indices..."
    if [[ ! -f "$TEST_DATA_DIR/test_genome.fa.fai" ]]; then
        samtools faidx "$TEST_DATA_DIR/test_genome.fa"
    fi
    
    # Test alignment module (minimap2 for quick test)
    log "Testing alignment module..."
    local align_test_dir="$TEST_DATA_DIR/align_test"
    mkdir -p "$align_test_dir"
    
    if "$SCRIPT_DIR/scripts/alignment_module.sh" \
        --aligner minimap2 \
        --genome-fasta "$TEST_DATA_DIR/test_genome.fa" \
        --read1 "$TEST_DATA_DIR/test_sample_1.fq.gz" \
        --read2 "$TEST_DATA_DIR/test_sample_2.fq.gz" \
        --output "$align_test_dir/test_sample.bam" \
        --sample-id test_sample \
        > /dev/null 2>&1; then
        test_passed "Alignment module (minimap2)"
    else
        test_failed "Alignment module (minimap2)"
    fi
    
    # Test quantification module
    if [[ -f "$align_test_dir/test_sample.bam" ]]; then
        log "Testing quantification module..."
        local quant_test_dir="$TEST_DATA_DIR/quant_test"
        mkdir -p "$quant_test_dir"
        
        if "$SCRIPT_DIR/scripts/quantification_module.sh" \
            --quantifier featurecounts \
            --bam-file "$align_test_dir/test_sample.bam" \
            --genome-gtf "$TEST_DATA_DIR/test_genes.gtf" \
            --output "$quant_test_dir/test_sample.counts.txt" \
            --sample-id test_sample \
            > /dev/null 2>&1; then
            test_passed "Quantification module (featureCounts)"
        else
            test_failed "Quantification module (featureCounts)"
        fi
    else
        test_failed "Quantification module (no BAM file from alignment)"
    fi
}

# Test end-to-end pipeline
test_pipeline() {
    header "Testing End-to-End Pipeline"
    
    # Ensure test data exists
    if [[ ! -d "$TEST_DATA_DIR" ]] || [[ ! -f "$TEST_DATA_DIR/test_sample_1.fq.gz" ]]; then
        if [[ "$GENERATE_TEST_DATA" == "true" ]]; then
            generate_test_data
        else
            error "Test data not found. Use -g to generate test data."
            return 1
        fi
    fi
    
    # Create test genome configuration
    log "Creating test genome configuration..."
    local test_config="$TEST_DATA_DIR/test_config.txt"
    cat > "$test_config" << EOF
GENOME_NAME=test
GENOME_FASTA=$TEST_DATA_DIR/test_genome.fa
GENOME_GTF=$TEST_DATA_DIR/test_genes.gtf
BOWTIE2_INDEX=  # Not used for minimap2
GENOME_SIZES=$TEST_DATA_DIR/test_genome.sizes
EOF
    
    # Create genome sizes file
    cut -f1,2 "$TEST_DATA_DIR/test_genome.fa.fai" > "$TEST_DATA_DIR/test_genome.sizes"
    
    # Run pipeline test
    log "Running end-to-end pipeline test..."
    local pipeline_test_dir="$TEST_DATA_DIR/pipeline_test"
    
    if "$SCRIPT_DIR/rna_preprocessing_local.sh" \
        --data-dir "$TEST_DATA_DIR" \
        --output-dir "$pipeline_test_dir" \
        --genome-config "$test_config" \
        --aligner minimap2 \
        --threads 1 \
        --memory 2 \
        > /dev/null 2>&1; then
        test_passed "End-to-end pipeline"
        
        # Check output files
        if [[ -f "$pipeline_test_dir/counts/combined_counts.txt" ]]; then
            test_passed "Count matrix generation"
        else
            test_failed "Count matrix generation"
        fi
        
        if [[ -f "$pipeline_test_dir/qc_reports/multiqc_report.html" ]]; then
            test_passed "MultiQC report generation"
        else
            test_failed "MultiQC report generation"
        fi
        
    else
        test_failed "End-to-end pipeline"
    fi
}

# Performance benchmark
run_benchmark() {
    if [[ "$QUICK_TEST" == "false" ]]; then
        header "Performance Benchmark"
        
        log "Running performance benchmark..."
        local start_time end_time duration
        start_time=$(date +%s)
        
        # Run a small benchmark
        test_pipeline > /dev/null 2>&1
        
        end_time=$(date +%s)
        duration=$((end_time - start_time))
        
        log "Benchmark completed in ${duration} seconds"
        
        # Memory usage test
        local memory_usage
        memory_usage=$(get_memory_usage)
        log "Current memory usage: ${memory_usage}GB"
    fi
}

# Cleanup test files
cleanup_test_files() {
    if [[ "$CLEANUP" == "true" && -d "$TEST_DATA_DIR" ]]; then
        log "Cleaning up test files..."
        rm -rf "$TEST_DATA_DIR"
        success "Test files cleaned up"
    fi
}

# Generate test report
generate_test_report() {
    header "Test Summary"
    
    local total_tests=$((TESTS_PASSED + TESTS_FAILED))
    local success_rate=0
    
    if [[ $total_tests -gt 0 ]]; then
        success_rate=$((TESTS_PASSED * 100 / total_tests))
    fi
    
    echo ""
    echo "Tests run: $total_tests"
    echo "Passed: $TESTS_PASSED"
    echo "Failed: $TESTS_FAILED"
    echo "Success rate: ${success_rate}%"
    echo ""
    
    if [[ $TESTS_FAILED -eq 0 ]]; then
        success "🎉 All tests passed!"
        echo ""
        echo "The RNA-seq pipeline is ready for use!"
        echo ""
        echo "Next steps:"
        echo "1. Set up reference genomes: ./setup_references.sh --genome GRCh38"
        echo "2. Run the pipeline: ./rna_preprocessing_local.sh -d your_data/"
        return 0
    else
        error "❌ Some tests failed"
        echo ""
        echo "Failed tests:"
        for test in "${FAILED_TESTS[@]}"; do
            echo "  - $test"
        done
        echo ""
        echo "Please check the installation and try again."
        return 1
    fi
}

# Main test function
main() {
    echo ""
    header "RNA-seq Pipeline Testing Suite"
    echo ""
    
    log "Test configuration:"
    log "  Test type: $TEST_TYPE"
    log "  Test data directory: $TEST_DATA_DIR"
    log "  Generate test data: $GENERATE_TEST_DATA"
    log "  Quick test: $QUICK_TEST"
    log "  Verbose: $VERBOSE"
    log "  Cleanup: $CLEANUP"
    echo ""
    
    # Run tests based on type
    case "$TEST_TYPE" in
        "env")
            test_environment
            ;;
        "modules")
            test_modules
            ;;
        "pipeline")
            test_pipeline
            ;;
        "all")
            test_environment
            test_modules
            test_pipeline
            run_benchmark
            ;;
        *)
            error "Unknown test type: $TEST_TYPE"
            exit 1
            ;;
    esac
    
    # Generate report
    if generate_test_report; then
        cleanup_test_files
        exit 0
    else
        cleanup_test_files
        exit 1
    fi
}

# Run main function
main