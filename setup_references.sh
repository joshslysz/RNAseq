#!/bin/bash

# Reference Genome Setup Script
# Author: Joshua Slysz, PhD - Dalhousie University
# Downloads and builds indices for RNA-seq pipeline
# Supports GRCh38, GRCm39, and custom genomes

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Default values
GENOME=""
BUILD_INDICES=true
THREADS=$(nproc)
REFERENCE_DIR="references"
FORCE_DOWNLOAD=false
CUSTOM_FASTA=""
CUSTOM_GTF=""

# Logging functions
log() { echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }

# Usage function
usage() {
    cat << EOF
Usage: $0 --genome GENOME [options]

Download and setup reference genomes for RNA-seq pipeline

Required arguments:
  --genome GENOME         Reference genome (GRCh38, GRCm39, or custom)

Optional arguments:
  --reference-dir DIR     Reference directory [default: references]
  --threads N             Number of threads for index building [default: $(nproc)]
  --no-build              Download only, don't build indices
  --force                 Force re-download even if files exist
  --custom-fasta FILE     Custom genome FASTA file (for --genome custom)
  --custom-gtf FILE       Custom annotation GTF file (for --genome custom)
  --help                  Show this help message

Supported genomes:
  GRCh38                 Human genome (Homo sapiens)
  GRCm39                 Mouse genome (Mus musculus)
  custom                 User-provided FASTA and GTF files

Examples:
  # Download human genome
  $0 --genome GRCh38
  
  # Download mouse genome with custom directory
  $0 --genome GRCm39 --reference-dir /data/references
  
  # Use custom genome files
  $0 --genome custom --custom-fasta genome.fa --custom-gtf genes.gtf
EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --genome)
            GENOME="$2"
            shift 2
            ;;
        --reference-dir)
            REFERENCE_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --no-build)
            BUILD_INDICES=false
            shift
            ;;
        --force)
            FORCE_DOWNLOAD=true
            shift
            ;;
        --custom-fasta)
            CUSTOM_FASTA="$2"
            shift 2
            ;;
        --custom-gtf)
            CUSTOM_GTF="$2"
            shift 2
            ;;
        --help)
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

# Validate required arguments
if [[ -z "$GENOME" ]]; then
    error "Genome must be specified with --genome"
    usage
    exit 1
fi

# Validate custom genome arguments
if [[ "$GENOME" == "custom" ]]; then
    if [[ -z "$CUSTOM_FASTA" || -z "$CUSTOM_GTF" ]]; then
        error "Custom genome requires --custom-fasta and --custom-gtf"
        exit 1
    fi
    if [[ ! -f "$CUSTOM_FASTA" ]]; then
        error "Custom FASTA file not found: $CUSTOM_FASTA"
        exit 1
    fi
    if [[ ! -f "$CUSTOM_GTF" ]]; then
        error "Custom GTF file not found: $CUSTOM_GTF"
        exit 1
    fi
fi

# Check if conda environment is activated
if [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "$CONDA_DEFAULT_ENV" != "rnaseqpipe" ]]; then
    error "Please activate the rnaseqpipe conda environment first:"
    error "  conda activate rnaseqpipe"
    error "  OR source activate_pipeline.sh"
    exit 1
fi

# Check required tools
check_tool() {
    if ! command -v "$1" &> /dev/null; then
        error "Required tool '$1' not found. Please check conda environment."
        exit 1
    fi
}

log "Checking required tools..."
check_tool hisat2-build
check_tool samtools
check_tool wget
success "All required tools found!"

# Create reference directory
mkdir -p "$REFERENCE_DIR"
REFERENCE_DIR=$(realpath "$REFERENCE_DIR")
log "Using reference directory: $REFERENCE_DIR"

# Genome-specific configurations
case "$GENOME" in
    GRCh38)
        GENOME_NAME="GRCh38"
        FASTA_URL="http://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        GTF_URL="http://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz"
        ;;
    GRCm39)
        GENOME_NAME="GRCm39"
        FASTA_URL="http://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
        GTF_URL="http://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.gtf.gz"
        ;;
    custom)
        GENOME_NAME="custom"
        ;;
    *)
        error "Unsupported genome: $GENOME"
        error "Supported genomes: GRCh38, GRCm39, custom"
        exit 1
        ;;
esac

# Create genome-specific directory
GENOME_DIR="$REFERENCE_DIR/$GENOME_NAME"
mkdir -p "$GENOME_DIR"

# Function to download file
download_file() {
    local url=$1
    local output=$2
    local description=$3
    
    if [[ -f "$output" && "$FORCE_DOWNLOAD" == "false" ]]; then
        log "$description already exists: $output"
        return 0
    fi
    
    log "Downloading $description..."
    if wget -O "$output" "$url"; then
        success "Downloaded: $output"
    else
        error "Failed to download: $url"
        return 1
    fi
}

# Function to check memory requirements
check_memory() {
    local available_mem_gb
    available_mem_gb=$(free -g | awk 'NR==2{print $7}')
    
    if [[ $available_mem_gb -lt 4 ]]; then
        warning "Available memory: ${available_mem_gb}GB"
        warning "HISAT2 index building may require 4-8GB RAM"
        warning "Consider using a machine with more memory or pre-built indices"
        
        read -p "Continue anyway? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            log "Exiting due to memory constraints"
            exit 1
        fi
    else
        log "Available memory: ${available_mem_gb}GB (sufficient for index building)"
    fi
}

# Download and setup genome files
if [[ "$GENOME" != "custom" ]]; then
    # Download FASTA and GTF files
    FASTA_FILE="$GENOME_DIR/genome.fa.gz"
    GTF_FILE="$GENOME_DIR/genes.gtf.gz"
    
    download_file "$FASTA_URL" "$FASTA_FILE" "genome FASTA"
    download_file "$GTF_URL" "$GTF_FILE" "gene annotations GTF"
    
    # Decompress files
    log "Decompressing genome files..."
    if [[ ! -f "$GENOME_DIR/genome.fa" || "$FORCE_DOWNLOAD" == "true" ]]; then
        gunzip -c "$FASTA_FILE" > "$GENOME_DIR/genome.fa"
        success "Decompressed genome FASTA"
    fi
    
    if [[ ! -f "$GENOME_DIR/genes.gtf" || "$FORCE_DOWNLOAD" == "true" ]]; then
        gunzip -c "$GTF_FILE" > "$GENOME_DIR/genes.gtf"
        success "Decompressed gene annotations"
    fi
    
    GENOME_FASTA="$GENOME_DIR/genome.fa"
    GENOME_GTF="$GENOME_DIR/genes.gtf"
else
    # Copy custom files
    log "Setting up custom genome files..."
    cp "$CUSTOM_FASTA" "$GENOME_DIR/genome.fa"
    cp "$CUSTOM_GTF" "$GENOME_DIR/genes.gtf"
    success "Custom genome files copied"
    
    GENOME_FASTA="$GENOME_DIR/genome.fa"
    GENOME_GTF="$GENOME_DIR/genes.gtf"
fi

# Build indices if requested
if [[ "$BUILD_INDICES" == "true" ]]; then
    log "Building HISAT2 indices..."
    
    # Check memory requirements
    check_memory
    
    # Create index directory
    INDEX_DIR="$GENOME_DIR/hisat2_index"
    mkdir -p "$INDEX_DIR"
    
    # Build HISAT2 index
    INDEX_PREFIX="$INDEX_DIR/$GENOME_NAME"
    
    if [[ ! -f "$INDEX_PREFIX.1.ht2" || "$FORCE_DOWNLOAD" == "true" ]]; then
        log "Building HISAT2 index with $THREADS threads..."
        
        # Memory-optimized HISAT2 build
        if hisat2-build \
            -p "$THREADS" \
            "$GENOME_FASTA" \
            "$INDEX_PREFIX"; then
            success "HISAT2 index built successfully!"
        else
            error "Failed to build HISAT2 index"
            exit 1
        fi
    else
        log "HISAT2 index already exists"
    fi
    
    # Create FASTA index for samtools
    log "Creating FASTA index..."
    if [[ ! -f "$GENOME_FASTA.fai" || "$FORCE_DOWNLOAD" == "true" ]]; then
        samtools faidx "$GENOME_FASTA"
        success "FASTA index created"
    else
        log "FASTA index already exists"
    fi
    
    # Create genome size file for later use
    log "Creating genome size file..."
    cut -f1,2 "$GENOME_FASTA.fai" > "$GENOME_DIR/genome.sizes"
    success "Genome size file created"
fi

# Create reference configuration file
log "Creating reference configuration..."
cat > "$GENOME_DIR/config.txt" << EOF
# Reference genome configuration
# Generated on $(date)

GENOME_NAME=$GENOME_NAME
GENOME_FASTA=$GENOME_FASTA
GENOME_GTF=$GENOME_GTF
HISAT2_INDEX=$INDEX_PREFIX
GENOME_SIZES=$GENOME_DIR/genome.sizes

# Memory recommendations
HISAT2_MEMORY_GB=4
MINIMAP2_MEMORY_GB=2
FEATURE_COUNTS_MEMORY_GB=2
EOF

success "Reference configuration saved: $GENOME_DIR/config.txt"

# Create usage script
cat > "$GENOME_DIR/usage_example.sh" << EOF
#!/bin/bash
# Example usage with this reference genome

# Activate environment
source activate_pipeline.sh

# Run pipeline with this genome
./rna_preprocessing_local.sh \\
    --data-dir raw_data/ \\
    --genome-config $GENOME_DIR/config.txt \\
    --output-dir results_${GENOME_NAME,,}/

# Or specify paths directly
./rna_preprocessing_local.sh \\
    --data-dir raw_data/ \\
    --genome-fasta $GENOME_FASTA \\
    --genome-gtf $GENOME_GTF \\
    --hisat2-index $INDEX_PREFIX
EOF

chmod +x "$GENOME_DIR/usage_example.sh"

# Print summary
echo ""
success "🎉 Reference genome setup complete!"
echo ""
echo "Genome: $GENOME_NAME"
echo "Location: $GENOME_DIR"
echo ""
echo "Files created:"
echo "  📄 Genome FASTA: $GENOME_FASTA"
echo "  📄 Gene annotations: $GENOME_GTF"
if [[ "$BUILD_INDICES" == "true" ]]; then
echo "  🗂️  HISAT2 index: $INDEX_PREFIX.*"
echo "  📄 FASTA index: $GENOME_FASTA.fai"
echo "  📄 Genome sizes: $GENOME_DIR/genome.sizes"
fi
echo "  ⚙️  Configuration: $GENOME_DIR/config.txt"
echo "  📝 Usage example: $GENOME_DIR/usage_example.sh"
echo ""
echo "Next steps:"
echo "1. Run the pipeline:"
echo "   ./rna_preprocessing_local.sh --data-dir your_data/ --genome-config $GENOME_DIR/config.txt"
echo ""
echo "2. Or check the usage example:"
echo "   cat $GENOME_DIR/usage_example.sh"
echo ""

# Disk usage summary
GENOME_SIZE=$(du -sh "$GENOME_DIR" | cut -f1)
warning "Disk space used: $GENOME_SIZE"