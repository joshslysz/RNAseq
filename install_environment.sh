#!/bin/bash

# RNA-seq Pipeline Environment Installation Script
# Author: Joshua Slysz, PhD - Dalhousie University
# Automated conda environment setup for memory-optimized RNA-seq pipeline
# Compatible with conda/mamba package managers

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Check if script is in correct directory
if [[ ! -f "environment.yml" ]]; then
    error "environment.yml not found in current directory!"
    error "Please run this script from the RNAseqpipe directory."
    exit 1
fi

log "Setting up RNA-seq pipeline conda environment..."

# Detect conda/mamba installation
CONDA_CMD=""
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    log "Using mamba for faster installation..."
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    log "Using conda for installation..."
else
    error "Neither conda nor mamba found!"
    error "Please install Miniconda or Anaconda first:"
    error "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Check conda version
CONDA_VERSION=$($CONDA_CMD --version | awk '{print $2}')
log "Found $CONDA_CMD version: $CONDA_VERSION"

# Remove existing environment if it exists
if $CONDA_CMD env list | grep -q "rnaseqpipe"; then
    warning "Existing 'rnaseqpipe' environment found."
    read -p "Remove existing environment? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        log "Removing existing rnaseqpipe environment..."
        $CONDA_CMD env remove -n rnaseqpipe -y
    else
        error "Cannot proceed with existing environment. Exiting."
        exit 1
    fi
fi

# Create environment
log "Creating rnaseqpipe environment from environment.yml..."
if $CONDA_CMD env create -f environment.yml; then
    success "Environment created successfully!"
else
    error "Failed to create conda environment!"
    error "Try running: $CONDA_CMD clean --all"
    error "Then re-run this script."
    exit 1
fi

# Activate environment for verification
log "Activating environment for verification..."
source $(conda info --base)/etc/profile.d/conda.sh
if conda activate rnaseqpipe; then
    success "Environment activated successfully!"
else
    error "Failed to activate environment!"
    exit 1
fi

# Verify critical tool installations
log "Verifying tool installations..."

# Function to check tool installation
check_tool() {
    local tool=$1
    local version_cmd=$2
    
    if command -v "$tool" &> /dev/null; then
        local version_output
        version_output=$(eval "$version_cmd" 2>&1 | head -1 || echo "Version check failed")
        success "$tool: $version_output"
        return 0
    else
        error "$tool: NOT FOUND"
        return 1
    fi
}

# Check essential tools
TOOLS_OK=true

check_tool "fastqc" "fastqc --version" || TOOLS_OK=false
check_tool "bowtie2" "bowtie2 --version | head -1" || TOOLS_OK=false
check_tool "minimap2" "minimap2 --version" || TOOLS_OK=false
check_tool "samtools" "samtools --version | head -1" || TOOLS_OK=false
check_tool "featureCounts" "featureCounts -v 2>&1 | head -1" || TOOLS_OK=false
check_tool "salmon" "salmon --version | head -1" || TOOLS_OK=false
check_tool "multiqc" "multiqc --version" || TOOLS_OK=false
check_tool "python" "python --version" || TOOLS_OK=false

# Check Python packages
log "Checking Python packages..."
python -c "import pandas, numpy, matplotlib, seaborn, pysam, psutil; print('All Python packages imported successfully')" || TOOLS_OK=false

if [[ "$TOOLS_OK" == "true" ]]; then
    success "All tools verified successfully!"
else
    error "Some tools failed verification. Check the output above."
    exit 1
fi

# Create directories
log "Creating project directories..."
mkdir -p {references,results,logs,config,scripts,test_data}
success "Project directories created!"

# Create activation script
log "Creating environment activation script..."
cat > activate_pipeline.sh << 'EOF'
#!/bin/bash
# Activate RNA-seq pipeline environment
# Source this script before running the pipeline

source $(conda info --base)/etc/profile.d/conda.sh
conda activate rnaseqpipe

echo "RNA-seq pipeline environment activated!"
echo "Python version: $(python --version)"
echo "Available tools:"
echo "  - FastQC: $(fastqc --version)"
echo "  - Bowtie2: $(bowtie2 --version | head -1)"
echo "  - Samtools: $(samtools --version | head -1)"
echo ""
echo "Ready to run pipeline!"
echo "Usage: ./rna_preprocessing_local.sh -d <data_directory>"
EOF

chmod +x activate_pipeline.sh
success "Environment activation script created: activate_pipeline.sh"

# Create test script
log "Creating environment test script..."
cat > test_environment.sh << 'EOF'
#!/bin/bash
# Test RNA-seq pipeline environment
# Run this script to verify the environment works correctly

set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate rnaseqpipe

echo "Testing RNA-seq pipeline environment..."
echo "=================================="

# Test each tool
echo "Testing FastQC..."
fastqc --help > /dev/null && echo "✓ FastQC OK"

echo "Testing Bowtie2..."
bowtie2 --help > /dev/null && echo "✓ Bowtie2 OK"

echo "Testing minimap2..."
minimap2 --help > /dev/null && echo "✓ minimap2 OK"

echo "Testing Samtools..."
samtools --help > /dev/null && echo "✓ Samtools OK"

echo "Testing featureCounts..."
featureCounts 2>&1 | head -1 > /dev/null && echo "✓ featureCounts OK"

echo "Testing Salmon..."
salmon --help > /dev/null && echo "✓ Salmon OK"

echo "Testing MultiQC..."
multiqc --help > /dev/null && echo "✓ MultiQC OK"

echo "Testing Python packages..."
python -c "
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import psutil
print('✓ All Python packages OK')
"

echo ""
echo "🎉 All tests passed! Environment is ready for RNA-seq analysis."
EOF

chmod +x test_environment.sh
success "Environment test script created: test_environment.sh"

# Final instructions
echo ""
success "🎉 RNA-seq pipeline environment setup complete!"
echo ""
echo "Next steps:"
echo "==========="
echo "1. Test the environment:"
echo "   ./test_environment.sh"
echo ""
echo "2. Activate the environment:"
echo "   source activate_pipeline.sh"
echo "   # OR manually:"
echo "   conda activate rnaseqpipe"
echo ""
echo "3. Download reference genomes:"
echo "   ./setup_references.sh --genome GRCh38"
echo ""
echo "4. Run the pipeline:"
echo "   ./rna_preprocessing_local.sh -d your_data_directory/"
echo ""
echo "For help:"
echo "   ./rna_preprocessing_local.sh --help"
echo ""
warning "Remember to activate the environment before running any pipeline commands!"