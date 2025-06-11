#!/bin/bash
# Test RNA-seq pipeline environment
# Author: Joshua Slysz, PhD - Dalhousie University
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
