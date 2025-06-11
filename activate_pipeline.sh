#!/bin/bash
# RNA-seq Pipeline Environment Activation Script
# Author: Joshua Slysz, PhD - Dalhousie University
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
