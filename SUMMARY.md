# RNA-seq Pipeline Architecture Summary

## Overview

This memory-optimized RNA-seq pipeline provides a complete analysis workflow from raw FASTQ files to gene expression count matrices. The pipeline is designed to run efficiently on local systems with limited computational resources (2-32GB RAM) while maintaining professional-quality results.

## High-Level Architecture

The pipeline follows a modular design with intelligent resource management:

```
Raw FASTQ Files → Preprocessing → Alignment → Quantification → QC Reports → Count Matrices
       ↓               ↓            ↓             ↓            ↓            ↓
   Quality Check   Adapter Trim   RNA-seq      Gene-level   Comprehensive  Analysis-ready
   File Detection  Quality Filter  Splice-aware  Counting     QC Reports    Results
```

## Core Pipeline Components

### 1. Main Pipeline Controller (`rna_preprocessing_local.sh`)
- **Purpose**: Orchestrates the entire analysis workflow
- **Key Features**:
  - Automatic resource detection and optimization
  - Sample discovery and validation
  - Sequential processing to minimize memory usage
  - Resume capability for interrupted runs
  - Progress tracking and logging

### 2. Resource Management System (`scripts/resource_manager.sh`)
- **Purpose**: Dynamically manages system resources and optimizes tool selection
- **Intelligence**:
  - Auto-detects available memory, CPU cores, and disk space
  - Selects optimal RNA-seq aligner based on system capabilities
  - Configures thread counts and memory limits
  - Monitors resource usage during execution

### 3. Preprocessing Module (`scripts/preprocessing_module.sh`)
- **Purpose**: Prepares raw reads for alignment
- **Operations**:
  - Adapter trimming (fastp/trimmomatic/cutadapt)
  - Quality filtering and read length optimization
  - Optional rRNA removal (placeholder for future implementation)
  - Preprocessing statistics generation

### 4. RNA-seq Alignment Module (`scripts/rna_alignment_module.sh`)
- **Purpose**: Performs splice-aware alignment optimized for RNA sequencing
- **Supported Aligners**:
  - **HISAT2**: Memory-efficient splice-aware aligner (4-16GB RAM)
  - **STAR**: High-performance aligner (30GB+ RAM recommended)
  - **minimap2**: Ultra-low memory fallback (<4GB RAM)
- **Features**:
  - Automatic index building from reference genome
  - Splice site detection using gene annotations
  - Strand-specific alignment support
  - Troubleshooting mode for difficult samples

### 5. Gene Quantification Module (`scripts/quantification_module.sh`)
- **Purpose**: Counts reads mapping to genes/transcripts
- **Quantification Methods**:
  - **featureCounts**: Gene-level counting from BAM files
  - **Salmon**: Transcript-level quantification (future enhancement)
- **Optimizations**:
  - Bulk RNA-seq parameter tuning
  - Multi-mapping read handling
  - Strand-specific counting options

### 6. Quality Control Module (`scripts/qc_module.sh`)
- **Purpose**: Comprehensive quality assessment
- **QC Components**:
  - **FastQC**: Raw read quality analysis
  - **RSeQC**: RNA-seq specific quality metrics (memory permitting)
  - **MultiQC**: Integrated quality reports
- **Memory Adaptation**: Disables memory-intensive QC on low-resource systems

### 7. Legacy Alignment Module (`scripts/alignment_module.sh`)
- **Purpose**: General-purpose alignment (primarily for DNA-seq)
- **Note**: Superseded by RNA-specific alignment module for RNA-seq data

## Intelligent Tool Selection Logic

The pipeline automatically selects optimal tools based on available resources:

### Aligner Selection Strategy
```
Available Memory >= 32GB + STAR available → Use STAR
Available Memory >= 16GB + HISAT2 available → Use HISAT2  
Available Memory >= 4GB + HISAT2 available → Use HISAT2 (low-mem mode)
Available Memory < 4GB → Use minimap2
```

### Memory Optimization Features
- **Sequential Sample Processing**: Prevents memory conflicts
- **Dynamic Thread Allocation**: Scales with available CPU cores
- **Temporary File Management**: Automatic cleanup to save disk space
- **Progress Resumption**: Skip completed samples on restart

## Data Flow and Processing Steps

### Step 1: Input Validation and Setup
1. Validate FASTQ file naming conventions
2. Check system requirements (memory, disk space)
3. Load genome configuration
4. Initialize output directory structure

### Step 2: Per-Sample Processing Loop
For each sample:
1. **Preprocessing**:
   - FastQC quality assessment
   - Adapter trimming and quality filtering
   - Generate preprocessing statistics

2. **Alignment**:
   - Build/validate reference indices
   - Perform splice-aware alignment
   - Convert and sort BAM files
   - Generate alignment statistics

3. **Quantification**:
   - Count reads mapping to genes
   - Generate quantification statistics
   - Quality assessment of count data

### Step 3: Integration and Reporting
1. **Count Matrix Generation**:
   - Combine individual sample counts
   - Create analysis-ready matrices
   - Validate data integrity

2. **Quality Control**:
   - Generate comprehensive MultiQC reports
   - Assess overall pipeline performance
   - Flag potential issues

3. **Final Output**:
   - Pipeline summary report
   - Log file consolidation
   - Temporary file cleanup

## Output Structure

```
results/
├── aligned/              # BAM files and alignment statistics
├── counts/              # Gene expression count matrices
│   ├── combined_counts.txt    # Main output: all samples combined
│   └── individual_*.counts.txt # Per-sample count files
├── qc_reports/          # Quality control reports
│   ├── multiqc_report.html    # Comprehensive QC dashboard
│   ├── fastqc/               # FastQC outputs
│   └── rseqc/                # RNA-seq specific QC (if run)
├── logs/                # Processing logs and progress tracking
└── pipeline_summary.txt # Overall pipeline report
```

## Key Design Principles

### 1. Memory Efficiency
- Processes samples sequentially rather than in parallel
- Uses memory-mapped I/O where possible
- Automatic memory monitoring and adjustment
- Graceful degradation on low-memory systems

### 2. RNA-seq Optimization
- Prioritizes splice-aware aligners over DNA aligners
- Uses RNA-seq specific parameters for all tools
- Optimizes for bulk RNA-seq experimental designs
- Includes strand-specificity detection and handling

### 3. Robustness and Reliability
- Comprehensive error handling and logging
- Resume capability for long-running analyses
- Extensive input validation
- Quality checkpoints throughout the pipeline

### 4. User Accessibility
- Single-command execution for most use cases
- Automatic tool selection and parameter optimization
- Clear progress reporting and error messages
- Comprehensive documentation and troubleshooting guides

## System Adaptability

The pipeline dynamically adapts to different system configurations:

- **High-performance systems** (32GB+ RAM): Uses STAR aligner with full QC suite
- **Standard systems** (8-16GB RAM): Uses HISAT2 with comprehensive analysis
- **Limited systems** (4-8GB RAM): Uses HISAT2 in low-memory mode
- **Minimal systems** (2-4GB RAM): Uses minimap2 with essential QC only

This adaptive approach ensures the pipeline can run effectively across a wide range of hardware configurations while maximizing the quality of results given available resources.