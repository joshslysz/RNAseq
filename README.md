# RNA-seq Pipeline

**Author:** Joshua Slysz, PhD  
**Institution:** Dalhousie University

This pipeline takes your raw RNA-seq data and performs a complete RNA-seq analysis workflow:
1. **Quality Control** - FastQC analysis of raw reads, adapter detection, quality assessment
2. **Preprocessing** - Adapter trimming, quality filtering, read length optimization
3. **Splice-Aware Alignment** - Maps reads to genome using RNA-seq specific aligners (HISAT2, STAR, or minimap2)
4. **Gene Quantification** - Counts reads mapping to genes using featureCounts or Salmon
5. **Quality Assessment** - Comprehensive QC reports with MultiQC, alignment statistics
6. **Results Integration** - Combines counts across samples into analysis-ready matrices


## Input Data Requirements

### File Naming Convention

Your FASTQ files must follow these exact naming patterns:

**For paired-end sequencing (most common):**
```
sample1_1.fq.gz    # Read 1 for sample1
sample1_2.fq.gz    # Read 2 for sample1
sample2_1.fq.gz    # Read 1 for sample2
sample2_2.fq.gz    # Read 2 for sample2
```

**Important:**
- Files must end with `.fq.gz` (compressed FASTQ)
- Use descriptive sample names (e.g., `control_rep1`, `treatment_24h`)
- Avoid spaces and special characters in sample names
- Both read files for paired-end samples must have identical names except for `_1` vs `_2`

### Example File Structure
```
raw_data/
├── control_rep1_1.fq.gz
├── control_rep1_2.fq.gz
├── control_rep2_1.fq.gz
├── control_rep2_2.fq.gz
├── treatment_rep1_1.fq.gz
├── treatment_rep1_2.fq.gz
├── treatment_rep2_1.fq.gz
└── treatment_rep2_2.fq.gz
```

## 🛠️ Installation & Setup

### Step 1: Prerequisites

**Install Conda/Miniconda** (if you don't have it):
```bash
# Download Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install (follow the prompts)
bash Miniconda3-latest-Linux-x86_64.sh

# Restart your terminal or run:
source ~/.bashrc
```

### Step 2: Download the Pipeline

```bash
# Download the pipeline
git clone https://github.com/joshslysz/RNAseq.git
cd RNAseq

# Make sure you're in the right directory
ls  # You should see install_environment.sh
```

### Step 3: Install Required Software

```bash
# Run the automatic installer
./install_environment.sh

# This will take 10-30 minutes depending on your internet speed
# The installer will:
# - Create a conda environment with all required tools
# - Test that everything works
# - Create helpful scripts
```

### Step 4: Test Installation

```bash
# Test that everything is working
./test_pipeline.sh

# If you see "🎉 All tests passed!" you're ready to go!
```

### Step 5: Download Reference Genome

```bash
# Activate the environment
conda activate rnaseqpipe

# For human samples
./setup_references.sh --genome GRCh38

# For mouse samples
./setup_references.sh --genome GRCm39

# For rat samples
./setup_references.sh --genome GRCr8

# This downloads ~3-4GB of reference files
# Takes 15-30 minutes depending on internet speed
```

## Running the Pipeline

### Basic Usage (Most Common)

```bash
# 1. Activate the environment (do this every time)
conda activate rnaseqpipe

# 2. Run the pipeline (optimized for RNA-seq!)
# For human samples:
./rna_preprocessing_local.sh \
    -d raw_data/ \
    --genome-config references/GRCh38/config.txt

# For mouse samples:
./rna_preprocessing_local.sh \
    -d raw_data/ \
    --genome-config references/GRCm39/config.txt

# For rat samples:
./rna_preprocessing_local.sh \
    -d raw_data/ \
    --genome-config references/GRCr8/config.txt

# That's it! The pipeline will:
# - Automatically trim adapters and filter low-quality reads (fastp/trimmomatic)
# - Use splice-aware alignment (HISAT2/STAR/minimap2) optimized for RNA-seq
# - Perform gene-level quantification with featureCounts
# - Generate comprehensive quality reports with FastQC and MultiQC
# - Create a results/ folder with ready-to-analyze count matrices
```

## Understanding The Results

After the pipeline finishes, you'll find your results in the `results/` folder:

```
results/
├── aligned/                 # Alignment files
├── counts/                  # Gene expression data (MAIN OUTPUT)
├── qc_reports/             # Quality control reports
├── logs/                   # Processing logs
└── pipeline_summary.txt    # Summary of what happened
```

### Key Output Files

#### 1. Gene Expression Data (`counts/`)
- **`combined_counts.txt`** - Main result file with gene expression for all samples
- **`sample1.counts.txt`** - Individual sample gene counts
- **Format:** Tab-separated with genes as rows, samples as columns

#### 2. Quality Control Reports (`qc_reports/`)
- **`multiqc_report.html`** - Open this in your web browser for a comprehensive quality report
- **Shows:** Sequence quality, alignment rates, gene detection rates

#### 3. Alignment Files (`aligned/`)
- **`sample1.sorted.bam`** - Aligned sequences (for advanced users)
- **`sample1.flagstat`** - Alignment statistics

### Reading the Gene Count File

Your main result (`combined_counts.txt`) looks like this:
```
gene_id          control_rep1  control_rep2  treatment_rep1  treatment_rep2
ENSG00000000003  742           856           423             445
ENSG00000000005  0             2             1               0
ENSG00000000419  1247          1356          2134            2445
```

Each number represents how many RNA molecules from that gene were detected in that sample.

