#!/bin/bash
# Example usage with this reference genome

# Activate environment
source activate_pipeline.sh

# Run pipeline with this genome
./rna_preprocessing_local.sh \
    --data-dir raw_data/ \
    --genome-config /home/joshslysz/projects/RNAseqpipe/references/GRCh38/config.txt \
    --output-dir results_grch38/

# Or specify paths directly
./rna_preprocessing_local.sh \
    --data-dir raw_data/ \
    --genome-fasta /home/joshslysz/projects/RNAseqpipe/references/GRCh38/genome.fa \
    --genome-gtf /home/joshslysz/projects/RNAseqpipe/references/GRCh38/genes.gtf \
    --hisat2-index /home/joshslysz/projects/RNAseqpipe/references/GRCh38/hisat2_index/GRCh38
