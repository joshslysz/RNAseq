#!/bin/bash
# Example usage with this reference genome

# Activate environment
source activate_pipeline.sh

# Run pipeline with this genome
./rna_preprocessing_local.sh \
    --data-dir raw_data/ \
    --genome-config /home/joshslysz/projects/RNAseqpipe/references/GRCr8/config.txt \
    --output-dir results_grcr8/

# Or specify paths directly
./rna_preprocessing_local.sh \
    --data-dir raw_data/ \
    --genome-fasta /home/joshslysz/projects/RNAseqpipe/references/GRCr8/genome.fa \
    --genome-gtf /home/joshslysz/projects/RNAseqpipe/references/GRCr8/genes.gtf \
    --hisat2-index /home/joshslysz/projects/RNAseqpipe/references/GRCr8/hisat2_index/GRCr8
