#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Base Quality Score Recalibration (BQSR) using GATK
# This step corrects systematic errors made by the sequencing machine
# using a set of known variants

# directories
ref="/RAID1/working/A423/aakanksha/SNPcalling/reference_genome/TAIR10_chr_all.fas"
known_sites="/RAID1/working/A423/aakanksha/SNPcalling/knownvcf_ensembl58/arabidopsis_thaliana.vcf"
aligned="/RAID1/working/A423/aakanksha/SNPcalling/aligned"
results="/RAID1/working/A423/aakanksha/SNPcalling/results"

file_extension="_markdup.bam"

# Loop through all BAM files after duplicate marking
for bam_file in "${aligned}"/*"${file_extension}"; do

    # Extract sample name
    sample_name=$(basename "$bam_file" $file_extension)

    echo "Running BQSR for: ${sample_name}"

    # Step 1: Build recalibration model
    # GATK uses known variant sites to detect systematic errors in base quality scores
    gatk-4.5.0.0/gatk BaseRecalibrator \
        -I "$bam_file" \
        -R "$ref" \
        --known-sites "$known_sites" \
        -O "${results}/${sample_name}_recal.table"

    # Step 2: Apply recalibration
    # Adjust base quality scores based on the model generated above
    gatk-4.5.0.0/gatk ApplyBQSR \
        -I "$bam_file" \
        -R "$ref" \
        --bqsr-recal-file "${results}/${sample_name}_recal.table" \
        -O "${results}/${sample_name}_bqsr.bam"

    echo "Processed sample: ${sample_name}"

done

# Note:
# Accuracy of BQSR depends on quality of known variant sites
