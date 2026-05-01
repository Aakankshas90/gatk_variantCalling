#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Variant calling using GATK HaplotypeCaller
# This step identifies SNPs and indels from recalibrated BAM files
# Output is generated in GVCF format for downstream joint genotyping

# directories
ref="/RAID1/working/A423/aakanksha/SNPcalling/reference_genome/TAIR10_chr_all.fas"
results="/RAID1/working/A423/aakanksha/SNPcalling/results"

file_extension="_bqsr.bam"

# Loop through all recalibrated BAM files
for bam_file in "${results}"/*"${file_extension}"; do

    # Extract sample name
    sample_name=$(basename "$bam_file" $file_extension)

    echo "Calling variants for: ${sample_name}"

    # Perform variant calling in GVCF mode
    # -ERC GVCF allows joint genotyping across multiple samples later
    gatk-4.5.0.0/gatk HaplotypeCaller \
        -ERC GVCF \
        -R "$ref" \
        -I "$bam_file" \
        -O "${results}/${sample_name}.g.vcf"

    echo "Processed sample: ${sample_name}"

done
