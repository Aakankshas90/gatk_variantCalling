#!/bin/bash

# directories
ref="/RAID1/working/A423/aakanksha/SNPcalling/reference_genome/TAIR10_chr_all.fas"
known_sites="/RAID1/working/A423/aakanksha/SNPcalling/knownvcf_ensembl58/arabidopsis_thaliana.vcf"
aligned_reads="/RAID1/working/A423/aakanksha/SNPcalling/aligned"
results="/RAID1/working/A423/aakanksha/SNPcalling/results"
file_extension=".bam"

# Loop through all files with the specified extension
for bam_file in "${aligned_reads}"/*"${file_extension}"; do
    # Check if the file exists
    if [ -e "$bam_file" ]; then
        # Extract sample name without extension
        sample_name=$(basename "$bam_file" $file_extension)

        # 1. build the model
        gatk-4.5.0.0/gatk BaseRecalibrator -I $bam_file -R ${ref} --known-sites ${known_sites} -O ${aligned_reads}/${sample_name}_recal.table

        # 2. Apply the model to adjust the base quality scores
        gatk-4.5.0.0/gatk ApplyBQSR -I $bam_file -R ${ref} --bqsr-recal-file ${aligned_reads}/${sample_name}_recal.table -O ${results}/${sample_name}_bqsr.bam 

        echo "Processed sample: ${sample_name}"
    else
        echo "File not found: ${bam_file}"
    fi
done
