#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# directories
ref="/RAID1/working/A423/aakanksha/SNPcalling/reference_genome/TAIR10_chr_all.fas"
vcf_files="/RAID1/working/A423/aakanksha/SNPcalling/results/vcf"
file_extension=".vcf"
results="/RAID1/working/A423/aakanksha/SNPcalling/results/vcf/split"

# Loop through all files with the specified extension
for vcf_file in "${vcf_files}"/*"${file_extension}"; do
    # Check if the file exists
    if [ -e "$vcf_file" ]; then
        # Extract sample name without extension
        sample_name=$(basename "$vcf_file" $file_extension)

        # 1. SelectVariants SNP
        gatk-4.5.0.0/gatk SelectVariants -R ${ref} -V "$vcf_file" --select-type-to-include SNP -O ${results}/"${sample_name}"_snp.vcf.gz 2> ${results}/"${sample_name}"_snp.txt

        # 1. SelectVariants INDELS
        gatk-4.5.0.0/gatk SelectVariants -R ${ref} -V "$vcf_file" --select-type-to-include INDEL -O ${results}/"${sample_name}"_indel.vcf.gz 2> ${results}/"${sample_name}"_indel.txt

        echo "Processed sample: ${sample_name}"
    else
        echo "File not found: ${vcf_file}"
    fi
done
