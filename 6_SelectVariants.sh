#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Split combined VCF into SNPs and INDELs using GATK SelectVariants

# directories
ref="/RAID1/working/A423/aakanksha/SNPcalling/reference_genome/TAIR10_chr_all.fas"
vcf_files="/RAID1/working/A423/aakanksha/SNPcalling/results/vcf"
results="/RAID1/working/A423/aakanksha/SNPcalling/results/vcf/split"

file_extension=".vcf"

# Loop through all VCF files
for vcf_file in "${vcf_files}"/*"${file_extension}"; do

    if [ -e "$vcf_file" ]; then

        # Extract sample name
        sample_name=$(basename "$vcf_file" $file_extension)

        echo "Splitting variants for: ${sample_name}"

        # 1. Extract SNPs
        gatk-4.5.0.0/gatk SelectVariants \
            -R ${ref} \
            -V "$vcf_file" \
            --select-type-to-include SNP \
            -O ${results}/${sample_name}_snp.vcf.gz \
            2> ${results}/${sample_name}_snp.log

        # 2. Extract INDELs
        gatk-4.5.0.0/gatk SelectVariants \
            -R ${ref} \
            -V "$vcf_file" \
            --select-type-to-include INDEL \
            -O ${results}/${sample_name}_indel.vcf.gz \
            2> ${results}/${sample_name}_indel.log

        echo "Processed sample: ${sample_name}"

    else
        echo "File not found: ${vcf_file}"
    fi

done
