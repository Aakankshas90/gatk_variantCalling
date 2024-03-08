#!/bin/bash

# directories
ref="/RAID1/working/A423/aakanksha/SNPcalling/reference_genome/TAIR10_chr_all.fas"
results="/RAID1/working/A423/aakanksha/SNPcalling/results"
file_extension="_gatk_bqsr.bam"

# Loop through all files with the specified extension
for bam_file in "${results}"/*"${file_extension}"; do
    # Check if the file exists
    if [ -e "$bam_file" ]; then
        # Extract sample name without extension
        sample_name=$(basename "$bam_file" $file_extension)

        # variant calling using gatk HaplotypeCaller
        gatk-4.5.0.0/gatk HaplotypeCaller -R ${ref} -I $bam_file -O ${results}/${sample_name}.vcf

        echo "Processed sample: ${sample_name}"
    else
        echo "File not found: ${bam_file}"
    fi
done
