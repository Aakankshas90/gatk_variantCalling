#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# directories
vcf_files="/Users/balgovindyadav/Downloads/Aakanksha/snp/results/varscan_snp"
file_extension=".vcf"
results="/Users/balgovindyadav/Downloads/Aakanksha/snp/results/varscan_snp/annotated"
snpEff="/Users/balgovindyadav/Downloads/Aakanksha/snp/results/snpEff"

# Loop through all files with the specified extension
for vcf_file in "${vcf_files}"/*"${file_extension}"; do
    # Check if the file exists
    if [ -e "$vcf_file" ]; then
        # Extract sample name without extension
        sample_name=$(basename "$vcf_file" $file_extension)

	# annotate
	java -Xmx8g -jar ${snpEff}/snpEff.jar \
	Arabidopsis_thaliana "${vcf_file}" \
	> ${results}/"${sample_name}"_annot.vcf \
	-stats ${results}/"${sample_name}"_annot

        echo "Processed sample: ${vcf_file}"
    else
        echo "File not found: ${vcf_file}"
    fi
done
