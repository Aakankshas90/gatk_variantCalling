#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# directories
ref="/RAID1/working/A423/aakanksha/SNPcalling/reference_genome/TAIR10_chr_all.fas"
snp_files="/RAID1/working/A423/aakanksha/SNPcalling/results/combinedVCF/SNP"
indel_files="/RAID1/working/A423/aakanksha/SNPcalling/results/combinedVCF/indel"
snp_extension="_snp.vcf.gz"
indel_extension="_indels.vcf.gz"
results="/RAID1/working/A423/aakanksha/SNPcalling/results/VariantFiltration"
log_file="${results}/error_log.txt"

# Loop through all snp files with the specified extension
for snp_file in "${snp_files}"/*"${snp_extension}"; do
    # Check if the file exists
    if [ -e "$snp_file" ]; then
        # Extract sample name without extension
        snp_name=$(basename "$snp_file" $snp_extension)

        # 1.Filter SNPs
        gatk-4.5.0.0/gatk VariantFiltration \
        -R ${ref} \
        -V "${snp_file}" \
        -O ${results}/"${snp_name}"_snp_filt.vcf  \
	2>> "${log_file}" \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "QUAL30"  -filter "QUAL < 30.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
        -genotype-filter-expression "DP < 10" \
        -genotype-filter-name "DP_filter" \
        -genotype-filter-expression "GQ < 10" \
        -genotype-filter-name "GQ_filter"

	# 2.Select SNPs that PASS filters
	gatk-4.5.0.0/gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/"${snp_name}"_snp_filt.vcf \
	-O ${results}/"${snp_name}"_snp_filtPASS.vcf \
	2>> "${log_file}"

	# 3.Exclude SNPs that failed genotype filters
	cat ${results}/"${snp_name}"_snp_filtPASS.vcf|grep \
	-v -E "DP_filter|GQ_filter" \
	> ${results}/"${snp_name}"_snp_filtPASS_filtGT.vcf \
	2>> "${log_file}"
      
        echo "Processed sample: ${snp_name}"
    else
        echo "File not found: ${snp_file}"
    fi
done

# Loop through all files with the specified extension
for indel_file in "${indel_files}"/*"${indel_extension}"; do
    # Check if the file exists
    if [ -e "$indel_file" ]; then
        # Extract sample name without extension
        indel_name=$(basename "$indel_file" $indel_extension)

        # 1.Filter INDELS
        gatk-4.5.0.0/gatk VariantFiltration \
        -R ${ref} \
        -V "${indel_file}" \
        -O ${results}/"${indel_name}"_indel_filt.vcf  \
	2>> "${log_file}" \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "QUAL30"  -filter "QUAL < 30.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -genotype-filter-expression "DP < 10" \
        -genotype-filter-name "DP_filter" \
        -genotype-filter-expression "GQ < 10" \
        -genotype-filter-name "GQ_filter"

	# 2.Select Indels that PASS filters
	gatk-4.5.0.0/gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/"${indel_name}"_indel_filt.vcf \
	-O ${results}/"${indel_name}"_indel_filtPASS.vcf \
	2>> "${log_file}"

	# 3.Exclude Indels that failed genotype filters
	cat ${results}/"${indel_name}"_indel_filtPASS.vcf|grep \
	-v -E "DP_filter|GQ_filter" \
	> ${results}/"${indel_name}"_indel_filtPASS_filtGT.vcf \
	2>> "${log_file}"

        echo "Processed sample: ${indel_name}"
    else
        echo "File not found: ${indel_file}"
    fi
done
