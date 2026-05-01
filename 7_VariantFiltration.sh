#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Variant filtering using GATK VariantFiltration
# Applies hard filters to SNPs and INDELs separately
# Followed by selection of PASS variants and removal of low-quality genotypes

# directories
ref="/RAID1/working/A423/aakanksha/SNPcalling/reference_genome/TAIR10_chr_all.fas"
snp_files="/RAID1/working/A423/aakanksha/SNPcalling/results/combinedVCF/SNP"
indel_files="/RAID1/working/A423/aakanksha/SNPcalling/results/combinedVCF/indel"
results="/RAID1/working/A423/aakanksha/SNPcalling/results/VariantFiltration"

log_file="${results}/error_log.txt"

snp_extension="_snp.vcf.gz"
indel_extension="_indel.vcf.gz"

########################################
# SNP FILTERING
########################################

for snp_file in "${snp_files}"/*"${snp_extension}"; do

    if [ -e "$snp_file" ]; then

        snp_name=$(basename "$snp_file" $snp_extension)

        echo "Filtering SNPs for: ${snp_name}"

        # 1. Apply SNP filters
        gatk-4.5.0.0/gatk VariantFiltration \
            -R ${ref} \
            -V "${snp_file}" \
            -O ${results}/${snp_name}_snp_filt.vcf \
            2>> "${log_file}" \
            -filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "QUAL30" -filter "QUAL < 30.0" \
            -filter-name "FS_filter" -filter "FS > 60.0" \
            -filter-name "MQ_filter" -filter "MQ < 40.0" \
            -filter-name "SOR_filter" -filter "SOR > 4.0" \
            -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
            -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
            -genotype-filter-expression "DP < 10" \
            -genotype-filter-name "DP_filter" \
            -genotype-filter-expression "GQ < 10" \
            -genotype-filter-name "GQ_filter"

        # 2. Select only PASS variants
        gatk-4.5.0.0/gatk SelectVariants \
            --exclude-filtered \
            -V ${results}/${snp_name}_snp_filt.vcf \
            -O ${results}/${snp_name}_snp_filtPASS.vcf \
            2>> "${log_file}"

        # 3. Remove variants failing genotype filters
        cat ${results}/${snp_name}_snp_filtPASS.vcf | grep -v -E "DP_filter|GQ_filter" \
            > ${results}/${snp_name}_snp_filtPASS_filtGT.vcf \
            2>> "${log_file}"

        echo "Processed SNPs: ${snp_name}"

    else
        echo "File not found: ${snp_file}"
    fi

done

########################################
# INDEL FILTERING
########################################

for indel_file in "${indel_files}"/*"${indel_extension}"; do

    if [ -e "$indel_file" ]; then

        indel_name=$(basename "$indel_file" $indel_extension)

        echo "Filtering INDELs for: ${indel_name}"

        # 1. Apply INDEL filters
        gatk-4.5.0.0/gatk VariantFiltration \
            -R ${ref} \
            -V "${indel_file}" \
            -O ${results}/${indel_name}_indel_filt.vcf \
            2>> "${log_file}" \
            -filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "QUAL30" -filter "QUAL < 30.0" \
            -filter-name "FS_filter" -filter "FS > 60.0" \
            -filter-name "SOR_filter" -filter "SOR > 4.0" \
            -genotype-filter-expression "DP < 10" \
            -genotype-filter-name "DP_filter" \
            -genotype-filter-expression "GQ < 10" \
            -genotype-filter-name "GQ_filter"

        # 2. Select only PASS INDELs
        gatk-4.5.0.0/gatk SelectVariants \
            --exclude-filtered \
            -V ${results}/${indel_name}_indel_filt.vcf \
            -O ${results}/${indel_name}_indel_filtPASS.vcf \
            2>> "${log_file}"

        # 3. Remove variants failing genotype filters
        cat ${results}/${indel_name}_indel_filtPASS.vcf | grep -v -E "DP_filter|GQ_filter" \
            > ${results}/${indel_name}_indel_filtPASS_filtGT.vcf \
            2>> "${log_file}"

        echo "Processed INDELs: ${indel_name}"

    else
        echo "File not found: ${indel_file}"
    fi

done
