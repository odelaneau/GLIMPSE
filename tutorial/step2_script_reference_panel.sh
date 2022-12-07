#!/bin/bash

mkdir -p reference_panel
rm -f reference_panel/*

wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz{,.tbi}

echo -e "Checking MD5 sums of downloaded files..."
vcf_md5_to_test=aaf19d9c7ffcd86b34275899ddc898e7
vcf_md5_from_file=$(md5sum CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz | cut -d " " -f1)

if [[ $vcf_md5_to_test == $vcf_md5_from_file ]]; then
	echo -e "\n\e[92mSUCCESS\e[39m\nMD5 sum of downloaded vcf.gz file matches the expected hash!"
else
	echo -e "\n\e[91mFAILURE\e[39m\nMD5 sum of downloaded vcf.gz file does not match the expected hash. Please try downloading the file again!"
	exit 1
fi

tbi_md5_to_test=b7da8d96ae7e24f75133c41e9854a895
tbi_md5_from_file=$(md5sum CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz.tbi | cut -d " " -f1)

if [[ $tbi_md5_to_test == $tbi_md5_from_file ]]; then
	echo -e "\n\e[92mSUCCESS\e[39m\nMD5 sum of downloaded vcf.gz.tbi file matches the expected hash!"
else
	echo -e "\n\e[91mFAILURE\e[39m\nMD5 sum of downloaded vcf.gz.tbi file does not match the expected hash. Please try downloading the file again!"
	exit 1
fi

# Multiple processing commands piped together
CHR=22

bcftools norm -m -any CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps -s ^NA12878,NA12891,NA12892 --threads 4 -Ob -o reference_panel/1000GP.chr22.noNA12878.bcf
bcftools index -f reference_panel/1000GP.chr22.noNA12878.bcf --threads 4

#rm -f CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz*

