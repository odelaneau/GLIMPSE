#!/bin/bash

mkdir -p NA12878_1x_vcf
bcftools view -G -m 2 -M 2 -v snps reference_panel/1000GP.chr22.noNA12878.bcf -Oz -o reference_panel/1000GP.chr22.noNA12878.sites.vcf.gz
bcftools index -f reference_panel/1000GP.chr22.noNA12878.sites.vcf.gz 

bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' reference_panel/1000GP.chr22.noNA12878.sites.vcf.gz | bgzip -c > reference_panel/1000GP.chr22.noNA12878.sites.tsv.gz

##Bugfix for [E::_regions_match_alleles] Compressed and indexed targets file is required
#suggested by @alek0991 https://github.com/odelaneau/GLIMPSE/issues/3
#requres htslib in your path

tabix -s1 -b2 -e2 reference_panel/1000GP.chr22.noNA12878.sites.tsv.gz

BAM=NA12878_1x_bam/NA12878.bam
VCF=reference_panel/1000GP.chr22.noNA12878.sites.vcf.gz
TSV=reference_panel/1000GP.chr22.noNA12878.sites.tsv.gz
REFGEN=reference_genome/hs38DH.chr22.fa.gz
OUT=NA12878_1x_vcf/NA12878.chr22.1x.vcf.gz

bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} -r chr22 ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
bcftools index -f ${OUT}
