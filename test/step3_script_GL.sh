#!/bin/bash

mkdir -p NA12878_1x_vcf
bcftools view -m2 -M2 -v snps -G reference_panel/1000GP.chr22.noNA12878.bcf -Oz -o reference_panel/1000GP.chr22.noNA12878.sites.vcf.gz
bcftools index -f reference_panel/1000GP.chr22.noNA12878.sites.vcf.gz 

bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' reference_panel/1000GP.chr22.noNA12878.sites.vcf.gz | bgzip -c > reference_panel/1000GP.chr22.noNA12878.sites.tsv.gz

BAM=NA12878_1x_bam/NA12878.chr22.1x.bam
VCF=reference_panel/1000GP.chr22.noNA12878.sites.vcf.gz
TSV=reference_panel/1000GP.chr22.noNA12878.sites.tsv.gz
REFGEN=reference_genome/human_g1k_v37.chr22.fasta.gz
OUT=NA12878_1x_vcf/NA12878.chr22.1x.vcf.gz

bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} -r 22 ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
bcftools index -f ${OUT}
