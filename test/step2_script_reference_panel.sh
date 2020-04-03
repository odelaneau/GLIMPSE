#!/bin/bash

mkdir -p reference_panel
rm -f reference_panel/*

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz{,.tbi} 

bcftools view -s ^NA12878 ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -Ob -o reference_panel/1000GP.chr22.noNA12878.bcf
bcftools index -f reference_panel/1000GP.chr22.noNA12878.bcf

rm ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz*
