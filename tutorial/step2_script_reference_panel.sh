#!/bin/bash

mkdir -p reference_panel
rm -f reference_panel/*

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz{,.tbi}

# Generate a chromosome renaming file
for CHR in {1..23} X ; do 
    echo ${CHR} chr${CHR}
done >> chr_names.txt

# Multiple processing commands piped together
CHR=22
bcftools annotate --rename-chrs chr_names.txt \
    ALL.chr${CHR}_GRCh38.genotypes.20170504.vcf.gz -Ou | \
bcftools view -m 2 -M 2 -v snps -s ^NA12878 --threads 4 -Ob -o reference_panel/1000GP.chr22.noNA12878.bcf
bcftools index -f reference_panel/1000GP.chr22.noNA12878.bcf

rm ALL.chr22_GRCh38.genotypes.20170504.vcf.gz*
rm chr_names.txt
