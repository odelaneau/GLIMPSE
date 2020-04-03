#!/bin/bash

bin/GLIMPSE_chunk --input NA12878_1x_vcf/NA12878.chr22.1x.vcf.gz --reference reference_panel/1000GP.chr22.noNA12878.sites.vcf.gz --region 22 --window-size 2000000 --buffer-size 200000 --output chunks.chr22.txt
