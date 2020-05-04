#!/bin/bash

bin/GLIMPSE_chunk --input reference_panel/1000GP.chr22.noNA12878.sites.vcf.gz --region 22 --window-size 2000000 --buffer-size 200000 --output chunks.chr22.txt
