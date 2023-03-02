#!/bin/bash
MAP=../maps/genetic_maps.b38/chr22.b38.gmap.gz
bin/GLIMPSE2_chunk --input reference_panel/1000GP.chr22.noNA12878.sites.vcf.gz --region chr22 --sequential --output chunks.chr22.txt --map ${MAP}
