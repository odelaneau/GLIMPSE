#!/bin/bash

mkdir -p reference_panel/split

#VCF=NA12878_1x_vcf/NA12878.chr22.1x.vcf.gz
REF=reference_panel/1000GP.chr22.noNA12878.bcf
MAP=../maps/genetic_maps.b38/chr22.b38.gmap.gz
BAM=NA12878_1x_bam/NA12878.bam

while IFS="" read -r LINE || [ -n "$LINE" ]; 
do   
	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
	IRG=$(echo $LINE | cut -d" " -f3)
	ORG=$(echo $LINE | cut -d" " -f4)
	#OUT=GLIMPSE_impute/NA12878.chr22.imputed.${ID}.bcf

	./bin/GLIMPSE2_split_reference --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output reference_panel/split/1000GP.chr22.noNA12878
done < chunks.chr22.txt
