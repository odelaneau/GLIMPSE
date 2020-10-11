#!/bin/bash

mkdir -p GLIMPSE_sample
VCF=GLIMPSE_ligate/NA12878.chr22.merged.bcf

OUT=GLIMPSE_sample/NA12878.chr22.sampled.bcf

../sample/bin/GLIMPSE_sample --input ${VCF} --solve --output ${OUT}
bcftools index -f ${OUT}
