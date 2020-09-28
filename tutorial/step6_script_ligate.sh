#!/bin/bash

mkdir -p GLIMPSE_ligate
LST=GLIMPSE_ligate/list.chr22.txt
ls GLIMPSE_impute/NA12878.chr22.imputed.*.bcf > ${LST}

OUT=GLIMPSE_ligate/NA12878.chr22.merged.bcf

../ligate/bin/GLIMPSE_ligate --input ${LST} --output $OUT
bcftools index -f ${OUT}
