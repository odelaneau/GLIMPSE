#!/bin/bash

mkdir -p GLIMPSE_ligate
LST=GLIMPSE_ligate/list.chr22.txt
ls -1v GLIMPSE_impute/NA12878_imputed_*.bcf > ${LST}

OUT=GLIMPSE_ligate/NA12878_chr22_ligated.bcf

./bin/GLIMPSE2_ligate --input ${LST} --output $OUT
