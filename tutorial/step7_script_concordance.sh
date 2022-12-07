#!/bin/bash

mkdir -p GLIMPSE_concordance

./bin/GLIMPSE2_concordance --input concordance.lst --min-val-dp 8 --output GLIMPSE_concordance/output --min-val-gl 0.9999 --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000 --af-tag AF_nfe --thread 4

cd plot;
./concordance_plot.py


