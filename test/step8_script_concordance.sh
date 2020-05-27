#!/bin/bash

mkdir -p GLIMPSE_concordance

./bin/GLIMPSE_concordance --input concordance.lst --minDP 8 --output GLIMPSE_concordance/output --minPROB 0.9999 --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.02000 0.05000 0.10000 0.15000 0.20000 0.25000 0.30000 0.35000 0.40000 0.45000 0.50000 --info_af AF_nfe --thread 4

cd plot;
./concordance_plot.py


