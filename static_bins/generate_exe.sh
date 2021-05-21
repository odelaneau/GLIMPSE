#!/bin/bash

rm GLIMPSE*

for TOOL in chunk concordance ligate phase sample snparray stats; do
	cd ../${TOOL}
	make clean;
	make -j static_exe;
	cp bin/GLIMPSE*_static ../static_bins/
done

cd ../static_bins
