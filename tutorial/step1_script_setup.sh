#!/bin/bash

#compile
cd ../chunk; make clean; make -j desktop; chmod u+x bin/GLIMPSE2_chunk
cd ../split_reference; make clean; make -j desktop; chmod u+x bin/GLIMPSE2_split_reference
cd ../phase; make clean; make -j  desktop; chmod u+x bin/GLIMPSE2_phase
cd ../ligate; make clean; make -j desktop; chmod u+x bin/GLIMPSE2_ligate
cd ../concordance; make clean; make -j desktop; chmod u+x bin/GLIMPSE2_concordance

#symlinks
cd ../tutorial
rm -rf bin
mkdir -p bin

ln -s ../../chunk/bin/GLIMPSE2_chunk bin/GLIMPSE2_chunk
ln -s ../../split_reference/bin/GLIMPSE2_split_reference bin/GLIMPSE2_split_reference
ln -s ../../phase/bin/GLIMPSE2_phase bin/GLIMPSE2_phase
ln -s ../../ligate/bin/GLIMPSE2_ligate bin/GLIMPSE2_ligate
ln -s ../../concordance/bin/GLIMPSE2_concordance bin/GLIMPSE2_concordance

chmod u+x bin/GLIMPSE2_chunk
chmod u+x bin/GLIMPSE2_split_reference
chmod u+x bin/GLIMPSE2_phase
chmod u+x bin/GLIMPSE2_ligate
chmod u+x bin/GLIMPSE2_concordance
