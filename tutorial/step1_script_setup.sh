#!/bin/bash

#compile
cd ../chunk; make clean; make -j; chmod u+x bin/GLIMPSE_chunk
cd ../phase; make clean; make -j; chmod u+x bin/GLIMPSE_phase
cd ../ligate; make clean; make -j; chmod u+x bin/GLIMPSE_ligate
cd ../sample; make clean; make -j; chmod u+x bin/GLIMPSE_sample
cd ../concordance; make clean; make -j; chmod u+x bin/GLIMPSE_concordance

#symlinks
cd ../tutorial
rm -rf bin
mkdir -p bin

ln -s ../../chunk/bin/GLIMPSE_chunk bin/GLIMPSE_chunk
ln -s ../../phase/bin/GLIMPSE_phase bin/GLIMPSE_phase
ln -s ../../ligate/bin/GLIMPSE_ligate bin/GLIMPSE_ligate
ln -s ../../sample/bin/GLIMPSE_sample bin/GLIMPSE_sample
ln -s ../../concordance/bin/GLIMPSE_concordance bin/GLIMPSE_concordance

chmod u+x bin/GLIMPSE_chunk
chmod u+x bin/GLIMPSE_phase
chmod u+x bin/GLIMPSE_ligate
chmod u+x bin/GLIMPSE_sample
chmod u+x bin/GLIMPSE_concordance
