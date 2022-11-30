#Clean up
mkdir -p ../static_bins
mkdir -p resources
rm ../static_bins/GLIMPSE2_*
rm resources/GLIMPSE2_*


#Compile phas
cd ../phase
make clean
make -j static_exe
cp bin/GLIMPSE2_phase_static ../static_bins/.

#Compile ligate
cd ../ligate/
make clean
make -j static_exe
cp bin/GLIMPSE2_ligate_static ../static_bins/.

#Compile concordance
cd ../concordance/
make clean
make -j static_exe
cp bin/GLIMPSE2_concordance_static ../static_bins/.

#Compile chunk
cd ../chunk/
make clean
make -j static_exe
cp bin/GLIMPSE2_chunk_static ../static_bins/.

#Compile split_reference
cd ../split_reference/
make clean
make -j static_exe
cp bin/GLIMPSE2_split_reference_static ../static_bins/.

#Buld docker image
LAB=glimpse2_$(git log -1 --format=%cd --date=short)\_$(git rev-parse --short HEAD)

cd ../docker/
mkdir -p resources
cp ../static_bins/GLIMPSE2* resources/.

docker build -t $LAB -f Dockerfile .
docker save $LAB | gzip -c > $LAB\.tar.gz
