
#!/bin/bash

rm -rf build/*
mkdir build
cp CMakeLists.txt build/CMakeLists.txt
cd build
cmake . 
make

##./exe ../test.xyz traj.psf 1 0.4 15.56 6 2 dipol-pol

./exe /home/naveenk/temp/dippol/pos.xyz traj.psf 10 0.4 45 7200 2400 dipol-pol field


