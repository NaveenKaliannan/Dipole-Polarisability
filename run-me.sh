
#!/bin/bash

rm -rf build/*
cp CMakeLists.txt build/CMakeLists.txt
mkdir build
cd build
cmake . 
make

./exe ../test.xyz traj.psf 1 0.4 15.56 6 2 dipol-pol


