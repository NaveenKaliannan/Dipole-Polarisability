
#!/bin/bash

rm -rf build/*
cp CMakeLists.txt build/CMakeLists.txt
mkdir build
cd build
cmake . 
make

./exe ../traj.xyz traj.psf 1 0.4 15.56 389 131 dipol-pol


