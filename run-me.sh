
#!/bin/bash

rm -rf build/*
cp CMakeLists.txt build/CMakeLists.txt
cd build
cmake . 
make
./exe ../../NVE-pos-1.xyz traj.psf 5000 0.5 15.72 405 135


