
#!/bin/bash

rm -rf build/*
cp CMakeLists.txt build/CMakeLists.txt
cd build
cmake . 
make
./exe /home/naveenk/temp/NVE-pos-1.xyz traj.psf 1000 0.4 15.7202 405 135


