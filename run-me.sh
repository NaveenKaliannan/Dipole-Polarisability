
#!/bin/bash

rm -rf build/*
cp CMakeLists.txt build/CMakeLists.txt
cd build
cmake . 
make
##./exe ../traj.xyz traj.psf 1 0.4 15.6404 3 1
#./exe ../test.xyz traj.psf 1 0.4 15.56 3 1
./exe ../traj.xyz traj.psf 1 0.4 15.56 389 131


