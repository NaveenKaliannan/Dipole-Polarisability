
#!/bin/bash

rm -rf build/*
mkdir build
cp CMakeLists.txt build/CMakeLists.txt
cd build
cmake . 
make

./exe ../mgcl2/NVE-pos-1.xyz traj.psf 10000 1 15.56 3 3 ../field_atomicunit 


##./exe ../nacl/NVE-pos-1.xyz traj.psf 7000 1 15.56 2 2 ../field_atomicunit 

##./exe ../traj.xyz traj.psf 1 1 15.56 389 131 ../field_atomicunit 



##./exe ../traj.xyz traj.psf 9 1 15.56 2 2 ../field_atomicunit 


##./exe ../test.xyz traj.psf 1000 1 15.56 3 3 ../field_atomicunit 

##./exe ../test3.xyz traj.psf 1000 1 15.56 15 5 ../field_atomicunit  

##./exe ../test2.xyz traj.psf 1000 1 15.56 18 8 ../field_atomicunit 

##./exe ../test.xyz traj.psf 163 1 15.56 3 3 ../field_atomicunit 

##./exe ../../nacl.xyz traj.psf 200 1 15.56 2 2 ../field_atomicunit 

##./exe ../mgcl2-nacl.xyz traj.psf 5000 1 15.56 7 7 ../field_atomicunit 



## argumets
## 1. trajectory
## 2. PSF file (not mandatory)
## 3. No of Frames
## 4. time step
## 5. Box Length
## 6. No of atoms
## 7. Number of Molecules
## 8. External field (in atomic unit)

## Important Note: Partial charges, atomic polarisability for ions should be 
## described in parameters() function in  dipol.cpp file

## SO4 coordinates in xyz file should in the following order O S O O O
## H2O coordinates in xyz file should in the following order O H H

