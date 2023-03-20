#!/bin/bash

rm -rf build
mkdir build
cp CMakeLists.txt build/CMakeLists.txt
cd build
cmake . 
make

declare -i n
declare -i x=$1
declare -i y=$2

./exe /home/naveenk/temp/traj/AIMD-Equilibrium/purewater/blyp-d3-600/OHH.xyz  traj.psf 50000 0.4 15.6404 384 128 ../field_atomicunit 
#./exe /home/naveenk/temp/traj/200ps-AMOEBA-AMBER-SPC-trajectories-Equilibirum/200ps-PFFMD-FFMD-trajectories-equibirum/gromacs/water.xyz  traj.psf 10000 1 15.6404 384 128 ../field_atomicunit 
#./exe /home/naveenk/temp/traj/200ps-AMOEBA-AMBER-SPC-trajectories-Equilibirum/200ps-PFFMD-FFMD-trajectories-equibirum/tinker/water.xyz  traj.psf 10000 1 15.6404 384 128 ../field_atomicunit 

#./exe ../$0 traj.psf 50000 1 15.6404 384 128 ../field_atomicunit perm_polaniso  permplusDID_polaniso permplusXDID_polaniso permplusXXDID_polaniso permanet Induced total
#./exe ../$1 traj.psf 50000 1 15.5732 381 131 ../field_atomicunit perm_polaniso  permplusDID_polaniso permplusXDID_polaniso permplusXXDID_polaniso permanet Induced total
#./exe ../$1 traj.psf 50000 1 15.6287 378 136 ../field_atomicunit perm_polaniso  permplusDID_polaniso permplusXDID_polaniso permplusXXDID_polaniso permanet Induced total
#./exe ../$1 traj.psf 50000 1 15 327 125 ../field_atomicunit perm_polaniso  permplusDID_polaniso permplusXDID_polaniso permplusXXDID_polaniso permanet Induced total
#./exe ../$1 traj.psf 100000 1 15.6404 384 128 ../field_atomicunit translation_function
#./exe ../$1 traj.psf 100000 1 15 327 125 ../field_atomicunit translation_function
#./exe ../$1 traj.psf 120000 1 15 327 125 ../field_atomicunit 
#./exe ../$1 traj.psf 100000 1 15.6404 384 128 ../field_atomicunit FFMD_water
#./exe ../$1 traj.psf 100000 1 15 327 125 ../field_atomicunit FFMD_mgcl2_4mol
#./exe ../$1 traj.psf 100000 1 15.6287 378 136 ../field_atomicunit FFMD_mgcl2_2mol
#./exe ../$1 traj.psf 100000 1 15.5732 381 131 ../field_atomicunit FFMD_mgcl2_1mol  


## argumets
## 1. trajectory
## 2. PSF file (not mandatory)
## 3. No of Frames
## 4. time step
## 5. Box Length
## 6. No of atoms
## 7. Number of Molecules
## 8. External field (in atomic unit)

## Important Note: Partial charges, atomic and molecular polarisability for ions should be 
## described in parameters() function in  dipol.cpp file
## SO4 coordinates in xyz file should in the following order O S O O O
## H2O coordinates in xyz file should in the following order O H H
##check io.h file for deltat, dipol.cpp for ncell for replica
