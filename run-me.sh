
#!/bin/bash

rm -rf build/*
mkdir build
cp CMakeLists.txt build/CMakeLists.txt
cd build
cmake . 
make

taskset --cpu-list 1 ./exe ../../mgcl2-ho/NVE-pos-1.xyz traj.psf 10000 1 15.56 18 8 ../field_atomicunit 




declare -i n
##for i in {1..1..1}
for i in {1000..1000..400}
do
n+=1 
##echo $n $i
##taskset --cpu-list 1 ./exe ../../water/blyp-d3-$i/OHH.xyz traj.psf 400 0.4 15.6404 384 128 ../field_atomicunit water$n cosine$n
##./exe ../../test/water$i.xyz traj.psf 11400 0.4 15.6404 384 128 ../field_atomicunit water$n cosine$n
##./exe ../../test/water$i.xyz traj.psf 11400 0.4 16.099 414 148 ../field_atomicunit water$n cosine$n
done


##//./exe ../../test/NVE-pos-1.xyz traj.psf 1600 1 16.096 417 143 ../field_atomicunit tinker


##./exe ../../mgcl2-h2o/NVE-pos-1.xyz traj.psf 340 1 16.099 414 148 ../field_atomicunit tinker

##./exe ../../mgcl2-h2o/NVE-pos-1.xyz traj.psf 340 1 16.099 414 148 ../field_atomicunit tinker

##./exe  /home/naveenk/temp_soft/Tinker/exe/OHH.xyz traj.psf 11400 0.4 15.6404 384 128 ../field_atomicunit tinker

#./exe  /home/naveenk/temp_soft/Tinker/exe/OHH.xyz traj.psf 11400 0.4 15.6404 18 8 ../field_atomicunit tinker

##./exe  /home/naveenk/temp_soft/Tinker/exe/OHH.xyz traj.psf 5000 0.4 15.6404 384 128 ../field_atomicunit tinker




declare -i n
##for i in {1..2..1}
for i in {1000..183000..400}
do
n+=1 
##echo $n
##./exe ../../water/blyp-d3-$i/OHH.xyz traj.psf 11400 0.4 15.6404 384 128 ../field_atomicunit water$n cosine$n
##./exe ../../test/water$i.xyz traj.psf 11400 0.4 15.6404 384 128 ../field_atomicunit water$n cosine$n
done


##./exe ../../mgcl/OHH.xyz traj.psf 10 1 15.6404 2 2 ../field_atomicunit
  

##./exe ../../mgcl2-h2o/NVE-pos-1.xyz traj.psf 340 1 16.099 414 148 ../field_atomicunit

##./exe ../../../Tinker/exe/OHH.xyz traj.psf 11400 0.4 15.6404 384 128 ../field_atomicunit
  
##./exe ../../frame.xyz traj.psf 1 1 15.6404 384 128 ../field_atomicunit

##./exe ../../mgcl2-h2o/NVE-pos-1.xyz traj.psf 2100 1 16.099 414 148 ../field_atomicunit 

##./exe ../../mgcl2-h2o/NVE-pos-1.xyz traj.psf 10000 1 15.56 18 8 ../field_atomicunit 

##./exe ../../mgcl2/NVE-pos-1.xyz traj.psf 2420 1 15.56 3 3 ../field_atomicunit 

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

### assigning executables for each cores in cpu
### for i in {0..39}; do echo "n+=1" ;  echo "taskset --cpu-list $i ./exe  ../../sim\$n/water.arc traj.psf 11400 0.4 15.6404 384 128 ../field_atomicunit water\$n watercosine\$n &"   ; done



