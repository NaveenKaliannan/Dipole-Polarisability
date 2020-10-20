
#!/bin/bash

##rm -rf build/*
mkdir build
cp CMakeLists.txt build/CMakeLists.txt
cd build
cmake . 
make


taskset --cpu-list 1 ./exe ../../water.arc traj.psf 10000 1 15.6404 384 128 ../field_atomicunit

##./exe /home/naveenk/temp/water.arc traj.psf 65 0.4 145 288000 96000 ../field_atomicunit KEnCos

##./exe ../naf.xyz traj.psf 2500 0.4 16.099 417 143 ../field_atomicunit water$n cosine$n
##./exe ../mgcl2.xyz traj.psf 2500 0.4 16.099 414 148 ../field_atomicunit water$n cosine$n

declare -i n
for i in {1000..183000..400}
do
n+=1 
##echo $n $i
##taskset --cpu-list 1 ./exe ../../water/blyp-d3-$i/OHH.xyz traj.psf 10000 0.4 15.6404 384 128 ../field_atomicunit water$n cosine$n dipol$n
done



##./exe ../../mgcl2-h2o/NVE-pos-1.xyz traj.psf 340 1 16.099 414 148 ../field_atomicunit tinker

##./exe ../../mgcl2-h2o/NVE-pos-1.xyz traj.psf 340 1 16.099 414 148 ../field_atomicunit tinker

##./exe  /home/naveenk/temp_soft/Tinker/exe/OHH.xyz traj.psf 11400 0.4 15.6404 384 128 ../field_atomicunit tinker

#./exe  /home/naveenk/temp_soft/Tinker/exe/OHH.xyz traj.psf 11400 0.4 15.6404 18 8 ../field_atomicunit tinker

##./exe  /home/naveenk/temp_soft/Tinker/exe/OHH.xyz traj.psf 5000 0.4 15.6404 384 128 ../field_atomicunit tinker


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



