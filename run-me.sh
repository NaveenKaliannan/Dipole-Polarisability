
#!/bin/bash



##rm -rf build/*
mkdir build
cp CMakeLists.txt build/CMakeLists.txt
cd build
cmake . 
make

./exe ../$1 traj.psf 10000 0.4 7 18 18 ../field_atomicunit perm_polaniso  permplusDID_polaniso permplusXDID_polaniso permplusXXDID_polaniso permanet Induced total

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


#module load chem/CP2K/5.1-foss-2018b
#module load vis/GLib/2.64.1-GCCcore-9.3.0

#./exe ../$1 traj.psf 25000 0.4 15.56 389 131 ../field_atomicunit one two three four
##./exe ../$1 traj.psf 25000 0.4 16.099 414 148 ../field_atomicunit one two three four
##./exe ../$1 traj.psf 80000 0.4 15.6404 384 128 ../field_atomicunit one two three four



##for loop for many files
#declare -i n
#n=0
#for i in {1000..183000..400}
#do
#n+=1 
##echo $n $i
##./exe  ../../traj/water/blyp-d3-$i/OHH.xyz traj.psf 7000 1 15.6404 384 128 ../field_atomicunit  TBT$n.dat TBP$n.dat TBI$n.dat
##taskset --cpu-list 1 ./exe ../../water/blyp-d3-$i/OHH.xyz traj.psf 10000 0.4 15.6404 384 128 ../field_atomicunit water$n cosine$n dipol$n
#done

#for i in {1..8000..1}
#do
#a=$(wc -l sim$i/water.arc )
#b=$(echo $a | awk '{x=$1/385; print x}')
#if [ $b -lt 10000 ] ;  then
# echo $i  $b
#fi
#done
### assigning executables for each cores in cpu
### for i in {0..39}; do echo "n+=1" ;  echo "taskset --cpu-list $i ./exe  /scratch/hpc-prf-wcat/naveenk/withfield/water/sim\$n/water.arc traj.psf 10000 0.4 15.6404 384 128 ../field_atomicunit data/watercos\$n data/waterdipol\$n &"   ; done
### echo "declare -i n=\$1" ; for i in {0..39}; do echo "cd sim\$n" ;  echo " mpirun -np 1 taskset --cpu-list $i /scratch/hpc-prf-wcat/naveenk/withfield/bin/dynamic water 10000 0.4 0.0004 1 &" ; echo "cd .." ; echo "n+=1" ;  done
##declare -i n ; n=20000; for i in {10001..11360..1} ; do n+=1 ;  echo $i $n ; done
##for i in {0..39}; do echo "n+=1" ;  echo "taskset --cpu-list $i ./exe  /scratch/hpc-prf-wcat/naveenk/withfield/naf/sim\$n/water.arc traj.psf 10000 0.4 16.099 417 143 ../field_atomicunit data/nafcos\$n data/nafdipol\$n &"   ; done
##for i in {0..39}; do echo "n+=1" ;  echo "taskset --cpu-list $i ./exe  /scratch/hpc-prf-wcat/naveenk/withfield/mgcl2/sim\$n/water.arc traj.psf 10000 0.4 16.099 414 148 ../field_atomicunit data/mgcl2cos\$n data/mgcl2dipol\$n &"   ; done


