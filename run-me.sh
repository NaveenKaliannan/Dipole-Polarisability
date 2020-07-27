
#!/bin/bash

rm -rf build/*
mkdir build
cp CMakeLists.txt build/CMakeLists.txt
cd build
cmake . 
make

./exe ../test.xyz traj.psf 1 0.4 15.56 4 2 ../field_atomicunit  dipol-pol


declare -i n
n=0
for i in {1000..183000..400}
do
n+=1
##echo $n
##./exe /home/naveenk/temp/water/blyp-d3-$i/OHH.xyz traj.psf 11400 0.4 15.6404 384 128 ../field_atomicunit dipol-pol$n
done

cd ..


##./Average 456 9 456 build/dipol-pol mean.dat; xmgrace -block mean.dat -bxy 1:9 field_atomicunit_timestep

