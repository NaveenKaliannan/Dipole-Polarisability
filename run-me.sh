
#!/bin/bash

rm -rf build/*
cp CMakeLists.txt build/CMakeLists.txt
mkdir build
cd build
cmake . 
make

declare -i n
for i in {1000..183000..400}
do
n+=1
echo $n
./exe /home/naveenk/temp/water/blyp-d3-$i/OHH.xyz traj.psf 11400 0.4 15.6404 384 128 dipol-pol$n
done

