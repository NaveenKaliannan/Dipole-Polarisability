#!/bin/bash

declare -i n
declare -i x=$1
declare -i y=$2

n=0
###for (( i = $x; i <= $y; i+=400 ))
for i in {1400..183000..400}
do
n+=1 
tar -xvf donor/Traj$i.tar.gz
cd acceptor
tar -xvf Traj$i.tar.gz
cd ..
./exe /home/naveenk/temp/traj/AIMD-water-non-eqilibirum/blyp-d3-$i/OHH.xyz traj.psf 10000 0.4 15.6404 384 128 ../field_atomicunit donor/Traj$i-time acceptor/donor/Traj$i-time output/dataset$n
rm donor/Traj$i-time*.dat
rm acceptor/donor/Traj$i-time*.dat
done

