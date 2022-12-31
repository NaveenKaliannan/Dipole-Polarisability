ulimit -s unlimited
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64:$LD_LIBRARY_PATH


ion=tmgso4
tinkerexe=/home/naveenk/my_programs/tinkerhp/withoutfield
inputfileloc=/home/naveenk/energytransfer/input/$ion
restartfileloc=/home/naveenk/energytransfer/restarts/$ion
ensemble=nvt


declare -i n

for i in {11..90000..1}
do
n=i
rm water.dyn water.arc  water.xyz
cp $restartfileloc/frame$n.xyz  water.xyz
cp $restartfileloc/frame$n.dyn  water.dyn
cp $inputfileloc/md1.key  water.key
cp $inputfileloc/*.prm  .

taskset --cpu-list 0 $tinkerexe/dynamic.x  water 1000 1 1 2 300

n+=1
tail -407 water.arc |tee  $restartfileloc/frame$n.xyz
tail -1629 water.dyn |tee  $restartfileloc/frame$n.dyn

done




## timestep 1 fms = 1
## dump ps 0.1 = 0.001 means prints every 1 fms 
##               0.01 prints every every 10 fms
##               0.1 prints every 100 fms
##               1 prints every 1000 fms = 1 ps

