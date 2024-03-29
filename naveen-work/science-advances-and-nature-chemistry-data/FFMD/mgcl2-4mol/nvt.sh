
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/naveenk/my_programs/gcc-9/bin/lib64:$LD_LIBRARY_PATH
export LIBRARY_PATH=/home/naveenk/my_programs/gcc-9/bin/lib64:$LIBRARY_PATH
export CC=/home/naveenk/my_programs/gcc-9/bin/bin/gcc
export CXX=/home/naveenk/my_programs/gcc-9/bin/bin/g++


ion=4mol
gmxexe=/home/naveenk/my_programs/gromacs-2020.4/bin/bin
inputfileloc=/home/naveenk/energytransfer/input/$ion
ensemblefileloc=/home/naveenk/energytransfer/input/mdp
restartfileloc=/home/naveenk/energytransfer/restarts/$ion
ensemble=nve

declare -i n

for i in {50000..500000..1}
do
n=i
$gmxexe/gmx_d grompp -f $ensemblefileloc/$ensemble.mdp -c $restartfileloc/frame$i.gro -t $restartfileloc/frame$i.cpt -p $inputfileloc/topol.top -o $ensemble.tpr  -maxwarn 5
$gmxexe/gmx_d mdrun -v -deffnm $ensemble  -nt 1
##$gmxexe/gmx_d  make_ndx -f $ensemble.gro -o $ensemble.ndx
##$gmxexe/gmx_d energy -f $ensemble.edr -o $ensemble.xvg
n+=1
mv nve.gro $restartfileloc/frame$n.gro 
mv nve.cpt $restartfileloc/frame$n.cpt 
rm *trr* *gro* *cpt*
rm *log* *xtc*
rm *mdout* *edr* 
rm *tpr*
done


##gmx trjconv -f $ensemble.trr -s $ensemble.tpr  -o $ensemble.gro  -ndec 5 
##vmd -e PSF.tcl

#a=$ensemble.xyz

#sed -i 's/OW/O/g' $a
#sed -i 's/HW1/H/g' $a
#sed -i 's/HW2/H/g' $a
#sed -i '/MW /d' $a
#sed -i ':a;N;$!ba;s/547\n generated by VMD/414\n generated by VMD/g' $a


##rm -rf \#*
#rm -rf *.mdp *.top *.grp



