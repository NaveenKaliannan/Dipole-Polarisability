
module load mpi
##make

declare -i x=$1


for (( t  = 1; t  <= $x; t +=1 ))
do
##./exe $t $2/moments$t.data  $2/polar$t.data $2/charges.out.mulliken
./exe $t $2/dipole.out-moments-1_$t.dat  $2/polar.out-raman-1_$t.data $2/charges.out-1.mulliken
done


