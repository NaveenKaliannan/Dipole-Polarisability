
module load mpi
##make

declare -i x=$1


for (( t  = 1; t  <= $x; t +=1 ))
do
./exe $t $2/moments$t.data  $2/polar$t.data $2/charges.out.mulliken
done


