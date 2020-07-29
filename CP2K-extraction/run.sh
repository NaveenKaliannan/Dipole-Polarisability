
module load mpi
##make


for t in {1..200..1}
do
./exe $t /home/naveenk/temp/Pol/Frame$t/dipole.out-moments-1.dat /home/naveenk/temp/Pol/Frame$t/polar.out-raman-1.data /home/naveenk/temp/Pol/Frame$t/charges.out.mulliken
done


