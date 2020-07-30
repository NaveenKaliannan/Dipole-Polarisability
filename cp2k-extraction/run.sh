
module load mpi
##make


for t in {1..1000..1}
do
./exe $t /home/naveenk/temp/1-mgcl2/Frame$t/dipole.out-moments-1.dat /home/naveenk/temp/1-mgcl2/Frame$t/polar.out-raman-1.data /home/naveenk/temp/1-mgcl2/Frame$t/charges.out.mulliken
done


