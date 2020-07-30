
module load mpi
##make


for t in {1..1000..1}
do
./exe $t /home/naveenk/temp/5-H2O/Frame$t/dipole.out-moments-1.dat /home/naveenk/temp/5-H2O/Frame$t/polar.out-raman-1.data /home/naveenk/temp/5-H2O/Frame$t/charges.out.mulliken
done


