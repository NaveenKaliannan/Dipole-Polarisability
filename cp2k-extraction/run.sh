
module load mpi
##make


for t in {1..1000..1}
do
./exe $t /home/naveenk/temp/mgcl2-h2o/Frame$t/dipole.out-moments-1.dat /home/naveenk/temp/mgcl2-h2o/Frame$t/polar.out-raman-1.data /home/naveenk/temp/mgcl2-h2o/Frame$t/charges.out.mulliken
done


