##length of each processor simulation
declare -i len
len=100
##change in field.sh file y+=9

declare -i n
n=0
declare -i j
j=1
declare -i pr
pr=0
for i in {35001..35044..1}
do
	mkdir sim$i
	cd sim$i
	cp ../field.sh .
	 nohup taskset --cpu-list $pr  ./field.sh $i $n $j &
	 pr+=1
	 j+=len
	cd ..
	n+=1
done

#$dipolexe/exe /home/naveenk/energytransfer/restarts/$ion/traj$i.arc traj.psf 100 32 15 327 125 /home/naveenk/my_programs/Dipole-Polarisability5/field_atomicunit  $outputfileloc/one$i  $outputfileloc/two$i $outputfileloc/three$i  $outputfileloc/four$i
