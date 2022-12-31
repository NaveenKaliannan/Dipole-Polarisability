##length of each processor simulation
declare -i len
len=100
##change in field.sh file y+=9

declare -i n
n=0
declare -i j
j=1
for i in {1..50..1}
do
	mkdir sim$i
	cd sim$i
	cp ../field.sh .
	 nohup ./field.sh $i $n $j &
	 j+=len
	cd ..
	n+=1
done

