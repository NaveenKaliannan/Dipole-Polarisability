a=2500 ; b=7 ; c=10000  ; ./run-me.sh ../build/watercosine $a $b $c ; mv mean.dat watercosine.dat ; ./run-me.sh ../build/mgcl2cosine $a $b $c ; mv mean.dat mgcl2cosine.dat ; ./run-me.sh ../build/nafcosine $a $b $c ; mv mean.dat nafcosine.dat ; xmgrace -block watercosine.dat -bxy 1:2 -block  mgcl2cosine.dat -bxy 1:2 -block nafcosine.dat -bxy 1:2  ../pulse/100-times-stronger-x.data

xmgrace -block watercosine.dat -bxy 1:3 -block  mgcl2cosine.dat -bxy 1:3 -block nafcosine.dat -bxy 1:3  ../pulse/100-times-stronger-x.data

 b=16 ; ./run-me.sh ../build/water $a $b $c ; mv mean.dat water.dat ; ./run-me.sh ../build/mgcl2 $a $b $c ; mv mean.dat mgcl2.dat ; ./run-me.sh ../build/naf $a $b $c ; mv mean.dat naf.dat ; xmgrace -block water.dat -bxy 1:4 -block  mgcl2.dat -bxy 1:4 -block naf.dat -bxy 1:4  ../pulse/100-times-stronger-x.data

xmgrace -block water.dat -bxy 1:3 -block  mgcl2.dat -bxy 1:3 -block naf.dat -bxy 1:3  ../pulse/100-times-stronger-x.data


##rm *.dat
