
a=19200 ; b=15 ; c=500  ; ./run-me.sh ../build2/waterdipol $a $b $c ; mv mean.dat waterdipol.dat ; ./run-me.sh ../build2/mgcl2dipol $a $b $c ; mv mean.dat mgcl2dipol.dat ; ./run-me.sh ../build2/nafdipol $a $b $c ; mv mean.dat nafdipol.dat ; 

xmgrace -block waterdipol.dat -bxy 1:5 -block  mgcl2dipol.dat -bxy 1:5 -block nafdipol.dat -bxy 1:5  ../pulse/100-times-stronger-x.data

xmgrace -block waterdipol.dat -bxy 1:9 -block  mgcl2dipol.dat -bxy 1:9 -block nafdipol.dat -bxy 1:9  ../pulse/100-times-stronger-x.data

a=19200 ; b=9 ; c=380  ; ./run-me.sh ../build2/watercosine $a $b $c ; mv mean.dat watercosine.dat ; ./run-me.sh ../build2/mgcl2cosine $a $b $c ; mv mean.dat mgcl2cosine.dat ; ./run-me.sh ../build2/nafcosine $a $b $c ; mv mean.dat nafcosine.dat ; xmgrace -block watercosine.dat -bxy 1:2 -block  mgcl2cosine.dat -bxy 1:2 -block nafcosine.dat -bxy 1:2  ../pulse/100-times-stronger-x.data

xmgrace -block watercosine.dat -bxy 1:3 -block  mgcl2cosine.dat -bxy 1:3 -block nafcosine.dat -bxy 1:3  ../pulse/100-times-stronger-x.data

 b=22 ; ./run-me.sh ../build2/water $a $b $c ; mv mean.dat water.dat ; ./run-me.sh ../build2/mgcl2 $a $b $c ; mv mean.dat mgcl2.dat ; ./run-me.sh ../build2/naf $a $b $c ; mv mean.dat naf.dat ; xmgrace -block water.dat -bxy 1:4 -block  mgcl2.dat -bxy 1:4 -block naf.dat -bxy 1:4  ../pulse/100-times-stronger-x.data

xmgrace -block water.dat -bxy 1:3 -block  mgcl2.dat -bxy 1:3 -block naf.dat -bxy 1:3  ../pulse/100-times-stronger-x.data


xmgrace -block water.dat -bxy 1:3 -block  mgcl2.dat -bxy 1:6 -block naf.dat -bxy 1:6 


##rm *.dat
