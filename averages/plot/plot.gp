
set encoding utf8 
set terminal postscript enhanced 
set terminal postscript eps size 3.5,3 enhanced color \
    font 'Arial,9'  linewidth 1.5
set key box lw 1
set key width 0.5 height 0.5
set key font 'Arial,5'
set key spacing 2
set key right


set terminal postscript eps enhanced size 3.0in,2.2in
set output 'Fig1.eps'

set size 1,1

NZ=2000
NX=2000
SCALE=0.2

set lmargin 7
set rmargin 2

# Axes
set xr [0:NZ] #Time /16
set mytics 5
set mxtics

#----------------
#-  First plot  -
#----------------

# labels and axis
set tmargin at screen 0.88; set bmargin at screen 0.15


set label "{/Symbol m}_x" at 100,10  font 'Arial,20' textcolor rgb "black"
set xrange [0:4500]
set yrange [-0.15:0.15]
set xlabel "t (fms)"
set ylabel "water bisector cosine thetha"


#plotting mgcl2_cosine.dat
plot "../../../../Ions-THz-data-figures/cosine-theta/total.dat" using ($1):($2-0.021) title "AIMD (pure water)"  with line ls 1 lc rgb "black" lw 1, "../../../../Ions-THz-data-figures/cosine-KE-pol/mgcl2_cosine.dat" using ($1):($2-0.023 + 0.05) title "AIMD (mgcl2)"  with line ls 1 lc rgb "red" lw 1, "../water.dat" using ($1):($2 - 0.035) title "AMOEBA pure water"  with line ls 0 lc rgb "black" lw 1, "../mgcl2.dat" using ($1):($2- 0.03) title "AMOEBA mgcl2"  with line ls 0 lc rgb "red" lw 1, "../naf.dat" using ($1):($2+0.016) title "AMOEBA NaF"  with line ls 0 lc rgb "green" lw 1

unset label 
unset arrow 
unset ylabel
unset ytics
unset key

unset output

