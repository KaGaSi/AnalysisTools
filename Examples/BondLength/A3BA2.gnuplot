reset
set terminal pngcairo size 900,700 fontscale 2 linewidth 2 enhanced

set output 'A3BA2.png'
set multiplot

set lmargin 6.0
set rmargin 2.0
set tmargin 0.5
set bmargin 2.0

set xr  [ 0.2 to 1.1]
set yr  [ 0   to 0.1]

set xtics  0.2  scale 2,1 offset 0.0,0.4
set ytics  0.02 scale 2,1 offset 0.8,0.0

set mxtics 2
set mytics 2

set xlabel 'Distance' offset 0,1
set ylabel 'Distribution' offset 2.5,0
set title 'bond lengths in A3BA2 molecules' offset 0,-3

x=0.99
y=0.8
set key at graph x,y samplen 2 spacing 1.1
plot 'bonds.txt' u 1:2 w l lw 2 lc rgb '#FF0000' t 'A-A bond', \
     'bonds.txt' u 1:3 w l lw 2 lc rgb '#0000FF' t 'A-B bond', \
     'bonds.txt' u 1:4 w l lw 2 dt (5,5) lc rgb '#00AA00' t '1-2 bond', \
     'bonds.txt' u 1:5 w l lw 2 dt (5,5) lc rgb '#000000' t '2-3 bond', \
     'bonds.txt' u 1:6 w l lw 2 dt (5,5) lc rgb '#998800' t '3-4 bond', \
     'bonds.txt' u 1:7 w l lw 2 dt (5,5) lc rgb '#00CCFF' t '4-5 bond', \
     'bonds.txt' u 1:8 w l lw 2 dt (5,5) lc rgb '#BB00BB' t '5-6 bond'
unset multiplot
unset output
