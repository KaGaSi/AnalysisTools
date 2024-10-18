reset
set terminal pngcairo size 900,700 fontscale 2 linewidth 2 enhanced

set output 'end_to_end.png'
set multiplot

set lmargin 6.0
set rmargin 2.0
set tmargin 0.5
set bmargin 2.0

set xr  [ 0.25 to 2.75]
set yr  [ 0    to 0.1]

set xtics 0.5  scale 2,1 offset 0.0,0.4
set ytics 0.02 scale 2,1 offset 0.8,0.0

set mxtics  2
set mytics  2

set xlabel 'Distance' offset 0,1
set ylabel 'Distribution' offset 2.5,0
set title 'distance between the first and last beads' offset 0,-3

x=0.5
y=0.8
set key at graph x,y
plot 'end_to_end.txt' u 1:2 w l lw 2 lc rgb '#FF0000' t 'A3BA2', \
     'end_to_end.txt' u 1:3 w l lw 2 lc rgb '#0000FF' t 'A2B2'
unset multiplot
unset output
