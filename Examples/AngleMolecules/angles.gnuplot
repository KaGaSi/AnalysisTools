reset
set terminal pngcairo size 900,700 fontscale 2 linewidth 2 enhanced

set output 'angles.png'
set multiplot

set lmargin 6.0
set rmargin 2.0
set tmargin 0.5
set bmargin 2.0

set xr  [30 to 180]
set yr  [ 0 to   0.1]

set xtics 30    scale 2,1 offset 0.0,0.4
set ytics  0.02 scale 2,1 offset 0.8,0.0

set mxtics  2
set mytics  2

set xlabel 'Angle [degrees]' offset 0,1
set ylabel 'Distribution' offset 2.5,0
set title 'angles in A3BA2 molecule' offset 0,-3

x=0.5
y=0.8
mx=0.4
set key at graph x,y
plot 'angles.txt' u 1:2 w l lw 2 lc rgb '#FF0000' t 'A-A-A angle', \
     'angles.txt' u 1:3 w l lw 2 lc rgb '#0000FF' t 'A-A-B angle', \
     'angles.txt' u 1:4 w l lw 2 lc rgb '#00AA00' t 'A-B-A angle'
set key at graph x+mx,y
plot 'angles.txt' u 1:5 w l lw 2 dt (5,5) lc rgb '#00CCFF' t '1-2-3 angle', \
     'angles.txt' u 1:6 w l lw 2 dt (5,5) lc rgb '#998800' t '2-3-4 angle', \
     'angles.txt' u 1:7 w l lw 2 dt (5,5) lc rgb '#000000' t '3-4-5 angle', \
     'angles.txt' u 1:8 w l lw 2 dt (5,5) lc rgb '#CC00CC' t '4-5-6 angle'
unset multiplot
unset output
