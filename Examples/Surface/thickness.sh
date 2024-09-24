#!/bin/bash

# directory with script helpers
h="${HOME}/Scripting"

source "${h}/shell/func.sh"

# axes ranges
xr=( 4.2 11 2 )
yr=( 0 0.04 0.01 )
# graph size (x[cm] y[cm] dpi)
size=(10 7.2 600)
# margins - left right top bottom
mar=(6.5 1.6 0.5 2.0)
# output file name
file=${0/.sh/}

f=./thickness.txt

# calculate ranges
source "${h}/gnuplot/ranges.sh"
# define axis label positions (optional values 0-1 for [xl] [yl] [x2l] [y2l])
source "${h}/gnuplot/axis_labels.sh"

gnuplot << EOF
  reset
  set terminal cairolatex standalone lw 2 size ${size[0]}cm,${size[1]}cm \
    header '\\usepackage{amsmath}'

  set output '${file}.tex'
  set multiplot

  set lmargin ${mar[0]}
  set rmargin ${mar[1]}
  set tmargin ${mar[2]}
  set bmargin ${mar[3]}

  set xr  [${xr[0]}  to ${xr[1]}]
  set yr  [${yr[0]}  to ${yr[1]}]

  set xtics  ${xr[2]}  scale 2,1 offset 0.0,0.4
  set ytics  ${yr[2]}  scale 2,1 offset 0.8,0.0

  set mxtics  2
  set mytics  2

  # set style line 100 lt 1 lw 1 lc rgb 'grey' dt (10,10)
  # set style line 200 lt 1 lw 1 lc rgb 'grey' dt (5,5)
  # set grid x y mx my ls 100, ls 200

  # axes labels
  set label 1000 'Thickness' at ${xl}
  set label 1001 'Distribution' at ${yl}

  unset key
  plot '${f}' u 1:2 w lp pt 7 ps 0.3 lw 1 t 'sin(\$\\alpha\$)'
  unset multiplot
  unset output
EOF

source "${h}/gnuplot/latex.sh" "${size[2]}" "${size[0]}"
