#!/usr/bin/env gnuplot

set term epslatex standalone size 3.2, 2.4 font 9 color
set out "<name>.tex"
set key center bottom maxrows 6 width -10

set xlabel '$\mathrm{Kn}$' offset graph 0.48, 0.18
set ylabel '$\displaystyle\int_0^\frac12\!\! <expr> \mathrm{d}x$' \
    offset graph <xcoord>, 0.45 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 5
set rmargin 1
set bmargin 2
set nokey
set bars 2

#set xrange [0:0.135]
set xrange [3e-3:0.15]
set log x
set yrange [<ymin>:<ymax>]

heat(x) = <heat>
snit(x) = <snit>
filter(x,y,xmax) = x > xmax ? 1/0 : y
errT(x,y,xmax) = x > xmax ? 1/0 : y/1.004
errUk(x,y,xmax) = x > xmax ? 1/0 : 3e-4*y/x

set macros
dummy = "NaN title ' ' lt -3"

plot \
    filter(x, heat(x), <is_heat>) title "Heat-conduction equation" w l lw 3 lc 6 dt 3, \
    filter(x, snit(x), 0.05) title "KGF equations$" w l lw 3 lc 1 dt 2, \
    "asym1.txt" using 1:<column> title "KGF equations with d1, k0" w l lw 3 lc 2 dt 4, \
    "asym2.txt" using 1:<column> title "KGF equations with d1, k0, a4, d3" w l lw 3 lc 4 dt 1, \
    "data.txt"  using 1:<column> title "Uniform velocity grid" w p lw 3 lc 7 dt 1 pt 6, \
    "high.txt"  using 1:<column> title "Hermite velocity grid" w p lw 3 lc 7 dt 1 pt 4, \
    "data.txt"  using 1:<column>:(errUk($1,$<column>, 1-<is_temp>)) notitle w yerrorbars lw 1 lc 8 dt 1 pt -1

#    "high2.txt" using 1:<column> title "Squared velocity grid" w p lw 3 lc 7 dt 1 pt 8 , \
#    "high2.txt" using 1:<column>:(errT($1,$<column>, <is_temp>)):<column> notitle w yerrorbars lw 1 lc 8 dt 1 pt -1, \


#    "data.txt"  using 1:($<column>/1.04) title "Uniform velocity grid" w lp lw 3 lc 2 dt 2 pt 7, \
#    filter(x, snit(x)*(1-5*x), 0.05) title "KGF equations$" w l lw 1 lc 1 dt 1, \
#    filter(x, snit(x)*(1-20*x**2), 0.07) title "KGF equations$" w l lw 1 lc 1 dt 1, \
#    filter(x, snit(x)*(1-5*x**3), 0.1) title "KGF equations$" w l lw 1 lc 1 dt 1, \
