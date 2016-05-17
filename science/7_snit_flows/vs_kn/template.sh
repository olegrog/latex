#!/usr/bin/env gnuplot

set term epslatex standalone size 3.2, 2.4 font 8 mono #color dashed #mono
set out "<name>.tex"
set key center bottom maxrows 6 width -10

set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.16
set ylabel '$\displaystyle\int_0^\frac12\!\! <expr> \mathrm{d}x$' \
    offset graph <xcoord>, 0.45 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set nokey

#set xrange [0:0.135]
set xrange [3e-3:0.15]
set log x
set yrange [<ymin>:<ymax>]

heat(x) = <heat>
snit(x) = <snit>
base(x) = x
filter(x,y,xmax) = x > xmax ? 1/0 : y

set macros
dummy = "NaN title ' ' lt -3"

plot \
    filter(x, heat(x), <is_heat>) title "Heat-conduction equation" w l lw 2 lc 4, \
    filter(x, snit(x), 0.05) title "KGF equations$" w l lw 2 lc 3, \
    "asym1.txt" using 1:(base($<column>)) title "KGF equations with d1, k0" w l lw 2 lc 5, \
    "asym2.txt" using 1:(base($<column>)) title "KGF equations with d1, k0, a4, d3" w l lw 2 lc 1, \
    "data.txt" using 1:(base($<column>)) title "Uniform velocity grid" w lp lw 3 lc 2 pt 6, \
    "high.txt" using 1:(base($<column>*<corr>)) title "Nonuniform velocity grid" w lp lw 3 lc 4 pt 4

#    "asym-old.txt" using 1:(base($<column>)) title "KGF equations with $T_{B1}\\ne0$" w l lw 1 lt 1 lc 3, \
