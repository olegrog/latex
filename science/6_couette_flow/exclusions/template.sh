#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set out "<file>.tex"
set key left bottom maxrows 3 width 0

set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.15
set ylabel '$\displaystyle\frac{|\mathcal{M}|}{|\mathcal{N}|}$' offset graph 0.22, 0.45 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy

set xrange [0.007:150]
set yrange [5e-4:2e-1]

set macros
dummy = "NaN title ' ' lt -3"

plot \
    "exclusions-2e4.txt" u 1:<column> title '$20011$'  w lp lw 3 lt 1 lc 1, \
    "exclusions-5e4.txt" u 1:<column> title '$50021$'  w lp lw 3 lt 2 lc 2, \
    "exclusions-1e5.txt" u 1:<column> title '$100003$' w lp lw 3 lt 3 lc 3, \
    "exclusions-2e5.txt" u 1:<column> title '$200003$' w lp lw 3 lt 4 lc 4, \
    "exclusions-5e5.txt" u 1:<column> title '$500009$' w lp lw 3 lt 5 lc 5
