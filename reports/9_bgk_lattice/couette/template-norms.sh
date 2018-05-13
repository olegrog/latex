#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set out "<name>.tex"
set key center bottom

set xlabel '$y$' offset graph 0.48, 0.15
set colors classic

set border 3
set xtics nomirror
set ytics nomirror
set grid

set lmargin 6
set bmargin 2

set xrange [0:.5]
set format y '$10^{%T}$'
set logscale y

U = 0.02
our_width = 2
our_point = 6
exact_width = 1
exact_point = 5
exact_color = 0

x_0 = <buffer>
if (<hybrid>) set arrow from x_0, graph 0 to x_0, graph 1 nohead dt 4 lc 0 lw 3

plot \
    "<name>.txt" using 1:2 title '$L^\infty$' w lp lw our_width pt our_point, \
    "<name>.txt" using 1:3 title '$L^2$' w lp lw our_width pt our_point, \
    "<name>.txt" using 1:4 title '$L^1$' w lp lw our_width pt our_point, \
    "<name>.txt"[7:] using 1:(-$4) title '$-q_x$' w lp lw our_width pt our_point
