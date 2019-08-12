#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set out "<name>.tex"
set key center bottom

set xlabel '$N_0$' offset graph 0.47, 0.15
set ylabel '$\displaystyle\frac{t_\mathrm{DV}}{t_\mathrm{hyb}}$' offset graph 0.13, 0.44 rotate by 0
set colors classic

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2

set xrange [2:1e3]
set yrange [1:14]
set logscale xy
set grid

our_width = 2
our_point = 6
k = 1 # 1.0479

plot \
    "hilbert.txt" using 1:($2/$3) title 'Linux, Intel Xeon E5-2620v4' w lp lw our_width pt our_point, \
    "mvs.txt" using 1:($2/$3) title 'Linux, Intel Xeon E5-2690' w lp lw our_width pt our_point, \
    "macrog.txt" using 1:($2/$3) title 'Darwin, Intel Core i7-7820HQ' w lp lw our_width pt our_point
