#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set out "<name>.tex"
set key left top

set xlabel '$y$' offset graph 0.5, 0.15
set colors classic

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2

set xrange [0:.5]
set yrange [0:.6]

U = 0.02
our_width = 1
our_point = 6
exact_point = 5
exact_color = 0

plot \
    "<name>.txt" using 1:($2/U) title '$v_x/\Delta v$' w lp lw our_width pt our_point, \
    "<name>.txt" using 1:(-$3/U) title '$-p_{xy}/\Delta v$' w lp lw our_width pt our_point, \
    "<name>.txt" using 1:(-$4/U/U) title '$-q_x/(\Delta v)^2$' w lp lw our_width pt our_point, \
    "k1e-1.txt" using 1:2 notitle w p pt exact_point lc exact_color, \
    "k1e-1.txt" using 1:(-2*$3) notitle w p pt exact_point lc exact_color
