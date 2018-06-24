#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set out "<name>.tex"
set key left top

set xlabel '$y$' offset graph 0.48, 0.15
set colors classic

set border 3
set xtics nomirror
set ytics nomirror
set grid

set lmargin 6
set bmargin 2

set xrange [0:.5]
set yrange [0:.5]

U = 0.02
our_width = 2
our_point = 6
exact_width = 1
exact_point = 5
exact_color = 0
factor = sqrt(2)
k = 40

x_0 = <buffer>
if (<hybrid>) set arrow from x_0, graph 0 to x_0, graph 1 nohead dt 4 lc 0 lw 3

plot \
    "<name>.txt" using 1:($2/U) title '$v_x/\Delta v$' w lp lw our_width pt our_point, \
    "<name>.txt" using 1:(-$3/U/factor) title '$-p_{xy}/\Delta v$' w lp lw our_width pt our_point, \
    "<name>.txt" using 1:(-k*$4/U) title sprintf('$-%dq_x/\Delta v$', k) w lp lw our_width pt our_point, \
    1/0 title 'benchmark' w lp lc exact_color lw exact_width pt exact_point, \
    "k1e-1-my.txt" using 1:2 notitle w l lc exact_color lw exact_width, \
    "k1e-1-my.txt" using 1:($4/factor) notitle w l lc exact_color lw exact_width, \
    "k1e-1-my.txt" using 1:(k*($5)) notitle w l lc exact_color lw exact_width, \
    "k1e-1.txt" using 1:2 notitle w p pt exact_point lc exact_color, \
    "k1e-1.txt" using 1:(-2*$3/factor) notitle w p pt exact_point lc exact_color
