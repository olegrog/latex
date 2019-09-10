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
set clip two

set lmargin 6
set bmargin 2

set xrange [0:.5]
set yrange [*:.5]

U = 0.02
our_width = 2
our_point = 6
exact_width = 1
exact_point = 5
exact_color = 0
factor = sqrt(2)

x_0 = <buffer>
if (<hybrid>) set arrow from x_0, graph 0 to x_0, graph 1 nohead dt 4 lc 0 lw 3

plot \
    "<name>.txt" using 1:($2/U) title '$v_x/\Delta v$' w lp lw our_width pt our_point, \
    "<name>.txt" using 1:(-<shear_factor>*$3/U/<kn>) title sprintf('$-%gp_{xy}/k\Delta v$', <shear_factor>) w lp lw our_width pt our_point, \
    "<name>.txt" using 1:(-<qflow_factor>*$4/U/<kn>) title sprintf('$-%gq_x/k\Delta v$', <qflow_factor>) w lp lw our_width pt our_point, \
    1/0 title 'benchmark' w lp lc exact_color lw exact_width pt exact_point, \
    "k<kn>-my.txt" using 1:2 notitle w l lc exact_color lw exact_width, \
    "k<kn>-my.txt" using 1:(<shear_factor>*$3/<kn>/factor) notitle w l lc exact_color lw exact_width, \
    "k<kn>-my.txt" using 1:(<qflow_factor>*$4/<kn>) every 1 notitle w l lc exact_color lw exact_width, \
    "k<kn>.txt" using 1:2 notitle w p pt exact_point lc exact_color, \
    "k<kn>.txt" using 1:(-<shear_factor>*$3/<kn>*factor) notitle w p pt exact_point lc exact_color
