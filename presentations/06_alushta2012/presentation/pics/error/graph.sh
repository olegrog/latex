#!/usr/bin/env gnuplot

set term epslatex standalone size 4,3 font 11 color dashed
set out "graph.tex"
set key left bottom

set xlabel '$V_\Omega$' offset graph 0.49,0.18
set ylabel '$\varepsilon$' offset graph 0.28,0.5 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror
#set mxtics 10

set lmargin 5
set bmargin 2
set logscale xy
#set grid

set xrange [4/pi*9**3:4/pi*45**3]
set yrange [5e-4:0.1]

plot \
	"data.txt" using (4/pi*$1**3):2 title "Heat flow error" w lp lw 2 lt 1 pt 6, \
	"data.txt" using (4/pi*$1**3):3 title "Heat conductivity error" w lp lw 2 lt 1 lc 3 pt 3
#	"data.txt" using 1:4 notitle w p lw 2 lt 1 lc 3 pt 3
