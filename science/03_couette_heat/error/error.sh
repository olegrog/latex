#!/usr/bin/env gnuplot

set term epslatex standalone size 4,3 font 11 color dashed
set out "error.tex"
set key left bottom
set colors classic

set xlabel '$N_R$' offset graph 0.48,0.17
set ylabel '$\varepsilon$' offset graph 0.25,0.5 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 5
set bmargin 2
set logscale xy
set xtics add (16, 26, 40)
#set grid

set xrange [9:45]
set yrange [5e-4:0.1]

plot \
	"data.txt" using 1:2 title "Heat flow error" w lp lw 2 lt 1 pt 6, \
	"data.txt" using 1:3 title "Heat conductivity error" w lp lw 2 lt 1 lc 3 pt 3
#	"data.txt" using 1:4 notitle w p lw 2 lt 1 lc 3 pt 3
