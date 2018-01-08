#!/usr/bin/env gnuplot

set term epslatex standalone size 5,3.75 font 11 color dashed
set out "qflow.tex"
set key center bottom
set colors classic
set format y "%g"

set xlabel '$\mathrm{Kn}$' offset graph 0.48,0.14
set ylabel '$-\displaystyle2\int_0^{\frac12}\frac{q_x}{\Delta{v}}\:\mathrm{d}y$' offset graph 0.68,0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 7
set bmargin 2
set logscale xy
#set grid

set xrange [0.005:180]
set yrange [2e-5:0.03]

eta = 0.562773
k0 = -1.2540
free = 1/sqrt(pi) 
error = 1.0069

plot \
	"data.txt" using 1:($3) title "Projection DVM" w lp lw 2 dt 1 pt 6, \
	"sone.txt" using ($1*2/sqrt(pi)):($4) title "LBE" w lp lw 2 dt 1 lc 3 pt 3
