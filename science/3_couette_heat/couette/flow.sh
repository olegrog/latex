#!/usr/bin/env gnuplot

set term epslatex standalone size 5,3.75 font 11 color dashed
set out "flow.tex"
set key center bottom
set colors classic

set xlabel '$\mathrm{Kn}$' offset graph 0.48,0.14
set ylabel '$\displaystyle2\int_0^{\frac12}\frac{v_y}{\Delta{v}}\:\mathrm{d}y$' offset graph 0.58,0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 5
set bmargin 2
set logscale xy
#set grid

set xrange [0.005:180]
set yrange [0.002:0.75]

eta = 0.562773
k0 = -1.2540
free = 1/sqrt(pi) 
error = 1.0069

plot \
	"data.txt" using 1:(0.5-$2) title "Projection DVM" w lp lw 2 dt 1 pt 6, \
	"sone.txt" using ($1*2/sqrt(pi)):($3) title "LBE" w lp lw 2 dt 1 lc 3 pt 3
