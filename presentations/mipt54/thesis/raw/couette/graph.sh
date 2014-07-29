#!/usr/bin/env gnuplot

set term epslatex standalone size 4,3 font 11 dashed
set out "graph.tex"
set key center bottom

set xlabel '$\mathrm{Kn}$' offset graph 0.48,0.18
set ylabel '$-\displaystyle\frac{\hat{p}_{xy}}{\hat{U}}$' offset graph 0.28,0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 4
set bmargin 2
set logscale xy
#set grid

set xrange [0.005:180]
set yrange [0.006:1.1]

set label "free molecular" at 1,.7
set label "fluid dynamics" at .05,.08 rotate by 57

eta = 0.562773
k0 = -1.2540
free = 1/sqrt(pi) 

asympt(x) = x<5 ? 2*eta*x/(1-sqrt(pi)*k0*x) : 1/0
hydro(x)  = x<3 ? 2*eta*x : 1/0
plot \
	"data.txt" using 1:(-$4) title 'projection method' w p lw 2 pt 6, \
 	free notitle lt 1 lw 2, \
 	asympt(x) title 'asymptotic solution' lt 2 lw 2, \
 	hydro(x) notitle lt 1 lw 2 , \
 	"sone.txt" using ($1*2/sqrt(pi)):($2/2) title 'LBE solution' w lp lw 2 lt 1 pt 3
