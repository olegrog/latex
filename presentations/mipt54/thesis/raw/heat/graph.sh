#!/usr/bin/env gnuplot

set term epslatex standalone size 4,3 font 11 dashed
set out "graph.tex"
set key center bottom

set xlabel '$\mathrm{Kn}$' offset graph 0.48,0.18
set ylabel '$-\displaystyle\frac{\hat{q}_x}{\hat{T}_1-\hat{T}_2}$' offset graph 0.32,0.46 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 4
set bmargin 2
set logscale xy
#set grid

set xrange [0.006:180]
set yrange [0.01:1.2]

set label "free molecular" at .5,.7
set label "fluid dynamics" at .03,.1 rotate by 58

kappa = 2.12947
eta = 0.562773
d = 2.4001
free = 1/sqrt(pi) 

asympt(x) = x<5 ? kappa*x/(1+sqrt(pi)*d*x) : 1/0
hydro(x)  = x<3 ? kappa*x : 1/0
plot \
	"qflow.txt" using ($1/$2):($3/2) title 'projection method' w p lw 2 pt 6, \
	free notitle lt 1 lw 2, \
	asympt(x) title 'asymptotic solution' lt 2 lw 2, \
	hydro(x) notitle lt 1 lw 2, \
	"sone.txt" using 1:($2/2) title 'LBE solution' w lp lw 2 lt 1 pt 3
