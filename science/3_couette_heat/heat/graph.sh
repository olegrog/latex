#!/usr/bin/env gnuplot

set term epslatex standalone size 5,3.75 font 11 color dashed
set out "graph.tex"
set key center bottom

set xlabel '$\mathrm{Kn}$' offset graph 0.48,0.14
set ylabel '$-\displaystyle\frac{\hat{q}_x}{\hat{T}_1-\hat{T}_2}$' offset graph 0.24,0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 4
set bmargin 2
set logscale xy
#set grid

set xrange [0.006:180]
set yrange [0.01:1.2]

kappa = 2.12947
eta = 0.562773
d = 2.4001
free = 1/sqrt(pi) 
error = 1.0069

asympt(x) = x<3 ? kappa*x/(1+sqrt(pi)*d*x) : 1/0
hydro(x)  = x<3 ? kappa*x : 1/0
plot \
	"qflow.txt" using ($1/$2):($3/2) title "Tchremissine solution" w lp lw 2 lt 1 pt 6, \
	free title "Free molecular solution" lt 2 lc 0 lw 2, \
	asympt(x) title "Asymptotic solution" lt 2 lc 3 lw 2, \
	hydro(x) title "Hydrodynamic solution" lt 1 lc 0 lw 2, \
	"sone.txt" using 1:($2/2) title "LBE solution" w lp lw 2 lt 1 lc 3 pt 3, \
	"bgkw.txt" using (1./(3*eta)/$1):($2*free) title "BGKW solution" w p lw 2 lt 1 lc 2 pt 1 ps 2
