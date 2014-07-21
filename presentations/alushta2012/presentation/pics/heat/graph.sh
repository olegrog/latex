#!/usr/bin/gnuplot

set term epslatex standalone size 4,3 font 11 color dashed
set out "graph.tex"
set key center bottom

set xlabel '$\mathrm{Kn}$' offset graph 0.48,0.17
set ylabel '$\displaystyle\frac{\hat{q}_x}{\hat{T}_1-\hat{T}_2}$' offset graph 0.32,0.46 rotate by 0

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

asympt(x) = x<100 ? kappa*x/(1+sqrt(pi)*d*x) : 1/0
hydro(x)  = x<3 ? kappa*x : 1/0
plot \
	"qflow.txt" using ($1/$2):($3/2) title "PMDO (HS)" w lp lw 2 lt 1 pt 6, \
	free title "Free molecular" lt 2 lc 0 lw 2, \
	asympt(x) title "Asymptotic linear theory (HS)" lt 2 lc 3 lw 2, \
	hydro(x) title "Navier--Stokes" lt 4 lc 0 lw 2, \
	"sone.txt" using 1:($2/2) title "Exact solution (HS)" w lp lw 2 lt 1 lc 3 pt 3, \
	"bgkw.txt" using (1./(3*eta)/$1):($2*free) title "BGK kinetics model" w p lw 2 lt 1 lc 2 pt 1 ps 2
