#!/usr/bin/env gnuplot

set term epslatex standalone size 5,3.75 font 11 color dashed
set out "graph.tex"
set key center bottom
set colors classic

set xlabel '$\mathrm{Kn}$' offset graph 0.48,0.14
set ylabel '$-\displaystyle\frac{p_{xy}}{\Delta{v}}$' offset graph 0.34,0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 4
set bmargin 2
set logscale xy
#set grid

set xrange [0.005:180]
set yrange [0.006:1.1]

eta = 0.562773
k0 = -1.2540
free = 1/sqrt(pi) 
error = 1.0069

asympt(x) = x<3 ? 2*eta*x/(1-sqrt(pi)*k0*x) : 1/0
hydro(x)  = x<3 ? 2*eta*x : 1/0
plot \
	"data.txt" using 1:(-$4) title "Projection DVM" w lp lw 2 dt 1 pt 6, \
	free title "Free molecular" dt 2 lc 0 lw 2, \
	asympt(x) title "Asymptotic" dt 2 lc 3 lw 2, \
	hydro(x) title "Nonslip Navier--Stokes" dt 1 lc 0 lw 2 , \
	"sone.txt" using ($1*2/sqrt(pi)):($2/2) title "LBE" w lp lw 2 dt 1 lc 3 pt 3 , \
	"bgkw.txt" using (1./(2*eta)/$1):($2*free) title "BGK" w p lw 2 dt 1 lc 2 pt 1 ps 2
