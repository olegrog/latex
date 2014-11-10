#!/usr/bin/env gnuplot

set term epslatex standalone size 5,3.75 font 11 color dashed
set out "graph.tex"
set key center bottom

set xlabel '$\mathrm{Kn}$' offset graph 0.48,0.14
set ylabel '$-\displaystyle\frac{P_{xy}}{\Delta v}$' offset graph 0.22,0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 4
set bmargin 2
set logscale xy
#set grid

set xrange [0.005:180]
set yrange [0.006:1.1]

gamma1 = 1.270042427
k0_hs = -1.2540
k0_bkw = -1.01619
k = sqrt(pi)/2

free = 1/sqrt(pi) 

asym_hs(x) = x<3 ? gamma1*x/(1-2*k0_hs*x) : 1/0
asym_bkw(x) = x<3 ? x/(1-2*k0_bkw*x) : 1/0
fluid(x)  = x<3 ? gamma1*x : 1/0
near_free(x) = x>2 ? free * (1-k/x) : 1/0

filter(x) = x > 0.12 ? x : 1/0
f_asym(x) = x < 0.15 ? x : 1/0

plot \
	"data.txt" using 1:(-$2) title "Tchremissine solution" w lp lw 2 lt 1 pt 6, \
	free title "Free molecular solution" lt 2 lc 0 lw 2, \
	"sone.txt" using ($1/k):($2/2) title "Exact solution" w lp lw 2 lt 1 lc 3 pt 3 , \
	asym_hs(k*f_asym(x)) notitle lt 2 lc 3 lw 2, \
	fluid(k*x) title "Continuum solution" lt 5 lc 0 lw 2 , \
	"bkw.txt" using (filter($1/k)/gamma1):(-$2) title "BKW solution" w l lw 1 lt 1 lc 2, \
	asym_bkw(k*f_asym(x)*gamma1) notitle w l lw 1 lt 2 lc 2

#	"bgkw.txt" using (1./k/$1/gamma1):($2*free) title "BKW solution" w p lw 2 lt 1 lc 2 pt 1 ps 1.5, \
#    near_free(x) title "BKW high Kn" w l lw 2 lt 4 lc 5
