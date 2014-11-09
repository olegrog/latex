#!/usr/bin/env gnuplot

set term epslatex standalone size 5,3.75 font 11 color dashed
set out "qflow.tex"
set key center bottom

set xlabel '$\mathrm{Kn}$' offset graph 0.48,0.14
set ylabel '$-\displaystyle\int_0^\frac12\frac{q_y}{\Delta v}\mathrm{d}x$' offset graph 0.32,0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy

set xrange [0.005:180]
set yrange [2e-5:0.05]

free = 1/sqrt(pi)
k = sqrt(pi)/2
gamma2 = 1.922284066

filter(x) = x > 0.15 ? x : 1/0
f_asym(x) = x < 0.2 ? x : 1/0

plot \
	"data.txt" using 1:($4/2) title "Tchremissine solution" w lp lw 2 lt 1 pt 6, \
	"sone.txt" using ($1/k):($4/2) title "Exact solution" w lp lw 2 lt 1 lc 3 pt 3, \
	"asym-hs.txt" using (f_asym($1)):(-$3) notitle w l lw 2 lt 2 lc 3, \
	"bkw.txt" using (filter($1/k)/gamma2):(-$4+$3/2) title "BKW solution" w l lw 2 lt 1 lc 2, \
	"asym-bkw.txt" using (f_asym($1)/gamma2):(-$3) notitle w l lw 2 lt 2 lc 2, \

#    "free.txt" using ($1/gamma2):(-$3) title "BKW high Kn" w l lw 2 lt 4 lc 5
