#!/usr/bin/env gnuplot

set term epslatex standalone size 5,3.75 font 11 color dashed
set out "flow.tex"
set key left bottom

set xlabel '$\mathrm{Kn}$' offset graph 0.48,0.14
set ylabel '$\displaystyle\int_0^\frac12 \frac{v_x}{\Delta v} \mathrm{d}y$' offset graph 0.26,0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 5
set bmargin 2
set logscale xy

set xrange [0.005:180]
set yrange [0.002:0.3]

cont = 0.125
k = sqrt(pi)/2
gamma1 = 1.270042

filter(x) = x > 0.12 ? x : 1/0
f_asym(x) = x < 0.15 ? x : 1/0

plot \
    cont title "Continuum solution" w l lw 2 lc 0 lt 5, \
	"data.txt" using 1:((0.5-$3)/2) title "Tchremissine solution" w lp lw 2 lt 1 pt 6, \
	"sone.txt" using ($1/k):($3/2) title "Exact solution" w lp lw 2 lt 1 lc 3 pt 3, \
	"asym-hs.txt" using (f_asym($1)):2 notitle w l lw 2 lt 2 lc 3, \
	"bkw.txt" using (filter($1/k)/gamma1):3 title "BKW solution" w l lw 1 lt 1 lc 2, \
	"asym-bkw.txt" using (f_asym($1)/gamma1):2 notitle w l lw 1 lt 2 lc 2
    
#    "free.txt" using ($1/gamma1):2 title "BKW high Kn" w l lw 2 lt 4 lc 5
