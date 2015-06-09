#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set out "qflow.tex"
set key center bottom maxrows 6 width -10

set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.14
set ylabel '$-\displaystyle\int_0^\frac12 \frac{q_x}{\Delta v} \mathrm{d}y$' \
    offset graph 0.3, 0.45 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy

set xrange [0.007:150]
set yrange [2e-5:0.6]

free = 1/sqrt(pi)
k = sqrt(pi)/2
gamma2 = 1.922284066

filter(x) = x > 0.12 ? x : 1/0
f_asym(x) = x < 0.15 ? x : 1/0

set macros
dummy = "NaN title ' ' lt -3"

plot \
    @dummy, @dummy, @dummy, @dummy, @dummy, \
    "bkw.txt" using (filter($1/k)/gamma2):(-$4+$3/2) title "BKW" w l lw 3 lt 1 lc 5, \
    "asym-bkw.txt" using (f_asym($1)/gamma2):(-$3) notitle w l lw 3 lt 1 lc 5, \
    "sone.txt" using ($1/k):($4/2) title "Ohwada et al." w lp lw 3 lt 1 lc 0 pt 3 ps 1.5, \
    "asym-hs.txt" using (f_asym($1)):(-$3) notitle w l lw 3 lt 1 lc 0, \
    NaN title "Projection method"    w lp lw 2 lt 1 lc 1 pt 6, \
    NaN title "DSMC"             w lp lw 2 lt 1 lc 2 pt 4, \
    NaN title "Asymptotic"       w l lw 2 lt 1 lc 3, \
    NaN title "$\\Delta v\\to0$" w l lw 2 lt 1 lc 0, \
    NaN title "$\\Delta v=0.1$"  w l lw 2 lt 2 lc 0, \
    NaN title "$\\Delta v=1$"    w l lw 2 lt 3 lc 0, \
    NaN title "$\\Delta v=2$"    w l lw 2 lt 4 lc 0, \
    NaN title "$\\Delta v=5$"    w l lw 2 lt 5 lc 0, \
    "asym-0.1.txt" using 1:(-$4) notitle w l lw 2 lt 2 lc 3, \
    "asym-1.0.txt" using 1:(-$4) notitle w l lw 2 lt 3 lc 3, \
    "asym-2.0.txt" using 1:(-$4) notitle w l lw 2 lt 4 lc 3, \
    "asym-5.0.txt" using 1:(-$4) notitle w l lw 2 lt 5 lc 3, \
    "data-0.1.txt" using 1:(-$4) notitle w lp lw 2 lt 2 pt 6 lc 1, \
    "data-1.0.txt" using 1:(-$4) notitle w lp lw 2 lt 3 pt 6 lc 1, \
    "data-2.0.txt" using 1:(-$4) notitle w lp lw 2 lt 4 pt 6 lc 1, \
    "data-5.0.txt" using 1:(-$4) notitle w lp lw 2 lt 5 pt 6 lc 1, \
    "dsmc-0.1.txt" using 1:(-$4) notitle w lp lw 2 lt 2 pt 4 lc 2, \
    "dsmc-1.0.txt" using 1:(-$4) notitle w lp lw 2 lt 3 pt 4 lc 2, \
    "dsmc-2.0.txt" using 1:(-$4) notitle w lp lw 2 lt 4 pt 4 lc 2, \
    "dsmc-5.0.txt" using 1:(-$4) notitle w lp lw 2 lt 5 pt 4 lc 2

