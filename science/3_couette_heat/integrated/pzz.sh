#!/usr/bin/env gnuplot

set term epslatex standalone size 6, 4 font 11 color dashed
set out "pzz.tex"
set key center bottom maxrows 3 width -3

set xlabel '$\mathrm{Kn}$' offset graph 0.48, 0.13
set ylabel '$\displaystyle\int_0^\frac12 \frac{P_{zz}-P_{yy}}{(\Delta v)^2} \mathrm{d}y$' \
    offset graph 0.3, 0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy

set xrange [0.005:180]
set yrange [1e-5:5e-2]

k = sqrt(pi)/2
gamma1 = 1.270042427
gamma2 = 1.922284066
gamma3 = 1.947906335
gamma7 = 0.189201
gamma8 = 1.4963
gamma9 = 1.6357
alpha = 1.25*gamma2/gamma1

asym(x) = x < .12 ? 0.5*x**2*(gamma3/alpha+2*(gamma9-gamma8)) : 1/0

set macros
dummy = "NaN title ' ' lt -3"

plot \
    NaN title "Tcheremissine"    w lp lw 2 lt 1 lc 1 pt 6, \
    NaN title "DSMC"             w lp lw 2 lt 1 lc 2 pt 4, \
    NaN title "Asymptotic"       w l lw 2 lt 1 lc 3, \
    NaN title "$\\Delta v=1$"    w l lw 2 lt 3 lc 0, \
    NaN title "$\\Delta v=2$"    w l lw 2 lt 4 lc 0, \
    NaN title "$\\Delta v=5$"    w l lw 2 lt 5 lc 0, \
    "asym-0.1.txt" using 1:(($8-$7)/0.1)    notitle w l lw 2 lt 2 lc 3, \
    "asym-1.0.txt" using 1:(($8-$7)/1)      notitle w l lw 2 lt 3 lc 3, \
    "asym-2.0.txt" using 1:(($8-$7)/2)      notitle w l lw 2 lt 4 lc 3, \
    "asym-5.0.txt" using 1:(($8-$7)/5)      notitle w l lw 2 lt 5 lc 3, \
    "data-0.1.txt" using 1:(($8-$7)/0.1)    notitle w lp lw 2 lt 2 pt 6 lc 1, \
    "data-1.0.txt" using 1:(($8-$7)/1)      notitle w lp lw 2 lt 3 pt 6 lc 1, \
    "data-2.0.txt" using 1:(($8-$7)/2)      notitle w lp lw 2 lt 4 pt 6 lc 1, \
    "data-5.0.txt" using 1:(($8-$7)/5)      notitle w lp lw 2 lt 5 pt 6 lc 1, \
    "dsmc-0.1.txt" using 1:(($8-$7)/0.1/2)  notitle w lp lw 2 lt 2 pt 4 lc 2, \
    "dsmc-1.0.txt" using 1:(($8-$7)/1/2)    notitle w lp lw 2 lt 3 pt 4 lc 2, \
    "dsmc-2.0.txt.new" using 1:(($8-$7)/2/2)  notitle w lp lw 2 lt 4 pt 6 lc 4, \
    "dsmc-2.0.txt" using 1:(($8-$7)/2/2)    notitle w lp lw 2 lt 4 pt 4 lc 2, \
    "dsmc-5.0.txt" using 1:(($8-$7)/5/2)    notitle w lp lw 2 lt 5 pt 4 lc 2

#    asym(k*x) title "Asymptotic" lw 3 lt 7 lc 3, \
