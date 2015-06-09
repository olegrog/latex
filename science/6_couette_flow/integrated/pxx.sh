#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set out "pxx.tex"
set key center bottom maxrows 3 width -10

set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.14
set ylabel '$\displaystyle\int_0^\frac12 \frac{P_{yy}-P_{xx}}{(\Delta v)^2}\mathrm{d}y + \frac{k}{2+4k}$' \
    offset graph 0.36, 0.45 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy

set xrange [0.007:150]
set yrange [3e-3:0.1]

k = sqrt(pi)/2
gamma1 = 1.270042427
gamma2 = 1.922284066
gamma3 = 1.947906335
gamma7 = 0.189201
gamma8 = 1.4963
gamma9 = 1.6357
alpha = 1.25*gamma2/gamma1

free = 0.25
asym(x) = x < .12 ? 0.5*x**2*((gamma3+gamma7/4)/alpha+2*gamma9) : 1/0

base(x,y) = free*x/(.5 + x) - y
filter(x) = x < 0.2 ? x : 1/0

set macros
dummy = "NaN title ' ' lt -3"

plot \
    NaN title "Projection method"    w lp lw 2 lt 1 lc 1 pt 6, \
    NaN title "DSMC"             w lp lw 2 lt 1 lc 2 pt 4, \
    NaN title "Asymptotic"       w l lw 2 lt 1 lc 3, \
    NaN title "$\\Delta v=1$"    w l lw 2 lt 3 lc 0, \
    NaN title "$\\Delta v=2$"    w l lw 2 lt 4 lc 0, \
    NaN title "$\\Delta v=5$"    w l lw 2 lt 5 lc 0, \
    "asym-1.0.txt" using (filter($1)):(base($1*k,($6-$7)/1)) notitle w l lw 2 lt 3 lc 3, \
    "asym-2.0.txt" using (filter($1)):(base($1*k,($6-$7)/2)) notitle w l lw 2 lt 4 lc 3, \
    "asym-5.0.txt" using (filter($1)):(base($1*k,($6-$7)/5)) notitle w l lw 2 lt 5 lc 3, \
    "data-1.0.txt" using 1:(base($1*k,($6-$7)/1))   notitle w lp lw 2 lt 3 pt 6 lc 1, \
    "data-2.0.txt" using 1:(base($1*k,($6-$7)/2))   notitle w lp lw 2 lt 4 pt 6 lc 1, \
    "data-5.0.txt" using 1:(base($1*k,($6-$7)/5))   notitle w lp lw 2 lt 5 pt 6 lc 1, \
    "dsmc-1.0.txt" using 1:(base($1*k,($6-$7)/1/2)) notitle w lp lw 2 lt 3 pt 4 lc 2, \
    "dsmc-2.0.txt" using 1:(base($1*k,($6-$7)/2/2)) notitle w lp lw 2 lt 4 pt 4 lc 2, \
    "dsmc-5.0.txt" using 1:(base($1*k,($6-$7)/5/2)) notitle w lp lw 2 lt 5 pt 4 lc 2

#    "asym-0.1.txt" using 1:(base($1*k,($6-$7)/0.1))     notitle w l lw 2 lt 2 lc 3, \
#    "data-0.1.txt" using 1:(base($1*k,($6-$7)/0.1))     notitle w lp lw 2 lt 2 pt 6 lc 1, \
#    "dsmc-0.1.txt" using 1:(base($1*k,($6-$7)/0.1/2))   notitle w lp lw 2 lt 2 pt 4 lc 2, \
#    base(k*x, asym(k*x)) title "Asymptotic" lw 3 lt 7 lc 3, \
