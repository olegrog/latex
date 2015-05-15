#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set out "qflowy.tex"
set key center bottom maxrows 3 width -10

set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.14
set ylabel '$\displaystyle\int_0^\frac12 \frac{q_y\mathrm{d}y}{(\Delta v)^2}$' \
    offset graph 0.26, 0.45 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy

set xrange [0.007:150]
set yrange [1e-3:4e-2]

gamma1 = 1.270042427
k = sqrt(pi)/2
k0 = -1.2540
Q0 = 1.5878152e-02
Q1 = 1.6134940e-01
Q2 = 3.3764590e-01
Q5 = 1.0630691e+00

fluid(x) = x < 1 ? gamma1*x/8 : 1/0
#asym(x) = x < 0.15 ? fluid(x) * (1 + 4*k0*x) : 1/0
asym(x) = x < 0.15 ? fluid(x) / (1 - 2*k0*x)**2 : 1/0
base(x,y) = gamma1*x/8/( 1+2*x**(5./3) ) - y
#base2(x,y,v,c) = gamma1*x/8/(1+2*x**(5./3) ) - y*v/c/4
base2(x,y,v,c) = y/v # - c*gamma1/2/v/v * x/(1-2*k0*x)**2

set macros
dummy = "NaN title ' ' lt -3"

plot \
    NaN title "Projection method"    w lp lw 2 lt 1 lc 1 pt 6, \
    NaN title "DSMC"             w lp lw 2 lt 1 lc 2 pt 4, \
    NaN title "Asymptotic"       w l lw 2 lt 1 lc 3, \
    NaN title "$\\Delta v=1$"    w l lw 2 lt 3 lc 0, \
    NaN title "$\\Delta v=2$"    w l lw 2 lt 4 lc 0, \
    NaN title "$\\Delta v=5$"    w l lw 2 lt 5 lc 0, \
    "asym-1.0.txt" using 1:(base2($1*k,$5,1,Q1)) notitle w l lw 2 lt 3 lc 3, \
    "asym-2.0.txt" using 1:(base2($1*k,$5,2,Q2)) notitle w l lw 2 lt 4 lc 3, \
    "asym-5.0.txt" using 1:(base2($1*k,$5,5,Q5)) notitle w l lw 2 lt 5 lc 3, \
    "data-1.0.txt" using 1:(base2($1*k,$5,1,Q1)) notitle w lp lw 2 lt 3 pt 6 lc 1, \
    "data-2.0.txt" using 1:(base2($1*k,$5,2,Q2)) notitle w lp lw 2 lt 4 pt 6 lc 1, \
    "data-5.0.txt" using 1:(base2($1*k,$5,5,Q5)) notitle w lp lw 2 lt 5 pt 6 lc 1, \
    "dsmc-1.0.txt" using 1:(base2($1*k,$5,1,Q1)) notitle w lp lw 2 lt 3 pt 4 lc 2, \
    "dsmc-2.0.txt" using 1:(base2($1*k,$5,2,Q2)) notitle w lp lw 2 lt 4 pt 4 lc 2, \
    "dsmc-5.0.txt" using 1:(base2($1*k,$5,5,Q5)) notitle w lp lw 2 lt 5 pt 4 lc 2

#    "asym-0.1.txt" using 1:(base2($1*k,$5,0.1,Q0)) notitle w l lw 2 lt 2 lc 3, \
#    "data-0.1.txt" using 1:(base2($1*k,$5,0.1,Q0)) notitle w lp lw 2 lt 2 pt 6 lc 1, \
#    "dsmc-0.1.txt" using 1:(base2($1*k,$5,0.1,Q0)) notitle w lp lw 2 lt 2 pt 4 lc 2, \
#    base(k*x, fluid(k*x)) title "Navier--Stokes" lt 1 lc 0 lw 2 , \
#    base(k*x, asym(k*x)) title "Asymptotic" lw 3 lt 7 lc 3, \
