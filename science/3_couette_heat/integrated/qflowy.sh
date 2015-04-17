#!/usr/bin/env gnuplot

set term epslatex standalone size 6, 4 font 11 color dashed
set out "qflowy.tex"
set key center bottom maxrows 3 width -3

set xlabel '$\mathrm{Kn}$' offset graph 0.48,0.12
set ylabel '$\displaystyle\frac{\gamma_1 k}{8(1+2k^{5/3})} - \frac{Q_{N\!S}}{Q_{N\!S}^*}\int_0^\frac12 \frac{q_y}{(\Delta v)^2}\mathrm{d}y$' \
    offset graph 0.4,0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy

set xrange [0.005:180]
set yrange [8e-5:5e-2]

gamma1 = 1.270042427
k = sqrt(pi)/2
k0 = -1.2540 
fluid(x) = x < 1 ? gamma1*x/8 : 1/0
#asym(x) = x < 0.15 ? fluid(x) * (1 + 4*k0*x) : 1/0
asym(x) = x < 0.15 ? fluid(x) / (1 - 2*k0*x)**2 : 1/0
base(x,y) = gamma1*x/8/( 1+2*x**(5./3) ) - y
base2(x,y,v,c) = gamma1*x/8/(1+2*x**(5./3) ) - y*v/c/4

set macros
dummy = "NaN title ' ' lt -3"

plot \
    NaN title "Tcheremissine"    w lp lw 2 lt 1 lc 1 pt 6, \
    NaN title "DSMC"             w lp lw 2 lt 1 lc 2 pt 4, \
    NaN title "Asymptotic"       w l lw 2 lt 1 lc 3, \
    NaN title "$\\Delta v=1$"    w l lw 2 lt 3 lc 0, \
    NaN title "$\\Delta v=2$"    w l lw 2 lt 4 lc 0, \
    NaN title "$\\Delta v=5$"    w l lw 2 lt 5 lc 0, \
    "asym-1.0.txt" using 1:(base2($1*k,$5,1,1.0217401*0.24867867))  notitle w l lw 2 lt 3 lc 3, \
    "asym-2.0.txt" using 1:(base2($1*k,$5,2,2.167795*0.49055138))   notitle w l lw 2 lt 4 lc 3, \
    "asym-5.0.txt" using 1:(base2($1*k,$5,5,7.1917*1.163882))       notitle w l lw 2 lt 5 lc 3, \
    "data-1.0.txt" using 1:(base2($1*k,$5,1,1.0217401*0.24867867))  notitle w lp lw 2 lt 3 pt 6 lc 1, \
    "data-2.0.txt" using 1:(base2($1*k,$5,2,2.167795*0.49055138))   notitle w lp lw 2 lt 4 pt 6 lc 1, \
    "data-5.0.txt" using 1:(base2($1*k,$5,5,7.1917*1.163882))       notitle w lp lw 2 lt 5 pt 6 lc 1, \
    "dsmc-1.0.txt" using 1:(base2($1*k,$5,1,1.0217401*0.24867867))  notitle w lp lw 2 lt 3 pt 4 lc 2, \
    "dsmc-2.0.txt" using 1:(base2($1*k,$5,2,2.167795*0.49055138))   notitle w lp lw 2 lt 4 pt 4 lc 2, \
    "dsmc-5.0.txt" using 1:(base2($1*k,$5,5,7.1917*1.163882))       notitle w lp lw 2 lt 5 pt 4 lc 2

#    "asym-0.1.txt" using 1:(base2($1*k,$5,0.1,0.1000220*0.0249986)) notitle w l lw 2 lt 2 lc 3, \
#    "data-0.1.txt" using 1:(base2($1*k,$5,0.1,0.1000220*0.0249986)) notitle w lp lw 2 lt 2 pt 6 lc 1, \
#    "dsmc-0.1.txt" using 1:(base2($1*k,$5,0.1,0.1000220*0.0249986)) notitle w lp lw 2 lt 2 pt 4 lc 2, \
#    base(k*x, fluid(k*x)) title "Navier--Stokes" lt 1 lc 0 lw 2 , \
#    base(k*x, asym(k*x)) title "Asymptotic" lw 3 lt 7 lc 3, \
