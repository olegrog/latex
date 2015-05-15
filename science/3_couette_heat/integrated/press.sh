#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set out "press.tex"
set key center bottom maxrows 3 width -10

set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.14
set ylabel '$\displaystyle\int_0^\frac12 \frac{P\mathrm{d}y}{(\Delta v)^2} - \frac{P_{N\!S}^*}{(\Delta v)^2} - \left(\frac1{12} - \frac{P_{N\!S}^*}{(\Delta v)^2}\right)\frac{k}{1+k}$' \
    offset graph 0.5, 0.45 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy

set xrange [0.007:150]
set yrange [4e-4:4e-2]

k = sqrt(pi)/2
free = 1./12
P1 = 2.19295e-02
P2 = 4.33624e-02
P5 = 1.03237e-01

base2(x,y,v,f) = y/v - f/v - (free - f/v)*x/(1+x)
base3(x,y,v,f) = (y*v-0.5)/v/v - f/v - (free - f/v)*x/(1+x)

filter(x) = x < 0.1 ? x : 1/0

set macros
dummy = "NaN title ' ' lt -3"

plot \
    NaN title "Projection method"    w lp lw 2 lt 1 lc 1 pt 6, \
    NaN title "DSMC"             w lp lw 2 lt 1 lc 2 pt 4, \
    NaN title "Asymptotic"       w l lw 2 lt 1 lc 3, \
    NaN title "$\\Delta v=1$"    w l lw 2 lt 3 lc 0, \
    NaN title "$\\Delta v=2$"    w l lw 2 lt 4 lc 0, \
    NaN title "$\\Delta v=5$"    w l lw 2 lt 5 lc 0, \
    "asym-1.0.txt" using (filter($1)):(base2($1*k,$10,1,P1)) notitle w l lw 2 lt 3 lc 3, \
    "asym-2.0.txt" using (filter($1)):(base2($1*k,$10,2,P2)) notitle w l lw 2 lt 4 lc 3, \
    "asym-5.0.txt" using (filter($1)):(base2($1*k,$10,5,P5)) notitle w l lw 2 lt 5 lc 3, \
    "data-1.0.txt" using 1:(base2($1*k,($6+$7+$8)/3,1,P1))   notitle w lp lw 2 lt 3 pt 6 lc 1, \
    "data-2.0.txt" using 1:(base2($1*k,($6+$7+$8)/3,2,P2))   notitle w lp lw 2 lt 4 pt 6 lc 1, \
    "data-5.0.txt" using 1:(base2($1*k,($6+$7+$8)/3,5,P5))   notitle w lp lw 2 lt 5 pt 6 lc 1, \
    "dsmc-1.0.txt" using 1:(base3($1*k,($6+$7+$8)/3/2,1,P1)) notitle w lp lw 2 lt 3 pt 4 lc 2, \
    "dsmc-2.0.txt" using 1:(base3($1*k,($6+$7+$8)/3/2,2,P2)) notitle w lp lw 2 lt 4 pt 4 lc 2, \
    "dsmc-5.0.txt" using 1:(base3($1*k,($6+$7+$8)/3/2,5,P5)) notitle w lp lw 2 lt 5 pt 4 lc 2

