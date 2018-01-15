#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set colors classic
set out "shear-overpeak.tex"
set key right bottom maxrows 5 width -7
set sample 101

set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.14
set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.14
set ylabel '$-\displaystyle\frac{P_{xy}}{\Delta v} - \frac1{\sqrt\pi}$' \
    offset graph 0.55, 0.45 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale x
set xzeroaxis

set xrange [1:150]
set yrange [-0.1:0.02]

gamma1 = 1.270042427
k = sqrt(pi)/2
Pxy0 = -6.35161e-01
Pxy1 = -6.48827e-01
Pxy2 = -6.88298e-01
Pxy5 = -9.13379e-01

free = -1/sqrt(pi)

base(x,y) = free - y
base2(x,y,f) = free - 2*y
base3(x,y,f) = free -y
filter(x) = x > 0.12 ? x : 1/0

set macros
dummy = "NaN title ' ' lc -3"

plot \
    "bkw.txt" using (filter($1/k)/gamma1):(base(filter($1)/gamma1,$2)) title "BGK" w l lw 3 dt 1 lc 5, \
    "sone.txt" using ($1/k):(base($1,-$2/2)) title "Sone et al." w lp lw 3 dt 1 lc 0 pt 3 ps 1.5, \
    NaN title "Projection DVM"    w lp lw 2 dt 1 lc 1 pt 6, \
    NaN title "DSMC"             w lp lw 2 dt 1 lc 2 pt 4, \
    @dummy, \
    NaN title "$\\Delta v\\to0$" w l lw 2 dt 1 lc 0, \
    NaN title "$\\Delta v=0.1$"  w l lw 2 dt 6 lc 0, \
    NaN title "$\\Delta v=1$"    w l lw 2 dt 2 lc 0, \
    NaN title "$\\Delta v=2$"    w l lw 2 dt 3 lc 0, \
    NaN title "$\\Delta v=5$"    w l lw 2 dt 4 lc 0, \
    "data-0.1.txt" using 1:(base2($1*k,$2,Pxy0)) notitle w lp lw 2 dt 6 pt 6 lc 1, \
    "data-1.0.txt" using 1:(base2($1*k,$2,Pxy1)) notitle w lp lw 2 dt 2 pt 6 lc 1, \
    "data-2.0.txt" using 1:(base2($1*k,$2,Pxy2)) notitle w lp lw 2 dt 3 pt 6 lc 1, \
    "data-5.0.txt" using 1:(base2($1*k,$2,Pxy5)) notitle w lp lw 2 dt 4 pt 6 lc 1, \
    "dsmc-0.1.txt" using 1:(base3($1*k,$2,Pxy0)) notitle w lp lw 2 dt 6 pt 4 lc 2, \
    "dsmc-1.0.txt" using 1:(base3($1*k,$2,Pxy1)) notitle w lp lw 2 dt 2 pt 4 lc 2, \
    "dsmc-2.0.txt" using 1:(base3($1*k,$2,Pxy2)) notitle w lp lw 2 dt 3 pt 4 lc 2, \
    "dsmc-5.0.txt" using 1:(base3($1*k,$2,Pxy5)) notitle w lp lw 2 dt 4 pt 4 lc 2

