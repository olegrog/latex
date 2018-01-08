#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set colors classic
set out "shear.tex"
set key right bottom maxrows 5 width -7
set sample 101

set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.14
set ylabel '$\displaystyle \int_0^\frac12 \frac{P_{xy}}{\Delta v}\mathrm{d}y + \frac{P_{N\!Sxy}^*k}{\Delta{v}+2\sqrt\pi P_{N\!Sxy}^*k}$' \
    offset graph 1.31, 0.45 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale x #y
set xzeroaxis

set xrange [0.007:150]
#set yrange [1e-3:1.5e-1]
set yrange [-0.015:0.015]

gamma1 = 1.270042427
k0_hs = -1.2540
k0_bkw = -1.01619
k = sqrt(pi)/2
Pxy0 = -6.35161e-01
Pxy1 = -6.48827e-01
Pxy2 = -6.88298e-01
Pxy5 = -9.13379e-01

free = -0.5/sqrt(pi)
fluid(x) = -gamma1/2*x
asym_hs(x) = x < 3 ? -gamma1*x/(1-2*k0_hs*x) : 1/0
asym_bkw(x) = x < 3 ? -x/(1-2*k0_bkw*x) : 1/0
near_free(x) = free * x/(x+k)

base(x,y) = -fluid(x)/(1-2*sqrt(pi)*fluid(x)) + y/2
base2(x,y,f) = -f*x/(1-2*sqrt(pi)*f*x) + y
base3(x,y,f) = -f*x/(1-2*sqrt(pi)*f*x) + y/2
#base(x,y) = free*x/(1+x) - y/2
#base2(x,y,f) = free*x/(1+x) - y
#base3(x,y,f) = free*x/(1+x) - y/2
#base(x,y) = fluid(x)/(1-2*k0_hs*x) + (free+fluid(x)/2/k0_hs/x)*x/(1+x) - y/2
#base2(x,y,f) = f*x/(1-2*k0_hs*x) + (free+f/2/k0_hs)*x/(1+x) - y
#base3(x,y,f) = f*x/(1-2*k0_hs*x) + (free+f/2/k0_hs)*x/(1+x) - y/2

filter(x) = x > 0.12 ? x : 1/0
f_asym(x) = x < 0.15 ? x : 1/0
filter2(x,c) = x < c ? x : 1/0

set macros
dummy = "NaN title ' ' lc -3"

plot \
    "bkw.txt" using (filter($1/k)/gamma1):(base(filter($1)/gamma1,$2)) title "BGK" w l lw 3 dt 1 lc 5, \
    base(k*x, asym_bkw(k*f_asym(x)*gamma1)) notitle w l lw 3 dt 1 lc 5, \
    "sone.txt" using ($1/k):(base($1,-$2/2)) title "Sone et al." w lp lw 3 dt 1 lc 0 pt 3 ps 1.5, \
    base(k*x,asym_hs(k*f_asym(x))) notitle lw 3 dt 1 lc 0, \
    NaN title "Projection DVM"    w lp lw 2 dt 1 lc 1 pt 6, \
    NaN title "DSMC"             w lp lw 2 dt 1 lc 2 pt 4, \
    NaN title "Asymptotic"       w l lw 2 dt 1 lc 3, \
    NaN title "$\\Delta v\\to0$" w l lw 2 dt 1 lc 0, \
    NaN title "$\\Delta v=0.1$"  w l lw 2 dt 6 lc 0, \
    NaN title "$\\Delta v=1$"    w l lw 2 dt 2 lc 0, \
    NaN title "$\\Delta v=2$"    w l lw 2 dt 3 lc 0, \
    NaN title "$\\Delta v=5$"    w l lw 2 dt 4 lc 0, \
    @dummy, @dummy, @dummy, @dummy, @dummy, \
    "asym-0.1.txt" using (filter2($1,0.50)):(base2($1*k,$2,Pxy0)) notitle w l lw 2 dt 6 lc 3, \
    "asym-1.0.txt" using (filter2($1,0.18)):(base2($1*k,$2,Pxy1)) notitle w l lw 2 dt 2 lc 3, \
    "asym-2.0.txt" using (filter2($1,0.12)):(base2($1*k,$2,Pxy2)) notitle w l lw 2 dt 3 lc 3, \
    "asym-5.0.txt" using (filter2($1,0.06)):(base2($1*k,$2,Pxy5)) notitle w l lw 2 dt 4 lc 3, \
    "data-0.1.txt" using 1:(base2($1*k,$2,Pxy0)) notitle w lp lw 2 dt 6 pt 6 lc 1, \
    "data-1.0.txt" using 1:(base2($1*k,$2,Pxy1)) notitle w lp lw 2 dt 2 pt 6 lc 1, \
    "data-2.0.txt" using 1:(base2($1*k,$2,Pxy2)) notitle w lp lw 2 dt 3 pt 6 lc 1, \
    "data-5.0.txt" using 1:(base2($1*k,$2,Pxy5)) notitle w lp lw 2 dt 4 pt 6 lc 1, \
    "dsmc-0.1.txt" using 1:(base3($1*k,$2,Pxy0)) notitle w lp lw 2 dt 6 pt 4 lc 2, \
    "dsmc-1.0.txt" using 1:(base3($1*k,$2,Pxy1)) notitle w lp lw 2 dt 2 pt 4 lc 2, \
    "dsmc-2.0.txt" using 1:(base3($1*k,$2,Pxy2)) notitle w lp lw 2 dt 3 pt 4 lc 2, \
    "dsmc-5.0.txt" using 1:(base3($1*k,$2,Pxy5)) notitle w lp lw 2 dt 4 pt 4 lc 2

#    @dummy, @dummy, @dummy, @dummy, @dummy, \
#    base(filter(k*x, 0.1),fluid(k*x)) title "Navier--Stokes" dt 1 lc 0 lw 2 , \
