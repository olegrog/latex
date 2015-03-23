#!/usr/bin/env gnuplot

set term epslatex standalone size 6, 4 font 11 color dashed
set out "shear.tex"
set key center bottom maxrows 5 width -4
set sample 101

set xlabel '$\mathrm{Kn}$' offset graph 0.48, 0.13
set ylabel '$\displaystyle\int_0^\frac12 \left(\frac{P_{xy}}{\Delta v} + \frac{\gamma_1^* k}{1+\gamma_1^*\sqrt\pi k}\right) \mathrm{d}y$' \
    offset graph 0.35, 0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale x
set xzeroaxis

set xrange [0.007:150]
set yrange [-0.015:0.022]

gamma1 = 1.270042427
k0_hs = -1.2540
k0_bkw = -1.01619
k = sqrt(pi)/2

free = -1/sqrt(pi) 
base(x,y) = .5*(gamma1*x/(1+sqrt(pi)*gamma1*x) + y)
base2(x,y,c) = .5*(c*gamma1*x/(1+sqrt(pi)*c*gamma1*x) + y)

asym_hs(x) = x < 3 ? -gamma1*x/(1-2*k0_hs*x) : 1/0
#asym_hs(x) = x < 3 ? -gamma1*x*(1+2*k0_hs*x) : 1/0
asym_bkw(x) = x < 3 ? -x/(1-2*k0_bkw*x) : 1/0
fluid(x)  = x < 0.07 ? -gamma1*x : 1/0
near_free(x) = free * x/(x+k)

filter(x) = x > 0.12 ? x : 1/0
f_asym(x) = x < 0.15 ? x : 1/0
filter2(x,c) = x < c ? x : 1/0

set macros
dummy = "NaN title ' ' lt -3"

plot \
    "bkw.txt" using (filter($1/k)/gamma1):(base(filter($1)/gamma1,$2)) title "BKW" w l lw 3 lt 1 lc 5, \
    base(k*x, asym_bkw(k*f_asym(x)*gamma1)) notitle w l lw 3 lt 1 lc 5, \
    "sone.txt" using ($1/k):(base($1,-$2/2)) title "Ohwada et al." w lp lw 3 lt 1 lc 0 pt 3 ps 1.5, \
    base(k*x,asym_hs(k*f_asym(x))) notitle lw 3 lt 1 lc 0, \
    NaN title "Tcheremissine"    w lp lw 2 lt 1 lc 1 pt 6, \
    NaN title "DSMC"             w lp lw 2 lt 1 lc 2 pt 4, \
    NaN title "Asymptotic"       w l lw 2 lt 1 lc 3, \
    NaN title "$\\Delta v\\to0$" w l lw 2 lt 1 lc 0, \
    NaN title "$\\Delta v=0.1$"  w l lw 2 lt 2 lc 0, \
    NaN title "$\\Delta v=1$"    w l lw 2 lt 3 lc 0, \
    NaN title "$\\Delta v=2$"    w l lw 2 lt 4 lc 0, \
    NaN title "$\\Delta v=5$"    w l lw 2 lt 5 lc 0, \
    @dummy, @dummy, @dummy, @dummy, @dummy, \
    "asym-0.1.txt" using (filter2($1,0.5)):(base2($1*k,$2,0.1000220/0.1))   notitle w l lw 2 lt 2 lc 3, \
    "asym-1.0.txt" using (filter2($1,0.2)):(base2($1*k,$2,1.0217401/1))     notitle w l lw 2 lt 3 lc 3, \
    "asym-2.0.txt" using (filter2($1,0.15)):(base2($1*k,$2,2.167795/2))      notitle w l lw 2 lt 4 lc 3, \
    "asym-5.0.txt" using (filter2($1,0.07)):(base2($1*k,$2,7.1917/5))        notitle w l lw 2 lt 5 lc 3, \
    "data-0.1.txt" using 1:(base2($1*k,$2,0.1000220/0.1))   notitle w lp lw 2 lt 2 pt 6 lc 1, \
    "data-1.0.txt" using 1:(base2($1*k,$2,1.0217401/1))     notitle w lp lw 2 lt 3 pt 6 lc 1, \
    "data-2.0.txt" using 1:(base2($1*k,$2,2.167795/2))      notitle w lp lw 2 lt 4 pt 6 lc 1, \
    "data-5.0.txt" using 1:(base2($1*k,$2,7.1917/5))        notitle w lp lw 2 lt 5 pt 6 lc 1, \
    "dsmc-0.1.txt" using 1:(base2($1*k,$2,0.1000220/0.1))   notitle w lp lw 2 lt 2 pt 4 lc 2, \
    "dsmc-1.0.txt" using 1:(base2($1*k,$2,1.0217401/1))     notitle w lp lw 2 lt 3 pt 4 lc 2, \
    "dsmc-2.0.txt.new" using 1:(base2($1*k,$2,2.167795/2))    notitle w lp lw 2 lt 4 pt 6 lc 4, \
    "dsmc-2.0.txt" using 1:(base2($1*k,$2,2.167795/2))      notitle w lp lw 2 lt 4 pt 4 lc 2, \
    "dsmc-5.0.txt" using 1:(base2($1*k,$2,7.1917/5))        notitle w lp lw 2 lt 5 pt 4 lc 2

#    base(k*x,fluid(k*x)) title "Navier--Stokes" lt 1 lc 0 lw 2 , \
