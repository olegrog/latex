#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set colors classic
set out "temp.tex"
set key center bottom maxrows 3 width -10

set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.14
set ylabel '$\displaystyle\int_0^\frac12 \frac{\tau\mathrm{d}y}{(\Delta v)^2} - \frac{\tau_{N\!S}^*}{(\Delta v)^2} - \left(\frac1{12} - \frac{\tau_{N\!S}^*}{(\Delta v)^2}\right)\frac{k}{1+k}$' \
    offset graph 1.76, 0.45 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy

set xrange [0.007:150]
set yrange [4e-4:4e-2]

k = sqrt(pi)/2
gamma1 = 1.270042427
gamma2 = 1.922284066
k0 = -1.2540
d1 = 2.4001
alpha = 1.25*gamma2/gamma1

free = 1./12
fluid = 1./24/alpha
tau1 = 2.21165e-02
tau2 = 4.47146e-02
tau5 = 1.16262e-01

#asym(x) = x < .1 ? fluid*(1+4*k0*x) + d1/4/alpha*x : 1/0
asym(x) = x < .1 ? (fluid + d1/4/alpha*x)/(1-2*k0*x)**2 : 1/0
f_asym(x) = x < 0.12 ? x : 1/0

base(x,y) = y - fluid - (free-fluid)*x/(1+x)
base2(x,y,v,f) = y/v - f/v - (free - f/v)*x/(1+x)
base3(x,y,v,f) = (y*v-0.5)/v/v - f/v - (free - f/v)*x/(1+x)

filter(x) = x < 0.1 ? x : 1/0

set macros
dummy = "NaN title ' ' lt -3"

plot \
    NaN title "Projection method"    w lp lw 2 dt 1 lc 1 pt 6, \
    NaN title "DSMC"             w lp lw 2 dt 1 lc 2 pt 4, \
    NaN title "Asymptotic"       w l lw 2 dt 1 lc 3, \
    NaN title "$\\Delta v=1$"    w l lw 2 dt 2 lc 0, \
    NaN title "$\\Delta v=2$"    w l lw 2 dt 3 lc 0, \
    NaN title "$\\Delta v=5$"    w l lw 2 dt 4 lc 0, \
    "asym-1.0.txt" using (filter($1)):(base2($1*k,$9,1,tau1)) notitle w l lw 2 dt 2 lc 3, \
    "asym-2.0.txt" using (filter($1)):(base2($1*k,$9,2,tau2)) notitle w l lw 2 dt 3 lc 3, \
    "asym-5.0.txt" using (filter($1)):(base2($1*k,$9,5,tau5)) notitle w l lw 2 dt 4 lc 3, \
    "data-1.0.txt" using 1:(base2($1*k,$9,1,tau1))   notitle w lp lw 2 dt 2 pt 6 lc 1, \
    "data-2.0.txt" using 1:(base2($1*k,$9,2,tau2))   notitle w lp lw 2 dt 3 pt 6 lc 1, \
    "data-5.0.txt" using 1:(base2($1*k,$9,5,tau5))   notitle w lp lw 2 dt 4 pt 6 lc 1, \
    "dsmc-1.0.txt" using 1:(base3($1*k,$9/2,1,tau1)) notitle w lp lw 2 dt 2 pt 4 lc 2, \
    "dsmc-2.0.txt" using 1:(base3($1*k,$9/2,2,tau2)) notitle w lp lw 2 dt 3 pt 4 lc 2, \
    "dsmc-5.0.txt" using 1:(base3($1*k,$9/2,5,tau5)) notitle w lp lw 2 dt 4 pt 4 lc 2

#    "asym-0.1.txt" using (filter($1)):(base2($1*k,$9,0.1,1.0004404781))    notitle w l lw 2 dt 6 lc 3, \
#    "data-0.1.txt" using 1:(base2($1*k,$9,0.1,1.0004404781))    notitle w lp lw 2 dt 6 pt 6 lc 1, \
#    "dsmc-0.1.txt" using 1:(base2($1*k,$9/2,0.1,1.0004404781))  notitle w lp lw 2 dt 6 pt 4 lc 2, \
#    base(k*x, asym(k*x)), \
#    (fluid-free)*x/(1+x) title "Navier--Stokes" dt 1 lc 0 lw 2 , \
#    "asym-hs.txt" using (f_asym($1)):(base($1*k,$4)) title "Asymptotic" w l lw 3 lt 7 lc 3, \

