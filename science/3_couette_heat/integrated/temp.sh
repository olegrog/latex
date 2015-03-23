#!/usr/bin/env gnuplot

set term epslatex standalone size 6, 4 font 11 color dashed
set out "temp.tex"
set key center bottom maxrows 3 width -3

set xlabel '$\mathrm{Kn}$' offset graph 0.48, 0.13
set ylabel '$\displaystyle\int_0^\frac12 \left[ \frac{\tau\mathrm{d}y}{(\Delta v)^2} - \frac{\tau_{N\!S}^*}{(\Delta v)^2} + \left(\frac{\tau_{N\!S}^*}{(\Delta v)^2}-\frac1{12}\right)\frac{k}{1+k} \right] \mathrm{d}y$' \
    offset graph 0.48, 0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy

set xrange [0.005:180]
set yrange [4e-4:5e-2]
#set yrange [0.02:0.04]

k = sqrt(pi)/2
gamma1 = 1.270042427
gamma2 = 1.922284066
k0 = -1.2540
d1 = 2.4001
alpha = 1.25*gamma2/gamma1

free = 1./12
fluid = 1./24/alpha
#asym(x) = x < .1 ? fluid*(1+4*k0*x) + d1/4/alpha*x : 1/0
asym(x) = x < .1 ? (fluid + d1/4/alpha*x)/(1-2*k0*x)**2 : 1/0
f_asym(x) = x < 0.12 ? x : 1/0

base(x,y) = y - fluid - (free-fluid)*x/(1+x)
#base(x,y) = y
#base2(x,y,v,f) = (y*v -.5)/v**2
#base2(x,y,v,f) = (y*v - .5)/v**2 - fluid - (free-fluid)*x/(1+x)
base2(x,y,v,f) = ( 2*y*v - f + (f-1)*x/(1+x) )/2/v**2 - free*x/(1+x)

filter(x) = x < 0.1 ? x : 1/0

set macros
dummy = "NaN title ' ' lt -3"

plot \
    NaN title "Tcheremissine"    w lp lw 2 lt 1 lc 1 pt 6, \
    NaN title "DSMC"             w lp lw 2 lt 1 lc 2 pt 4, \
    NaN title "Asymptotic"       w l lw 2 lt 1 lc 3, \
    NaN title "$\\Delta v=1$"    w l lw 2 lt 3 lc 0, \
    NaN title "$\\Delta v=2$"    w l lw 2 lt 4 lc 0, \
    NaN title "$\\Delta v=5$"    w l lw 2 lt 5 lc 0, \
    "asym-0.1.txt" using (filter($1)):(base2($1*k,$9,0.1,1.0004404781))    notitle w l lw 2 lt 2 lc 3, \
    "asym-1.0.txt" using (filter($1)):(base2($1*k,$9,1,1.04423265))        notitle w l lw 2 lt 3 lc 3, \
    "asym-2.0.txt" using (filter($1)):(base2($1*k,$9,2,1.1788573))         notitle w l lw 2 lt 4 lc 3, \
    "asym-5.0.txt" using (filter($1)):(base2($1*k,$9,5,2.1626155))         notitle w l lw 2 lt 5 lc 3, \
    "data-0.1.txt" using 1:(base2($1*k,$9,0.1,1.0004404781))    notitle w lp lw 2 lt 2 pt 6 lc 1, \
    "data-1.0.txt" using 1:(base2($1*k,$9,1,1.04423265))        notitle w lp lw 2 lt 3 pt 6 lc 1, \
    "data-2.0.txt" using 1:(base2($1*k,$9,2,1.1788573))         notitle w lp lw 2 lt 4 pt 6 lc 1, \
    "data-5.0.txt" using 1:(base2($1*k,$9,5,2.1626155))         notitle w lp lw 2 lt 5 pt 6 lc 1, \
    "dsmc-0.1.txt" using 1:(base2($1*k,$9/2,0.1,1.0004404781))  notitle w lp lw 2 lt 2 pt 4 lc 2, \
    "dsmc-1.0.txt" using 1:(base2($1*k,$9/2,1,1.04423265))      notitle w lp lw 2 lt 3 pt 4 lc 2, \
    "dsmc-2.0.txt.new" using 1:(base2($1*k,$9/2,2,1.1788573))     notitle w lp lw 2 lt 4 pt 6 lc 4, \
    "dsmc-2.0.txt" using 1:(base2($1*k,$9/2,2,1.1788573))       notitle w lp lw 2 lt 4 pt 4 lc 2, \
    "dsmc-5.0.txt" using 1:(base2($1*k,$9/2,5,2.1626155))       notitle w lp lw 2 lt 5 pt 4 lc 2

#    base(k*x, asym(k*x)), \
#    (fluid-free)*x/(1+x) title "Navier--Stokes" lt 1 lc 0 lw 2 , \
#    "asym-hs.txt" using (f_asym($1)):(base($1*k,$4)) title "Asymptotic" w l lw 3 lt 7 lc 3, \

