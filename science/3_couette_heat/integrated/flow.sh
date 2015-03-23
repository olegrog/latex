#!/usr/bin/env gnuplot

set term epslatex standalone size 6, 4 font 11 color dashed
set out "flow.tex"
set key center bottom maxrows 6 width -4

set xlabel '$\mathrm{Kn}$' offset graph 0.48, 0.13
set ylabel '$\displaystyle\int_0^\frac12 \frac{v_x}{\Delta v} \mathrm{d}y - \frac{M_{N\!S}^*}{1-2k_0k}$' \
    offset graph 0.32, 0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy
set xzeroaxis

set xrange [0.007:150]
set yrange [2e-4:2e-1]
#set yrange [-.01:.06]

k = sqrt(pi)/2
k0 = -1.2540
gamma1 = 1.270042
fluid = 0.125

base(x,y) = y - fluid/(1-2*k0*x)
base2(x,y,c) = y - c/2/(1-2*k0*x)
#base(x,y) = y #- fluid/(1-2*k0*x)
#base2(x,y,c) = y #- fluid/(1-2*k0*x)

filter(x) = x > 0.12*k ? x : 1/0
f_asym(x) = x < 0.19 ? x : 1/0

set macros
dummy = "NaN title ' ' lt -3"

plot \
    @dummy, @dummy, @dummy, @dummy, @dummy, \
    "bkw.txt" using (filter($1/k)/gamma1):(base(filter($1)/gamma1,$3)) title "BKW" w l lw 3 lt 1 lc 5, \
    "asym-bkw.txt" using (f_asym($1)/gamma1):(base(f_asym($1*k)/gamma1,$2)) notitle w l lw 3 lt 1 lc 5, \
    "sone.txt" using ($1/k):(base($1,$3/2)) title "Ohwada et al." w lp lw 3 lt 1 lc 0 pt 3 ps 1.5, \
    "asym-hs.txt" using (f_asym($1)):(base(f_asym($1*k),$2)) notitle w l lw 3 lt 1 lc 0, \
    NaN title "Tcheremissine"    w lp lw 2 lt 1 lc 1 pt 6, \
    NaN title "DSMC"             w lp lw 2 lt 1 lc 2 pt 4, \
    NaN title "Asymptotic"       w l lw 2 lt 1 lc 3, \
    NaN title "$\\Delta v\\to0$" w l lw 2 lt 1 lc 0, \
    NaN title "$\\Delta v=0.1$"  w l lw 2 lt 2 lc 0, \
    NaN title "$\\Delta v=1$"    w l lw 2 lt 3 lc 0, \
    NaN title "$\\Delta v=2$"    w l lw 2 lt 4 lc 0, \
    NaN title "$\\Delta v=5$"    w l lw 2 lt 5 lc 0, \
    "asym-0.1.txt" using 1:(base2($1*k,$3,0.0249986/0.1))   notitle w l lw 2 lt 2 lc 3, \
    "asym-1.0.txt" using 1:(base2($1*k,$3,0.24867867/1))    notitle w l lw 2 lt 3 lc 3, \
    "asym-2.0.txt" using 1:(base2($1*k,$3,0.49055138/2))    notitle w l lw 2 lt 4 lc 3, \
    "asym-5.0.txt" using 1:(base2($1*k,$3,1.163882/5))      notitle w l lw 2 lt 5 lc 3, \
    "data-0.1.txt" using 1:(base2($1*k,$3,0.0249986/0.1))   notitle w lp lw 2 lt 2 pt 6 lc 1, \
    "data-1.0.txt" using 1:(base2($1*k,$3,0.24867867/1))    notitle w lp lw 2 lt 3 pt 6 lc 1, \
    "data-2.0.txt" using 1:(base2($1*k,$3,0.49055138/2))    notitle w lp lw 2 lt 4 pt 6 lc 1, \
    "data-5.0.txt" using 1:(base2($1*k,$3,1.163882/5))      notitle w lp lw 2 lt 5 pt 6 lc 1, \
    "dsmc-0.1.txt" using 1:(base2($1*k,$3,0.0249986/0.1))   notitle w lp lw 2 lt 2 pt 4 lc 2, \
    "dsmc-1.0.txt" using 1:(base2($1*k,$3,0.24867867/1))    notitle w lp lw 2 lt 3 pt 4 lc 2, \
    "dsmc-2.0.txt.new" using 1:(base2($1*k,$3,0.49055138/2))  notitle w lp lw 2 lt 4 pt 6 lc 4, \
    "dsmc-2.0.txt" using 1:(base2($1*k,$3,0.49055138/2))    notitle w lp lw 2 lt 4 pt 4 lc 2, \
    "dsmc-5.0.txt" using 1:(base2($1*k,$3,1.163882/5))      notitle w lp lw 2 lt 5 pt 4 lc 2
    
#    base(k*x,fluid) title "Navier--Stokes" lt 1 lc 0 lw 2 , \
