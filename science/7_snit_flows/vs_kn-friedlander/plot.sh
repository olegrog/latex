#!/usr/bin/env gnuplot

set term epslatex standalone size 3.2, 2.4 font 8 mono # color dashed
set out "delta_p.tex"
set key center bottom maxrows 6 width -10

set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.17
set ylabel '$\displaystyle\frac{\Delta{p}}{\bar{p}}$' offset graph 0.3, 0.48 rotate by 0
# \times10^{5}

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set nokey

factor = 1
set log xy
set xrange [0.5e-3:2e-1]
set yrange [1e-7:1e-2*factor]

batt = 10*4.7
k = sqrt(pi)/2
base(delta,p) = delta*factor/p
asym(kn, delta, p) = (k*kn)**2 * base(delta, p)
mod_Kn(Kn,p) = Kn/p
mtorr2Pa(p) = 0.133322*p

k_B = 1.3806505e-23
T = 80 # 80 -- 187
L = 18e-3
amu = 1.660538921e-27

#  = <helium>  <argon>
d  = '0.219e-9 0.368e-9'
mu = '1.870e-5 2.096e-5'
m  = '4.003    39.95'

p2Kn(p,type) = (k_B*T) / (sqrt(2)*pi*word(d,type)**2*p*L)
#p2Kn(p,type) = 5./8*sqrt(2*pi*k_B*T/word(m,type)/amu)*word(mu,type) / (p*L)

plot \
    "asym.txt" using (mod_Kn($1,$3)):(asym($1,$2,$3)) title "SNIF equations with $T_{B1}\\ne0$" w l lw 3 lt 1 lc 3, \
    "asym-cyl.txt" using (mod_Kn($1,$3)):(asym($1,$2,$3)) title "SNIF equations with $T_{B1}\\ne0$" w l lw 5 lt 1 lc 3, \
    "snif.txt" using (mod_Kn($1,$3)):(asym($1,$2,$3)) title "SNIF equations with $T_{B1}=0$"    w l lw 3 lt 2 lc 2, \
    "data.txt" using (mod_Kn($1,$3)):(base($2,$3))             notitle                                   w l lw 3 lt 2 lc 1, \
    "data.txt" using (mod_Kn($1,$3)):(base($2,$3)):(base($4,$3))  title "Boltzmann equation"                w yerrorbars lw 2 lt 2 lc 1 pt 6, \
    "argon.txt" using (p2Kn(mtorr2Pa($1),2)):(base($2/batt,$1)) title "Experiment (Ar)"          w p lw 3 pt 1 ps 1, \
    "helium.txt" using (p2Kn(mtorr2Pa($1),1)):(base($2/batt,$1)) title "Experiment (He)"         w p lw 3 pt 2 ps 1


