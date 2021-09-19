#!/usr/bin/env gnuplot

set term epslatex standalone size 2.7, 2 font 8 color dashed header "\\usepackage{amsmath}\n"
set out "max.tex"
set key left bottom maxrows 3 width 0

set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.19
set ylabel '$\displaystyle\operatorname*{max}_{x_i,t}\frac{\sqrt{|\mathcal{M}|}}{|\mathcal{N}|}$' \
    offset graph 0.47, 0.45 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy

set xrange [0.007:150]
set yrange [1.2e-5:5e-3]

set macros
dummy = "NaN title ' ' lt -3"

plot \
    "exclusions-2e4.txt" u 1:(sqrt($3/20011)) title '$20011$'  w lp lw 3 lt 1 lc 1, \
    "exclusions-5e4.txt" u 1:(sqrt($3/50021)) title '$50021$'  w lp lw 3 lt 2 lc 2, \
    "exclusions-1e5.txt" u 1:(sqrt($3/100003)) title '$100003$' w lp lw 3 lt 3 lc 3, \
    "exclusions-2e5.txt" u 1:(sqrt($3/200003)) title '$200003$' w lp lw 3 lt 4 lc 4, \
    "exclusions-5e5.txt" u 1:(sqrt($3/500009)) title '$500009$' w lp lw 3 lt 5 lc 5
