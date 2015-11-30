#!/usr/bin/env gnuplot

set term epslatex standalone size 3.2, 2.4 font 8 mono # color dashed
set out "<name>.tex"
set key center bottom maxrows 6 width -10

set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.16
set ylabel '$\displaystyle\int_0^\frac12\!\! <expr> \mathrm{d}x$' \
    offset graph <xcoord>, 0.45 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set nokey

set xrange [0:0.135]
#set xrange [1e-3:0.15]
set yrange [<ymin>:<ymax>]

base(x) = x
asym(kn, y, p0) = <asym>

set macros
dummy = "NaN title ' ' lt -3"

plot \
    "heat.txt" using 1:(base($<column>)) title ( <is_heat> ? "Heat-conduction equation" : "" ) w l lw 3 lt 4 lc 5, \
    "asym.txt" using 1:(base($<column>)) title "SNIF equations with $T_{B1}\\ne0$" w l lw 3 lt 1 lc 3, \
    "snif.txt" using 1:(base($<column>)) title "SNIF equations with $T_{B1}=0$" w l lw 3 lt 2 lc 2, \
    "data.txt" using 1:(base($<column>*<corr>)) title "Uniform velocity grid" w lp lw 3 lt 2 lc 1 pt 6, \
    "high.txt" using 1:(base($<column>*<corr>)) title "Nonuniform velocity grid" w lp lw 3 lt 3 lc 4 pt 4


# heat lt 5 --> lt 4
# asym lw 3 --> lw 2
# snif lt 3 --> lt 2
# data lt 3 --> lt 2
# high lt 4 --> lt 5

