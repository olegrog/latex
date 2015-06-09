#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set out "alldiff.tex"
set key center bottom maxrows 6 width 0

set xlabel '$\mathrm{Kn}$' offset graph 0.5, 0.14

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy

set xrange [0.007:150]
set yrange [1e-6:1e-2]

set macros
dummy = "NaN title ' ' lt -3"

plot \
    NaN title "$P_{xy}$" w l lw 2 lt 1 lc 1, \
    NaN title "$v_x$"    w l lw 2 lt 1 lc 2, \
    NaN title "$q_x$"    w l lw 2 lt 1 lc 3, \
    NaN title "$\\tau$"  w l lw 2 lt 1 lc 4, \
    NaN title "$\\Delta v=0.1$" w l lw 2 lt 2 lc 0, \
    NaN title "$\\Delta v=1$"   w l lw 2 lt 3 lc 0, \
    NaN title "$\\Delta v=2$"   w l lw 2 lt 4 lc 0, \
    NaN title "$\\Delta v=5$"   w l lw 2 lt 5 lc 0, \
    "diff1-0.1.txt" using 1:2 notitle w l lw 2 lt 2 lc 1, \
    "diff1-1.0.txt" using 1:2 notitle w l lw 2 lt 3 lc 1, \
    "diff1-2.0.txt" using 1:2 notitle w l lw 2 lt 4 lc 1, \
    "diff1-5.0.txt" using 1:2 notitle w l lw 2 lt 5 lc 1, \
    "diff2-0.1.txt" using 1:2 notitle w lp lw 3 lt 2 lc 1 pt 6, \
    "diff2-1.0.txt" using 1:2 notitle w lp lw 3 lt 3 lc 1 pt 6, \
    "diff2-2.0.txt" using 1:2 notitle w lp lw 3 lt 4 lc 1 pt 6, \
    "diff2-5.0.txt" using 1:2 notitle w lp lw 3 lt 5 lc 1 pt 6, \
    "diff1-0.1.txt" using 1:3 notitle w l lw 2 lt 2 lc 2, \
    "diff1-1.0.txt" using 1:3 notitle w l lw 2 lt 3 lc 2, \
    "diff1-2.0.txt" using 1:3 notitle w l lw 2 lt 4 lc 2, \
    "diff1-5.0.txt" using 1:3 notitle w l lw 2 lt 5 lc 2, \
    "diff2-0.1.txt" using 1:3 notitle w lp lw 3 lt 2 lc 2 pt 6, \
    "diff2-1.0.txt" using 1:3 notitle w lp lw 3 lt 3 lc 2 pt 6, \
    "diff2-2.0.txt" using 1:3 notitle w lp lw 3 lt 4 lc 2 pt 6, \
    "diff2-5.0.txt" using 1:3 notitle w lp lw 3 lt 5 lc 2 pt 6, \
    "diff1-0.1.txt" using 1:4 notitle w l lw 2 lt 2 lc 3, \
    "diff1-1.0.txt" using 1:4 notitle w l lw 2 lt 3 lc 3, \
    "diff1-2.0.txt" using 1:4 notitle w l lw 2 lt 4 lc 3, \
    "diff1-5.0.txt" using 1:4 notitle w l lw 2 lt 5 lc 3, \
    "diff2-0.1.txt" using 1:4 notitle w lp lw 3 lt 2 lc 3 pt 6, \
    "diff2-1.0.txt" using 1:4 notitle w lp lw 3 lt 3 lc 3 pt 6, \
    "diff2-2.0.txt" using 1:4 notitle w lp lw 3 lt 4 lc 3 pt 6, \
    "diff2-5.0.txt" using 1:4 notitle w lp lw 3 lt 5 lc 3 pt 6, \
    "diff1-0.1.txt" using 1:9 notitle w l lw 2 lt 2 lc 4, \
    "diff1-1.0.txt" using 1:9 notitle w l lw 2 lt 3 lc 4, \
    "diff1-2.0.txt" using 1:9 notitle w l lw 2 lt 4 lc 4, \
    "diff1-5.0.txt" using 1:9 notitle w l lw 2 lt 5 lc 4, \
    "diff2-0.1.txt" using 1:9 notitle w lp lw 3 lt 2 lc 4 pt 6, \
    "diff2-1.0.txt" using 1:9 notitle w lp lw 3 lt 3 lc 4 pt 6, \
    "diff2-2.0.txt" using 1:9 notitle w lp lw 3 lt 4 lc 4 pt 6, \
    "diff2-5.0.txt" using 1:9 notitle w lp lw 3 lt 5 lc 4 pt 6

