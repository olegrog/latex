#!/usr/bin/gnuplot

set term epslatex standalone size 5,3.75 font 11 color dashed
set out "qflow.tex"
set key center bottom

set xlabel '$\mathrm{Kn}$' offset graph 0.48,0.14
set ylabel '$\displaystyle\frac{\int_0^{\hat{L}/2}\hat{q}_y\:\mathrm{d}x}{\hat{U}\hat{L}/2}$' offset graph 0.32,0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 6
set bmargin 2
set logscale xy
#set grid

set xrange [0.005:180]
set yrange [2e-5:0.03]

eta = 0.562773
k0 = -1.2540
free = 1/sqrt(pi) 
error = 1.0069

plot \
	"data.txt" using 1:($3) title "Tchremissine solution" w lp lw 2 lt 1 pt 6, \
	"sone.txt" using ($1*2/sqrt(pi)):($4) title "LBE solution" w lp lw 2 lt 1 lc 3 pt 3
