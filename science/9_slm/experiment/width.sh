#!/usr/bin/env gnuplot

set term epslatex standalone size 4, 3 font 9 color dashed
set out "width.tex"
set key top right

set xlabel '$U,\mathrm{mm/s}$' offset graph 0.0, 0.03
set ylabel '$d,\mathrm{\mu m}$'

file='trumpf.txt'

set log x
set xrange [90:3500]

lwe = 3
lws = 3
power = system("awk '!/^#/ { print $1 }' ".file." | sort -g | uniq")

plot for [p in power] \
    sprintf("<(grep '^%s' '%s')", p, file) u 2:3 w lp lw lwe pt 6 title p.' W (experiment)', \
    'heat_transfer-60W.txt'  u 1:2 w lp lw lws pt 5 dt 2 lt 1 title '60 W (simulation)', \
    'heat_transfer-90W.txt'  u 1:2 w lp lw lws pt 5 dt 2 lt 2 title '90 W (simulation)', \
    'heat_transfer-150W.txt' u 1:2 w lp lw lws pt 5 dt 2 lt 4 title '150 W (simulation)'

