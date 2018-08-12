#!/usr/bin/env gnuplot

set term epslatex standalone size 4, 3 font 9 color dashed
set out "width.tex"
set key top right

set xlabel '$U,\mathrm{mm/s}$' offset graph 0.0, 0.03
set ylabel '$d,\mathrm{\mu m}$'

file='single_track.txt'

set log x
set xrange [90:3500]

power = system("awk '!/^#/ { print $1 }' ".file." | sort -g | uniq")

plot for [p in power] \
    sprintf("<(grep '^%s' '%s')", p, file) u 2:3 w lp lw 2 pt 6 title p.' W', \
    'awesumm.txt' u 1:2 w lp lw 4 pt 5 lt 0 title 'simulation (90 W)'

