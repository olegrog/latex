#!/usr/bin/env gnuplot

set term epslatex standalone size 4, 3 font 9 color dashed
set out "<name>.tex"
set key <key>

set xlabel '$<xlabel>$' offset graph 0.0, 0.03
set ylabel '$<ylabel>$'

exper='trumpf.txt'
numer='heat_transfer-'
crude=''
fine='-accurate'

set xzeroaxis
set log x
#set xrange [90:3500]

lwe = 3
lwn = 3
powers = system("awk '!/^#/ { print $1 }' ".exper." | sort -g | uniq")
base(power, speed) = <base>
func(width) = width - 55

plot for [p in powers] \
    sprintf("<(grep '^%s' '%s')", p, exper) u (base(p,$2)):(func($3)) w lp lw lwe pt 6 title p.' W (experiment)', \
    numer.'60W'.crude.'.txt'  u (base(60,$1)) :(func($2)) w lp lw lwn pt 5 dt 2 lt 1 title '60 W (simulation)', \
    numer.'90W'.crude.'.txt'  u (base(90,$1)) :(func($2)) w lp lw lwn pt 5 dt 2 lt 2 title '90 W (simulation)', \
    numer.'90W'.fine.'.txt'   u (base(90,$1)) :(func($2)) w lp lw lwn pt 3 dt 3 title '90 W (fine mesh)', \
    numer.'150W'.crude.'.txt' u (base(120,$1)):(func($2)) w lp lw lwn pt 5 dt 2 lt 4 title '150 W (simulation)'

