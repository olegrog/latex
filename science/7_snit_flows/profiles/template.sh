#!/usr/bin/env gnuplot

set term epslatex standalone size 3.2, 2.4 font 9 mono #color dashed #mono
set out '<name>.tex'
set key center bottom maxrows 6 width -10

set xlabel '$<xlabel>$' offset graph 0.48, 0.18
set ylabel '$\displaystyle <ylabel>$' offset graph <xcoord>, 0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin <lmargin>
set rmargin 1
set bmargin 2
#set key default
set nokey
set ytics <ytics>

tc = 3

stats '<patch>-heat.txt' using 1 name 'x' nooutput
stats '< paste <patch>-kes.txt <patch>-heat.txt' using ($<column>-column(<column>+tc)) name 'y' nooutput
stats '< paste <patch>-kgf.txt <patch>-heat.txt' using ($<column>-column(<column>+tc)) name 'kgf' nooutput
#stats '< paste <patch>-curv.txt <patch>-asym.txt' using ($<column>-column(<column>+tc)) name 'curv' nooutput

set xrange [x_min:x_max]
if (y_max*y_min <= 0) set xzeroaxis #lt 1
is_kgf(x) = kgf_max-kgf_min > 0 ? x : 1/0
#is_curv(x) = curv_max-curv_min > 1e-3 ? x : 1/0

plot \
    '< paste <patch>-kes.txt <patch>-heat.txt'  using 1:($<column>/<corr>-column(<column>+tc)) title 'Boltzmann equation' w l lw 2 lc 1, \
    '< paste <patch>-asym-first.txt <patch>-heat.txt' using 1:($<column>-column(<column>+tc)) title 'KGF with additional BC' w l lw 2 lc 2, \
    '< paste <patch>-asym-second.txt <patch>-heat.txt' using 1:($<column>-column(<column>+tc)) title 'KGF with additional BC' w l lw 2 lc 5, \
    '< paste <patch>-kgf.txt <patch>-heat.txt'  using (is_kgf($1)):($<column>-column(<column>+tc)) title 'KGF equations' w l lw 2 lt 3 lc 3

#    '< paste <patch>-curv.txt <patch>-heat.txt' using (is_curv($1)):($<column>-column(<column>+tc)) title 'KGF with curvature' w l lw 3 lt 5 lc 5
