#!/usr/bin/env gnuplot

set term epslatex standalone size 4.5, 3 font 9 color dashed
set out "<name>.tex"
set key center top

set xlabel '$y$' offset graph 0.48, 0.15
set ylabel '$p_{xy} - p^{\mathrm{(ex)}}_{xy}$' offset graph 0.54, 0.46 rotate by 0
set colors classic

set border 3
set xtics nomirror
set ytics nomirror
set grid
#set xzeroaxis

set lmargin 7
set bmargin 2

set xrange [0:.5]
set yrange [-0.0004:0.0004]

our_width = 2
our_point = 6

U = 0.02
factor = sqrt(2)
exact = <exact>
func(x) = (x/U - 2*exact)/(2*exact)
func(x) = x/U - 2*exact

x_0 = <buffer>
if (<hybrid>) set arrow from x_0, graph 0 to x_0, graph 1 nohead dt 4 lc 0 lw 3

plot \
    "d3v96<suffix>.txt" using 1:(func($3)) title 'D3Q96' w lp lw our_width pt our_point, \
    "dvm<suffix>.txt" using 1:(func($3)) title 'DVM' w lp lw our_width pt our_point, \
    "hyb-d3q19<suffix>.txt" using 1:(func($3)) title 'hybrid D3Q19' w lp lw our_width pt our_point, \
    "hyb-d3v96<suffix>.txt" using 1:(func($3)) title 'hybrid D3Q96' w lp lw our_width pt our_point
