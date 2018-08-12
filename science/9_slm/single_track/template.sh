#!/usr/bin/env gnuplot

set term epslatex standalone size 4, 3 font 9 color dashed
set out "<name>.tex"
set key top left

set xlabel '$x$' offset graph 0.0, 0.03

set xzeroaxis
lw=3

plot \
    "<name>.txt" using 1:2 title "enthalpy"         w l lw lw, \
    "<name>.txt" using 1:4 title "temperature"      w l lw lw, \
    "<name>.txt" using 1:6 title "liquid fraction"  w l lw lw, \
    "<name>.txt" using 1:5 title "porosity"         w l lw lw, \
    "<name>.txt" using 1:3 title "diffusivity"      w l lw lw, \
