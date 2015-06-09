#!/usr/bin/env gnuplot

set term epslatex standalone size 2.8, 2.1 font 8 color dashed
set out "<name>.tex"
set key <place> maxrows 4 width -1

set xlabel '$y$' offset graph 0.0, 0.03
set ylabel '$<ylabel>$' offset graph 0.05, 0.0

set xzeroaxis

set xrange [0:0.5]
set yrange [<ymin>:<ymax>]

base(y,v) = <base>
flag = <filter>
filter(x) = x/flag

set macros
dummy = "NaN title (flag ? ' ' : '') lt -3"

plot @dummy, \
    NaN title "$\\mathrm{Kn}=0.1$"   w l lw 2 lt 1 lc 1, \
    NaN title "$\\mathrm{Kn}=1$"     w l lw 2 lt 1 lc 2, \
    NaN title "$\\mathrm{Kn}=10$"    w l lw 2 lt 1 lc 3, \
    NaN title ( flag ? "$\\Delta v=0.1$" : "" ) w l lw 2 lt 2 lc 0, \
    NaN title "$\\Delta v=1$"       w l lw 2 lt 3 lc 0, \
    NaN title "$\\Delta v=2$"       w l lw 2 lt 4 lc 0, \
    NaN title "$\\Delta v=5$"       w l lw 2 lt 5 lc 0, \
    "profile_<macro>-0.1.txt" using (filter($1)):(base($2,0.1)) notitle w l lw 2 lt 2 lc 1, \
    "profile_<macro>-0.1.txt" using 1:(base($3,1.0)) notitle w l lw 2 lt 3 lc 1, \
    "profile_<macro>-0.1.txt" using 1:(base($4,2.0)) notitle w l lw 2 lt 4 lc 1, \
    "profile_<macro>-0.1.txt" using 1:(base($5,5.0)) notitle w l lw 2 lt 5 lc 1, \
    "profile_<macro>-1.txt"   using (filter($1)):(base($2,0.1)) notitle w l lw 2 lt 2 lc 2, \
    "profile_<macro>-1.txt"   using 1:(base($3,1.0)) notitle w l lw 2 lt 3 lc 2, \
    "profile_<macro>-1.txt"   using 1:(base($4,2.0)) notitle w l lw 2 lt 4 lc 2, \
    "profile_<macro>-1.txt"   using 1:(base($5,5.0)) notitle w l lw 2 lt 5 lc 2, \
    "profile_<macro>-10.txt"  using (filter($1)):(base($2,0.1)) notitle w l lw 2 lt 2 lc 3, \
    "profile_<macro>-10.txt"  using 1:(base($3,1.0)) notitle w l lw 2 lt 3 lc 3, \
    "profile_<macro>-10.txt"  using 1:(base($4,2.0)) notitle w l lw 2 lt 4 lc 3, \
    "profile_<macro>-10.txt"  using 1:(base($5,5.0)) notitle w l lw 2 lt 5 lc 3

