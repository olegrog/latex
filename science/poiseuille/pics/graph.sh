#!/usr/bin/gnuplot

set term epslatex standalone size 5,3.75 font 11 color dashed
set out "graph.tex"
set key center top

set xlabel '$\mathrm{Kn}$' offset graph 0.48,0.14
set ylabel '$Q$' offset graph 0.13,0.48 rotate by 0

set border 3
set xtics nomirror
set ytics nomirror

set lmargin 4
set bmargin 2
set logscale xy
set grid

set xrange [0.015:180]
set yrange [0.35:5.8]

eta = 0.562773
sigma = 1.01619
Q = 1./(24*eta)
non(x) = Q/x
first(x)  = Q*(1/x+6*sqrt(pi)*eta)
second(x) = x<3 ? Q*(1/x + 6*(1+x)) : 1/0
bgk(x)    = x<.09 ? Q/x + sigma/2 + (2*(sigma**2)-1)*eta*x + (6-9*sigma+4*sigma**3)/6*(2*eta*x)**2 : 1/0
cerci(x)  = x<3 ? Q*(1/x + 6*1.1209+12*0.2335*x) : 1/0
new(x)    = Q*(1/x + 6 + 12/(1+sqrt(3)/2/x))
free(x)   = x>3 ? .5/pi**.5*log(2/sqrt(pi)*x) : 1/0

plot \
	non(x) title "Navier--Stokes without slip" lt 2 lc 7 lw 2, \
	first(x) title "Maxwell slip (MM), 1879" lt 4 lc 7 lw 2, \
	"bgk.txt" using ($1/(2*eta)):2 title "Cercignani (BGK), 1963" w p lw 2 lt 1 lc 2 pt 2, \
	"bgk.txt" using ($1/1.5/(2*eta)):($2+Q/2/($1/(2*eta))) title "Cercignani (ES), 1966" w l lw 1 lt 5 lc 2, \
	"ohwada.txt" title "Ohwada et al. (HS), 1989" w lp lw 2 lt 1 lc 3 pt 3, \
	"itakura.txt" using ($1/(2*eta)):($2/2) title "Sone et al. (BGK), 1998" w l lw 1 lt 1 lc 5, \
	"smodel.txt" using ($1/(2*eta)):2 title "Sharipov (S), 1999" w p lw 2 lc 4 ps 2 pt 1, \
	"ewart.txt" title "Ewart et al. (exper), 2007" w p lw 10 lc 7 ps 0 pt 6, \
	cerci(x) title "Cercignani slip (HS), 2010" lt 2 lc 3 lw 2, \
	"our.txt" using ($1/$2):3 title "Tcheremissine (HS), 2011" w lp lw 2 lt 1 pt 6

#	bgk(x) notitle lc 2 lw 2 lt 1, \
#	free(x) title "free BGK", \
#	new(x) title "new slip", \
#	second(x) title "second order slip", \
#	"loyalka.txt" using ($1/(2*eta)):2 title "Loyalka (BGK), 1979" w lp lw 1, \
#	"cerci.txt" title "Cercignani et al. (var), 2010" w lp lw 1, \
