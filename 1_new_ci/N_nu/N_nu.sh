#!/usr/bin/env gnuplot

reset
n=50	#number of intervals
max=0.277	#max value
min=0.2745	#min value
width=(max-min)/n	#interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set term epslatex standalone size 4,3 font 11 color dashed
set output "N_nu.tex"
set xrange [min:max]
set yrange [0:]
#to put an empty boundary around the
#data inside an autoscaled graph.
#set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
set boxwidth width*0.9
set style fill solid 0.5	#fillstyle
set tics out nomirror
set xlabel '$N_\nu/N_\mathrm{kor}$'
set ylabel 'Frequency'
#count and plot
plot "N_nu.txt" u (hist($1/2000003,width)):(1.0) smooth freq w boxes lc rgb"green" notitle
