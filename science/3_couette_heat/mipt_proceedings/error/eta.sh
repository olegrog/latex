#!/usr/bin/env gnuplot

k = -0.1
b = 0.1
FIT_LIMIT=1e-25

f(x) = b*x**k
fit f(x) "data.txt" using 1:3 via b,k
