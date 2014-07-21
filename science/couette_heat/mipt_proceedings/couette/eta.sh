#!/usr/bin/env gnuplot

eta = 0.562773
k0 = -1.2540
FIT_LIMIT=1e-25
error = 1.0069

q(x) = 2*eta*x/(1-sqrt(pi)*k0*x)
fit [0:0.2] q(x) "data.txt" using ($1):(-$4/1.0061) via eta, k0
