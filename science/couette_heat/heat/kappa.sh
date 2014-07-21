#!/usr/bin/env gnuplot

kappa = 2.12947
d = 2.4001
FIT_LIMIT=1e-25
error = 1.0069

q(x) = 2*kappa*x/(1+sqrt(pi)*d*x)
fit [0:0.2] q(x) "qflow.txt" using ($1/$2):3 via kappa, d
