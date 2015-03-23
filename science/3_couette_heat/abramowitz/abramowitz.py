#!/usr/bin/python
# Abramowitz function of zero order

import numpy as np
import sys
from scipy.integrate import quad
from functools import partial

if len(sys.argv) > 1:
    power = int(sys.argv[1]) 
else:
    power = 0

X = np.logspace(-5, 1, 100)
f = lambda x,t: np.exp(-t*t - x/t) * t**power
Y = [ quad(partial(f, x), 0, np.infty)[0] for x in X ]

np.savetxt(sys.stdout, np.transpose((X, Y)), fmt='%1.4e') 
