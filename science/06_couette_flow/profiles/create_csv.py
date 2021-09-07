#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys, os
from scipy.interpolate import interp1d

_, macro = sys.argv

N, L = 51, 0.5
U = [ 1, 2, 5 ]
Kn = [ 0.1, 1, 10 ]
labels = {
    "P":        "P/U^2",
    "Pxx":      "(P_{xx}-P_{yy})/U^2",
    "Pxy":      "P_{xy}/U",
    "Pyy":      "P_{yy}/U^2",
    "Pzz":      "(P_{zz}-P_{yy})/U^2",
    "omega":    "omega/U^2",
    "qx":       "q_x/U^2",
    "qy":       "q_y/U^2",
    "tau":      "tau/U^2",
    "vx":       "v_x/U"
}

X = np.linspace(0, L, N)
data = pd.DataFrame()
data['y'] = X

for kn in Kn:
    table = np.loadtxt("profile_%s-%s.txt" % (macro, kn)).T
    for i in xrange(3):
        func = interp1d(table[0], table[i+2], kind='cubic')
        values = func(X)
        values[np.fabs(values) < 1e-10] = 0
        data["%s (U=%d, Kn=%g)" % (labels[macro], U[i], kn)] = values

data.to_csv(sys.stdout, index=False, float_format="%.3e")
