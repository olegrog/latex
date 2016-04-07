#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys, os

_, path, part = sys.argv
eta_max = 5.4

labels = [{
    "B4":      "B^4",
    "QQ22":    "tilde{Q}_{22}^0",
    "QQ3":     "tilde{Q}_3^0",
    "Q2":      "Q_2^0",
    "Q3":      "Q_3^0"
},{
    "T1_1":     "T_1^1",
    "T1_2":     "T_2^1",
    "T2_1":     "T_1^2",
    "T2_2":     "T_2^2",
    "TT12":    "tilde{T}_{12}^0",
    "TT2":     "tilde{T}_2^0"
}]

data = pd.DataFrame()

for name, label in labels[int(part)].iteritems():
    eta, y = np.loadtxt(os.path.join(path, "%s.txt" % name)).T
    mask = eta <= eta_max
    data['eta'] = eta[mask]
    data[label] = y[mask]

data.to_csv(sys.stdout, index=False, float_format="%.3f")
