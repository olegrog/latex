#!/usr/bin/env python

import numpy as np
import sys
import StringIO
from scipy.interpolate import interp1d

_, U, sdata, sbase, sdiff = sys.argv

ncols = 9
columns = range(1, ncols)
ext = '-' + U + '.txt'

def interp(x, y):
    log_interp = interp1d(x, np.log(y), kind='cubic')
    return lambda x: np.exp(log_interp(x))

def diff(arr, base, column):
    mask = arr[0] < np.max(base[0])
    func = interp(base[0], np.abs(base[column]))
    return np.abs(np.abs(arr[column][mask]) - func(arr[0][mask]))

base = np.loadtxt(sbase + ext, usecols=range(ncols)).T
data = np.loadtxt(sdata + ext, usecols=range(ncols)).T
# some manual corrections
if sdata == 'dsmc':
    u = float(U)
    data[5:9] = (data[5:9]*u - 1)/u/2
    data[1] = data[1]/2

mask = data[0] < np.max(base[0])
Kn = data[0][mask]
result = [ diff(data, base, c) for c in columns ]

with open(sbase + ext, 'r') as f:
    first_line = f.readline()
names = first_line.split()[1:ncols+1]
sstream = StringIO.StringIO()
np.savetxt(sstream, [names], fmt='%10s')
header = sstream.getvalue().splitlines()[0]
sstream.close()

np.savetxt(sdiff + ext, np.vstack(([Kn], result)).T, fmt='%.4e', header=header[2:])

