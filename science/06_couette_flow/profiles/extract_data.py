#!/usr/bin/env python
# Usage: ./extract_data.py /Users/olegrog/kesolver/data/ .

import sys, math, os, glob
import numpy as np
import pylab as py
from scipy.optimize import curve_fit
sys.path.append('/Users/olegrog/kesolver/tools')

import out2

_, indir, outdir = sys.argv

Kn = [ 0.1, 1, 10 ]
U = [ 0.1, 1.0, 2.0, 5.0 ]
L = 0.5
s2 = np.sqrt(2)
skip = 0.75         # skip this amount of timesteps from the beginning

linestyle = {
    0.1:    '-',
    1.0:    '--',
    2.0:    ':',
    5.0:    '-.'
}

linecolor = {
    0.1:    'red',
    1:      'green',
    10:     'blue'
}

show = False

result = {}
macros = ['Pxy', 'vx', 'qx', 'qy', 'Pxx', 'Pyy', 'Pzz', 'tau', 'P', 'omega']
for m in macros:
    result[m] = {}

def parabola_fit(x, a, b):
    return a*x**2 + b

def spline_fit(x, a, b, c):
    return a*x**2 + b*x + c

def average_data(path, nrows):
    timesteps = sorted([int(f[1:].split('.')[0]) for f in os.listdir(path) if f.endswith('.txt')])
    first, last = int(skip*len(timesteps)), len(timesteps)-1
    read_macros = lambda ts: np.array(out2.readMacros(os.path.join(path, 'm%s.txt' % str(ts)), nrows))
    data = read_macros(timesteps[-1])
    for ts in timesteps[first:last]:
        data += read_macros(ts)
    print('average with %d files' % (last - first + 1))
    return data / (last - first + 1)

for kn in Kn:
    for u in U:
        dirname = os.path.join(indir, 'couette' + "%.1f" % u, "%.3f" % kn)
        nodes, cells = out2.readNodesCells(dirname + '/couette.kei')
        Y = [0.5*kn*(nodes[cell.nodes[0]] + nodes[cell.nodes[2]])[1] for cell in cells ]
        data = average_data(os.path.join(dirname, 'result'), len(cells))
        data = (np.vstack([Y, data]).T)[np.array(Y).argsort()].T
        XX, rho, vx, vy, vz, T, Txx, Tyy, Tzz, qx, qy, qz, Pyz, Pxz, Pxy, H = data
        Pxx, Pyy, Pzz = rho*np.array([Txx, Tyy, Tzz])
        ind = len(XX)/2
        X = np.append(np.append(0, XX[ind:]), L)
        def symm(A):
            B = 0.5*((A+A[::-1])[ind:])
            X = XX[ind:]
            f0, fL = parabola_fit, spline_fit
            fit0, cov = curve_fit(f0, X[:2], B[:2])
            fitL, cov = curve_fit(fL, X[-3:], B[-3:])
            return np.append(np.append(f0(0, *fit0), B), fL(L, *fitL))
        def asym(A):
            B = 0.5*((A-A[::-1])[ind:])
            fL = spline_fit
            fitL, cov = curve_fit(fL, XX[-3:], B[-3:])
            return np.append(np.append(0, B), fL(L, *fitL))
        
        result['vx'][u] = asym(vx)/u/s2
        result['qx'][u] = asym(qx)/u/u/s2
        result['qy'][u] = asym(qy)/u/u/s2
        result['Pxy'][u] = symm(Pxy)/u
        result['Pxx'][u] = symm(Pxx-Pyy)/u/u
        result['Pyy'][u] = (symm(Pyy)-1)/u/u
        result['Pzz'][u] = symm(Pzz-Pyy)/u/u
        result['tau'][u] = (symm(T) - 1)/u/u
        result['P'][u] = ((symm(Pxx) + symm(Pyy) + symm(Pzz)) / 3 - 1)/u/u
        result['omega'][u] = (symm(rho) - 1)/u/u
        qy -= Pxy*vx
        if show and u>=0.1:
            Y = symm(Pxy)/u
            py.plot(X, Y, label=r'$kn=%.1f, u=%.1f' % (kn,u), linestyle=linestyle[u], color=linecolor[kn], lw=1.2) #, marker='o')
    for m in macros:
        print(m, np.transpose(list(result[m].values())).shape, np.array([X]).T.shape)
        #np.savetxt(outdir + '/profile_' + m + '-' + str(kn) + '.txt', np.hstack((np.transpose([X]), np.transpose(result[m].values()))),
        np.savetxt(outdir + '/profile_' + m + '-' + str(kn) + '.txt', np.vstack([X, result[m][0.1], result[m][1.0], result[m][2.0], result[m][5.0]]).T,
            fmt='%.5e', header='        x       U=0.1       U=1.0       U=2.0       U=5.0')

if show:
#    py.legend()
    py.show()

