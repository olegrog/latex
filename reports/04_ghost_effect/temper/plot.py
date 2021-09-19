#!/usr/bin/env python

import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.style.use('classic')
matplotlib.use('pgf')
params = {
    'axes.labelsize': 11,
    'font.size': 11,
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'pgf.rcfonts': False,
    'figure.figsize': [5,4]
}
plt.rcParams.update(params)

files=['1.txt','2.txt']
for f in files:
    data = np.loadtxt(f)
    x,y = data.T
    p, = plt.plot(x[1:], y[1:], 'o-', lw=2, clip_on=False, zorder=10)
    c = plt.getp(p,'color')
    plt.plot(x[1], y[1], 'D', clip_on=False, zorder=10, color=c)
    plt.plot(x[0], y[0], 's', clip_on=False, zorder=10, color=c)

bbox = dict(boxstyle="round", fc="0.9")
arrow = dict(arrowstyle="->")
plt.annotate(r'$$x=0, z=0.5$$', xy=(.029, .986), xytext=(.02, .95), bbox=bbox, arrowprops=arrow)
plt.annotate(r'$$x=0, z=0.1696$$', xy=(.032, .885), xytext=(.03, .92), bbox=bbox, arrowprops=arrow)
plt.xlabel(r'$\mathrm{Kn}$', labelpad=-1)
plt.ylabel(r'$T$', y=.8, labelpad=-3, rotation=0)
plt.xlim(0,.05)

plt.tight_layout()
plt.savefig(sys.argv[1])
