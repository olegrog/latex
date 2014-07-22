#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 11,
          'text.fontsize': 11,
          'legend.fontsize': 12,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'figure.figsize': [5,4]}
py.rcParams.update(params)
import numpy as np

files=['1.txt','2.txt']
for f in files:
	data = py.loadtxt(f)
	x,y = data.T
	p, = py.plot(x[1:], y[1:], 'o-', lw=2, clip_on=False, zorder=10)
	c = py.getp(p,'color')
	py.plot(x[1], y[1], 'D', clip_on=False, zorder=10, color=c)
	py.plot(x[0], y[0], 's', clip_on=False, zorder=10, color=c)

bbox = dict(boxstyle="round", fc="0.9")
arrow = dict(arrowstyle="->")
py.annotate(r'$$x=0, z=0.5$$', xy=(.029, .986), xytext=(.02, .95), bbox=bbox, arrowprops=arrow)
py.annotate(r'$$x=0, z=0.1696$$', xy=(.032, .885), xytext=(.03, .92), bbox=bbox, arrowprops=arrow)
py.xlabel(r'$\mathrm{Kn}$', labelpad=-1)
py.ylabel(r'$T$', y=.8, labelpad=-3, rotation=0)
py.xlim(0,.05)

py.tight_layout()
py.savefig('temper.pdf')
