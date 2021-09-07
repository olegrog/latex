#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 11,
          'text.fontsize': 11,
          'legend.fontsize': 12,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'figure.figsize': [3,3]}
py.rcParams.update(params)
import numpy as np

data = py.loadtxt('Y1_2.txt')
x,y = data.T
p, = py.plot(x, y, '-', lw=2)

py.xlabel(r'$\eta$', labelpad=-1)
py.semilogy()
py.axes().set_xticks([0,5,10,15])
py.axes().set_yticks([1,1e-2,1e-4])
py.ylabel(r'$\displaystyle\frac{Y_1}2$', y=.8, labelpad=-3, rotation=0)
py.xlim(0,15)

py.tight_layout()
py.savefig('Y1.pdf')
