#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 9,
          'text.fontsize': 9,
          'legend.fontsize': 10,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [2.5,2.5]}
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

py.savefig('Y1.pdf', bbox_inches='tight')
