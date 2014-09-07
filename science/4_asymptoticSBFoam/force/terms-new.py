#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [3,3]}
py.rcParams.update(params)
import numpy as np
import sys

d = float(sys.argv[2])
data = {}
g3, g7 = 1.947906335, 1.758705

def read_data(name):
    x,y,z,s1,s2,s3 = py.loadtxt(sys.argv[1] + name + '.raw').T
    xnew, ynew = np.dstack((np.arctan2(x-d,y), s1))[0].T
    sort = np.argsort(xnew)
    data['x'], data[name] = xnew[sort], ynew[sort]

def plot(Fx, label, factor=1):
    x = data['x']
    py.plot(x, factor*np.poly1d(np.polyfit(x, Fx, 12))(x), '-', lw=1.5, label=label)
    py.plot(x, factor*Fx, 'k-', lw=0.1)

read_data('ShearStress1')
read_data('Pressure0')
read_data('Pressure7')
read_data('ShearStress3')

plot(data['Pressure0'], r'$\mathrm{term\;}p_2^\dag$')
plot((5*g7-g3)/g7*data['Pressure7'], r'$\mathrm{term\;}(\partial T_0/\partial x_k)^2$')
plot(data['ShearStress1'], r'$\mathrm{term\;}\partial u_{i1}/\partial x_j$')
plot(data['ShearStress3'] - g3/g7*data['Pressure7'], r'$\mathrm{term\;}\partial^2 T_0/\partial x_i\partial x_j$')

#read_data('Force')
#plot(data['Force'], r'$b$')

py.xlim(-np.pi/2, np.pi/2)
ax = py.axes()
ax.set_xticks([-np.pi/2, np.pi/2])
ax.set_xticklabels([r'$-\pi/2$', r'$\pi/2$'])
from matplotlib.ticker import MultipleLocator
ax.xaxis.set_minor_locator(MultipleLocator(np.pi/2))

if sys.argv[1] == 'cylinder':
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_minor_locator(MultipleLocator(50))
    lp = -5
if sys.argv[1] == 'sphere':
    ax.set_yticks([-300, 0, 200, 500])
    ax.yaxis.set_minor_locator(MultipleLocator(100))
    lp = -5
if sys.argv[1] == 'cylinder-inv':
    ax.yaxis.set_major_locator(MultipleLocator(40))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    lp = 0
if sys.argv[1] == 'sphere-inv':
    ax.yaxis.set_major_locator(MultipleLocator(40))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    lp = 0

py.xlabel(r'$\varphi$', labelpad=-5)
py.ylabel(r'$-p_0 F_{x2}$', y=.75, labelpad=lp, rotation=0)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=2)
ax.axhline(lw=.5, c='k', ls=':')

py.savefig('terms-' + sys.argv[1] + '.pdf', bbox_inches='tight')

