#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'text.fontsize': 8,
          'legend.fontsize': 8,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'font.size': 8,
          'text.usetex': True,
          'figure.figsize': [3,3]}
py.rcParams.update(params)
import numpy as np

def plot(filename, factor, fmt, lw, label, xytext, va):
    d, Q, F = py.loadtxt(filename).T
    F = F * factor
    py.plot(d, F, fmt, lw=lw, label=label)
    zero = np.where(np.diff(np.sign(F)))[0][1]
    x1, x2, y1, y2 = d[zero], d[zero+1], F[zero], F[zero+1]
    d0 = x1 - y1*(x2-x1)/(y2-y1)
    py.annotate(
        r'$d='+ "{0:.2f}".format(d0) + r'$',
        xy=(d0,0), xytext=xytext,
        textcoords = 'offset points', ha = 'center', va = va,
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-.2'))

plot('cylinders-inv.txt', -2, '-', 1.5, r'$\mathrm{cylinders}\;\alpha=-5$', (-20, 20), 'bottom')
plot('spheres-inv.txt', -100, '-', 1.5, r'$\mathrm{spheres}\;\alpha=-5$', (-20, -20), 'top')

py.xlabel(r'$d$', labelpad=-5)
py.ylabel(r'$-\displaystyle\oint p_0 F_{x2} \mathrm{d}S$', labelpad=-5)
py.ylim(-50, 50)
ax = py.axes()
ax.set_xticks([0, 1])
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=0)
ax.axhline(lw=.5, c='k', ls=':')

from matplotlib.ticker import MultipleLocator
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(MultipleLocator(10))

py.savefig('forces-inv.pdf', bbox_inches='tight')

