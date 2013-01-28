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

kappa = 2.129475
w=2
data = py.loadtxt("data.txt")
r,old,s1,new,s2 = data.T

py.plot(r, abs(new-kappa), color='b', linewidth=w, label=r'$\mathrm{new}$')
py.errorbar(r, abs(new-kappa), yerr=s2, color='b')
py.plot(r, abs(old-kappa), color='r', linewidth=w, label=r'$\mathrm{old}$')
py.errorbar(r, old-kappa, yerr=s1, color='r')

py.axhline(color='black')
py.xlabel(r'$N_R$')
py.ylabel(r'$\lambda-\lambda_\mathrm{ref}$')
py.xlim(9,62)
py.legend()
py.loglog()

py.savefig('conver_heat.pdf')
#py.show()
