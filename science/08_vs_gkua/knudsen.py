#!/usr/bin/env python

import numpy as np
import argparse

parser = argparse.ArgumentParser(description='calculator of nondimensional numbers')
parser.add_argument('-g', '--gas', default='Air', help='gas type')
args = parser.parse_args()

# constants
kB = 1.3806505e-23
NA = 6.0221415e+23
g1 = 1.25 # 1.270042427

class Ar:
    mu0 = 2.117e-5
    m = 66.3e-27
    mm = 39.948e-3
    omega = 0.81
    gamma = 5./3

class Air:
    mu0 = 1.719e-5
    m = 48.1e-27
    mm = 28.98e-3
    omega = 0.77
    gamma = 7./5

gas = eval(args.gas)
T0, P0 = 273.15, 101325
mu = lambda T: gas.mu0*(T/T0)**gas.omega
nu = lambda T: np.sqrt(2*kB*T*NA/gas.mm)

# parameters of problem
L, m = 0.178, 2.5
H = np.array([ 105.0011178e+3, 90.00174149e+3, 70.00208638e+3 ])
T = np.array([ 208.835, 186.867, 219.575 ])
p = np.array([ 0.0144771, 0.183593, 5.2180613 ])
Kn = np.array([ 1.8903, 1.3338e-1, 5.5130E-03 ])
Re = np.array([ 1.0489, 2.6862e+1, 8.2006E+02 ])
Ma = np.array([ 1.330, 2.403, 3.032 ])

print "mass of %s: " % args.gas, gas.m, gas.mm/NA

print "Knudsen number"
print " - given  :", Kn
print " - from_mu:", 4/g1/np.sqrt(np.pi)*mu(T)*nu(T)/p/L
print "Reynolds number"
print " - given  :", Re
print " - from_Re:", np.sqrt(2*gas.gamma)*Ma*p*L/nu(T)/mu(T)



