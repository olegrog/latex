#!/usr/bin/env python

import numpy as np
from scipy import interpolate

# -- Universal constants
R = 8.3145
T0 = 273.15

# -- External parameters
p0 = 101325

# -- Argon
MG = 39.95e-3
rhoG = lambda T: p0*MG/R/T
# from Kestin et al, 1984
muG = interpolate.interp1d(np.array([1000, 1500, 2000, 2500]) + T0,
    np.array([65.39, 81.33, 95.84, 109.45])*1e-6, kind = 'cubic')
kG = interpolate.interp1d(np.array([1000, 1500, 2000, 2500]) + T0,
    np.array([51.23, 63.73, 75.08, 85.73])*1e-3, kind = 'cubic')

# -- SS316L
ML = 55.95e-3
TM = 1700
TB = 3090
# from Kim, 1975
rhoL = lambda T: (7.4327 + 3.9338e-5*T - 1.8007e-10)*1e+3
muL = lambda T: np.exp((2385.2/T - 0.5958)*np.log(10))*1e-3
# from Schmidt-Hohagen & Egry, 2006
gamma = lambda T: 1.697 - (T-1402-T0)*8.89e-5

print(f'rhoG = {rhoG(TM):.2} at T = {TM}')
print(f'muG = {muG(TM):.2} at T = {TM}')
print(f'rhoL = {rhoL(TM):.3} at T = {TM}')
print(f'muL = {muL(np.inf):.2e}*exp({TM*(np.log(muL(TM))-np.log(muL(np.inf))):.3}/T)')
print(f'muL = {muL(TM):.2e}, {muL(TB):.2e} at T = {TM}, {TB}')
print(f'gamma = {gamma(0):.3} + {gamma(1)-gamma(0):.3}*T')

print('-- Estimations')
muL_mean = muL((TM+TB)/2)
dgamma = np.abs(gamma(1)-gamma(0))
dT = TB - TM
L = 1e-4    # size of the melt pool
U = dgamma*dT/muL_mean
print(f'muL = {muL((TM+TB)/2):.2e} at T = {(TM+TB)/2}')
print(f'U = {U:.2}')
print(f'ReL(L) = {rhoL(TM)*U*L/muL_mean:.2}')
print(f'ReG(10*L) = {10*rhoG(TM)*U*L/muG(TM):.2}')
