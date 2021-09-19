#!/usr/bin/env python

import sys, os, time
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

### function: ipow, coeff
problems = {
    'A': (6, 2),
    'B': (6, 1),
    'B1': (6, 2),
    'B2': (6, 1),
    'B3': (6, 1),
    'B4': (6, 1),
    'T0_1': (6, 5),
    'T0_2': (8, 1),
    'T1_1': (6, .625),
    'T1_2': (8, .125),
    'T2_1': (6, .625),
    'T2_2': (8, .125),
    'TT12': (6, .625),
    'TT2': (8, .125),
    'Q2': (6, 1),
    'Q3': (8, 1./7),
    'QQ22': (6, 1),
    'QQ3': (8, 1./7),
}
result, delta = {}, 0.1

def find_value(x, y):
    return ((x[0]**4*y[0]*(x[1]**2 - x[2]**2) - x[1]**4*y[1]*(x[0]**2 - x[2]**2) + x[2]**4*y[2]*(x[0]**2 - x[1]**2))
        /(x[0]**4*x[1]**2 - x[0]**4*x[2]**2 - x[0]**2*x[1]**4 + x[0]**2*x[2]**4 + x[1]**4*x[2]**2 - x[1]**2*x[2]**4))

calc_I = lambda n, xi, f: 8./15/np.sqrt(np.pi)*np.trapz(f*xi**n*np.exp(-xi**2), xi)

#dirs = filter(lambda d: d.isdigit(), os.walk('.').next()[1])[:3]
dirs = sys.argv
dirs.remove(dirs[0])
cases = list(map(int, dirs))
if len(dirs) != 3:
    raise RuntimeError()

result_dir='infty'
if not os.path.exists(result_dir):
    print("Create directory %s" % result_dir)
    os.mkdir(result_dir)

for problem, (ipow, coeff) in problems.items():
    func_filename = os.path.join(result_dir, problem + '.txt')
    func_exists = os.path.exists(func_filename)
    gammas, Ys = [], []
    for i, source_dir in enumerate(dirs):
        xi, y = np.loadtxt(os.path.join(source_dir, problem + '.txt')).T
        gammas.append(coeff*calc_I(ipow, xi, y))
        if not func_exists:
            func = interp1d(xi, y, kind='cubic')
            X = np.arange(xi.min(), xi.max()+delta, delta)
            Ys.append(func(X))
    gamma = find_value(cases, gammas)
    result[problem] = gamma
    print('%g*I_%d(%s) = %.10f' % (coeff, ipow, problem, gamma))
    if not func_exists:
        Y = find_value(cases, Ys)
        print('%g*I_%d(%s) = %.10f ---' % (coeff, ipow, problem, coeff*calc_I(ipow, X, Y)))
        np.savetxt(func_filename, np.transpose((X, Y)), fmt='%.10f')

print('gamma_3 =', result['T0_1'] + result['T0_2'])
print('gamma_7 =', -(result['B2'] + result['B3']))
print('gamma_8 =', (result['Q2'] - result['QQ22']) + (result['Q3'] - result['QQ3']))
print('gamma_10 =', (result['T1_1'] + result['T2_1'] - 2*result['TT12']) + (result['T1_2'] + result['T2_2'] - 2*result['TT2']))

