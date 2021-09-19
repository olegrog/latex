#!/usr/bin/env python

import sys
from sympy import *

x, y, t = symbols('x y t')
A0, A1, A2, A3 = symbols('A0 A1 A2 A3')

def improved_series(i, K, var, n=8):
    return 4*y/(3+2*i)*series(K/(4*y/(3+2*i)), var, n=n).removeO()

def my_series(F, var, n=8):
    return series(F, var, n=n).removeO()

def calc_kernels(sign):
    R_s, R_d = y**2 + x**2, sign*(y**2 - x**2)
    r_s, r_d = y + x, sign*(y - x)

    C0 = sqrt(pi)*erf((r_s-r_d)/2)*exp((R_s-R_d)/2)
    C1 = (2+R_d)*C0 - 2*(r_s-r_d)
    C2 = 6*C1 + R_d**2*C0 - 2*(r_s**3-r_d**3)
    C3 = 10*C2 + R_d**2*C1 - 2*(r_s**5-r_d**5)
    C4 = 14*C3 + R_d**2*C2 - 2*(r_s**7-r_d**7)

    G0 = C0 / (x*y)
    G1 = (R_s*C0 - C1) / (2*(x*y)**2)
    G2 = (R_s**2*C0 - 2*R_s*C1 + C2) / (4*(x*y)**3)
    G3 = (R_s**3*C0 - 3*R_s**2*C1 + 3*R_s*C2 - C3) / (8*(x*y)**4)
    G4 = (R_s**4*C0 - 4*R_s**3*C1 + 6*R_s**2*C2 - 4*R_s*C3 + C4) / (16*(x*y)**5)

    J0 = integrate(t**2, (t, r_d, r_s)) / (x*y)
    J1 = integrate((R_s-t**2)*t**2, (t, r_d, r_s)) / (2*(x*y)**2)
    J2 = integrate((R_s-t**2)**2*t**2, (t, r_d, r_s)) / (4*(x*y)**3)
    J3 = integrate((R_s-t**2)**3*t**2, (t, r_d, r_s)) / (8*(x*y)**4)
    J4 = integrate((R_s-t**2)**4*t**2, (t, r_d, r_s)) / (16*(x*y)**5)

    K1 = y**3/x * (4*G1 - 2*J1)
    K2 = y**4/x**2 * (6*G2 - 2*G0 - 3*J2 + J0)
    K3 = y**5/x**3 * (10*G3 - 6*G1 - 5*J3 + 3*J1)
    K4 = y**6/x**4 * (70*G4 - 60*G2 + 6*G0 - 35*J4 + 30*J2 - 3*J0) / 4
    return [ improved_series(i,K,x) for i,K in enumerate([K1, K2, K3, K4]) ]

nu = exp(-x**2) + (2*x+1/x)*sqrt(pi)/2*erf(x)
print('nu =', my_series(nu, x))

K_plus, K_minus = calc_kernels(1), calc_kernels(-1)
print('K(y>x) = K_+')
for i, K in enumerate(K_plus):
    print('K%d =' % (i+1), K)
print()

print('K(y<x) = K_-')
for i, K in enumerate(K_minus):
    print('K%d =' % (i+1), improved_series(i, K.subs(x, y/t), y))
print()

# series expansion of A(x) and B(x) for small x
print('int_0^x (K_+-K_-) A exp(-y**2)')
for i, (K1, K2) in enumerate(zip(K_minus, K_plus)):
    integral = integrate((K2-improved_series(i,K1.subs(x, y/t),y)).subs(t, y/x)*my_series(exp(-y**2), y)*(A0+A1*y**2+A2*y**4+A3*y**6), y)
    print('int%d =' % (i+1), my_series(integral.subs(y,x), x))


