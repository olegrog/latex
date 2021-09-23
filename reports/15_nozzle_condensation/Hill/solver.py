#!/usr/bin/env python3
import sys, argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from scipy.integrate import solve_ivp
from numpy import pi, sqrt, log, exp, tanh

parser = argparse.ArgumentParser(description='Solver for the 1D nozzle using Hill1966 model')
parser.add_argument('-L', '--length', type=float, default=2e-3, help='nozzle length (m)')
parser.add_argument('-P', '--pressure', type=float, default=0.6, help='initial pressure (atm)')
parser.add_argument('-T', '--temp', type=float, default=300, help='initial temperature (K)')
parser.add_argument('-m', '--dotm', type=float, default=1.86e-4, help='vapor mass flow rate (kg/s)')
parser.add_argument('-a', '--alpha', type=float, default=0.6, help='condensation coefficient')
parser.add_argument('-A', '--Astar', type=float, default=0.06e-4, help='minimal cross-sectional area (m^2)')
parser.add_argument('--phi', type=float, default=pi/120, help='nozzle diverging angle')
parser.add_argument('-p', '--plots', action='store_true', help='draw additional plots')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
args = parser.parse_args()

### Constants
class fixed:
    P0 = 101325             # Pa
    R = 8.3145              # J/mol/K
    I0 = 1e32               # Hz/m^3
    Tmin = 50               # K
    Tmax = 500              # K
    xmin = -0.1*args.length # m
    xmax = args.length      # m

Material = namedtuple('Material', 'M c_p gamma T_crit P_crit Omega')
water = Material(0.018, 1865, 4/3, 647.14, 22.12e6, 1.44)

# Nozzle properties
_A = lambda x: args.Astar + args.phi*np.abs(x)
_dAdt = lambda x: args.phi*np.sign(x)

# Material properties
_c = lambda T: sqrt(water.gamma*fixed.R*T/water.M)
_rho = lambda P, T: P/fixed.R/T
_H = lambda T: fixed.R/water.M*(7235 - 8.2*T + 5.71e-3*T**2)
_rhoL = lambda T: 1e3*(0.08*tanh((T-225)/46.2) + 0.7415*(1-T/water.T_crit)**0.33 + 0.32)
_sigma = lambda T: 1e-3*(93.6 + 9.13e-3*T - 2.75e-4*T**2)
_p_eq = lambda T: exp(77.34 - 7235/T - 8.2*log(T) + 5.71e-3*T)
_p_eq_hat = lambda T, r: _p_eq(T)*exp(2*_sigma(T)*water.M/_rhoL(T)/fixed.R/T/r)
_dotr = lambda P, T, r: 5*pi/16*args.alpha/_rhoL(T)*sqrt(water.M/2/pi/fixed.R/T)*(P - _p_eq_hat(T, r))
_I = lambda T, S: 0 if S <= 0 else fixed.I0*exp(-16*pi/3*water.Omega**3*(water.T_crit/T - 1)**3/log(S)**2)
_r_crit = lambda T, S: 2*_sigma(T)*water.M/_rhoL(T)/fixed.R/T/log(S)

# Mixed properties
_dotm = lambda x, P, T, M, mu: _rho(P, T)*_A(x)*M*_c(T)/(1-mu)

if args.plots:
    T = np.linspace(fixed.Tmin, fixed.Tmax)
    X = np.linspace(fixed.xmin, fixed.xmax)
    fig, axs = plt.subplots(2, 2, figsize=(15, 8))
    axs[0, 0].plot(T, _H(T))
    axs[0, 0].set_title('H(T), J/kg')
    axs[0, 1].plot(T, _rhoL(T))
    axs[0, 1].set_title('rhoL(T), kg/m^3')
    axs[1, 0].plot(T, _sigma(T))
    axs[1, 0].set_title('sigma(T), N/m')
    axs[1, 1].plot(X, np.vectorize(_A)(X)/args.Astar)
    axs[1, 1].set_title('A(x)/A*')
    plt.show()

def func(t, y):
    P, T, Ma, mu, Q1, Q2, Q3 = y
    g = water.gamma                         # heat capacity ratio
    c = _c(T)                               # speed of sound, m/s
    u = c*Ma                                # flow speed, m/s
    rhoL = _rhoL(T)                         # condensate density, kg/m^3
    A = _A(t)                               # nozzle area, m^2
    dAdt = _dAdt(t)                         # nozzle expansion rate, 1/m
    r_mean = sqrt(2*Q1/Q3)                  # mean particle radius, m
    S = P/_p_eq(T)                          # saturation
    I = _I(T, S)                            # nucleation rate, 1/m^3/s
    dotr = _dotr(P, T, r_mean)              # growth rate, m/s
    l = _H(T)/water.c_p/T                   # lambda = latent heat/c_p/T
    r0 = _r_crit(T, S)                      # critical radius, m
    dotm = _dotm(t, P, T, Ma, mu)           # mass flow rate, kg/s

    Coord.append(t)
    Oversat.append(S)
    Nucl.append(4/3*pi*rhoL*r0**3*I*A)
    Grow.append(4*pi*rhoL*dotr*y[4])
    R_crit.append(r0)
    R_mean.append(r_mean)

    #print(t, y)
    dydt = np.empty_like(y)
    dydt[3] = rhoL/dotm*(Q1*dotr/u + 4*pi/3*I*A*r0**3)                      # liquid/vapor
    dydt[4] = Q2*dotr/u + 4*pi*I*A*r0**2                                    # Q1
    dydt[5] = Q3*dotr/u + 8*pi*I*A*r0                                       # Q2
    dydt[6] = 8*pi*I*A                                                      # Q3
    dydt[0] = P*((l-1)*dydt[3] - dAdt/A)/(1 - (1-mu)*((g-1)/g + 1/g/Ma**2)) # pressure
    dydt[1] = T*(l*dydt[3] + (g-1)/g*(1-mu)*dydt[0]/P)                      # temperature
    dydt[2] = Ma*(-(1-y[3])*dydt[0]/P/g/Ma**2 - dydt[1]/2/T)                # Mach

    #print(' -- ', dydt)
    return dydt

# P, T, M, mu, Q1, Q2, Q3
n0, r0 = 1e-10, 1e-8
A = args.Astar
M = _A(fixed.xmin)/A
P0 = args.pressure*fixed.P0
y0 = [P0, args.temp, M, 0, n0*r0**2*A, n0*r0*A, n0*A]
Coord = []; Oversat = []; Nucl = []; Grow = []; R_crit = []; R_mean = []
sol = solve_ivp(func, [fixed.xmin, 0.6*fixed.xmax], y0, method='BDF')

X = sol.t/args.length; Coord = np.array(Coord)/args.length
fig, axs = plt.subplots(2, 5, figsize=(15, 8))
axs[0, 0].plot(X, sol.y[0]/P0)
axs[0, 0].set_title('P/P0')
axs[0, 1].plot(X, sol.y[1])
axs[0, 1].set_title('T, K')
axs[0, 2].plot(X, sol.y[2])
axs[0, 2].set_title('Ma')
axs[0, 3].plot(X, sol.y[3])
axs[0, 3].set_title('liquid/vapor')
axs[0, 4].plot(X, sol.y[6])
axs[0, 4].set_title('N')
axs[1, 0].plot(Coord, R_mean, label='mean')
axs[1, 0].plot(Coord, R_crit, '--', label='crit')
axs[1, 0].set_yscale('log')
axs[1, 0].legend(loc='upper center')
axs[1, 0].set_title('Particle radius, m')
axs[1, 1].plot(Coord, Oversat)
axs[1, 1].set_title('Oversaturation')
axs[1, 2].plot(Coord, Nucl)
axs[1, 2].set_title('Nucleation')
axs[1, 3].plot(Coord, Grow)
axs[1, 3].set_title('Growth')
axs[1, 4].plot(X, np.vectorize(_A)(sol.t)/args.Astar)
axs[1, 4].set_title('A(x)/A*')
fig.tight_layout()
plt.show()
