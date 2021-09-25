#!/usr/bin/env python3

import sys, argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from scipy.integrate import solve_ivp
from numpy import pi, sqrt, log, exp, tanh

parser = argparse.ArgumentParser(description='Solver for the 1D nozzle using the MoM model')
parser.add_argument('-n', '--nozzle', type=str, default='conic', help='type of nozzle (conic|trapezoid)')
parser.add_argument('-P', '--pressure', type=float, default=0.6, help='initial pressure (atm)')
parser.add_argument('-T', '--temp', type=float, default=300, help='initial temperature (K)')
parser.add_argument('-a', '--alpha', type=float, default=0.6, help='condensation coefficient')
parser.add_argument('-R', '--radius', type=float, default=2.5e-3, help='nozzle throat radius (m)')
parser.add_argument('-H', '--height', type=float, default=5e-3, help='nozzle throat height (m)')
parser.add_argument('-W', '--width', type=float, default=0.4, help='nozzle throat width/height')
parser.add_argument('-A', '--Aratio', type=float, default=1.5, help='max A/throat A')
parser.add_argument('--xmin', type=float, default=0, help='initial distance from throat/length')
parser.add_argument('--phi', type=float, default=1.5, help='nozzle divergence angle (degree)')
parser.add_argument('-w0', type=float, default=0.1, help='initial vapor mass fraction')
parser.add_argument('-n0', type=float, default=1e-18, help='initial concentration (1/m)')
parser.add_argument('-r0', type=float, default=1e-8, help='initial particle radius (m)')
parser.add_argument('-t', '--tol', type=float, default=1e-6, help='solution tolerance')
parser.add_argument('--maxS', type=float, default=1e8, help='maximum oversaturation')
parser.add_argument('--pdf', action='store_true', help='save PDF file instead')
parser.add_argument('-p', '--plots', action='store_true', help='draw additional plots')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
args = parser.parse_args()

### Constants
class fixed:
    P = 101325              # Pa
    R = 8.3145              # J/mol/K
    N_A = 6.022e23          # 1/mol
    Tmin = 100              # K
    Tmax = 500              # K
    Psmall = P*np.finfo(float).eps
    Plarge = P/np.finfo(float).eps
    rlarge = 1e-4

class vapor:        # Water
    M = 0.018       # kg/mol
    c_p = 1865      # J/kg/K
    T_crit = 647    # K
    P_crit = 22e6   # Pa
    J0 = 1e32       # Hz/m^3
    Omega = 1.44
    sigma = lambda T: 1e-3*(93.6 + 9.13e-3*T - 2.75e-4*T**2)
    P_eq = lambda T: exp(77.34 - 7235/T - 8.2*log(T) + 5.71e-3*T) + fixed.Psmall


class inert:        # Nitrogen
    M = 0.028       # kg/mol
    c_p = 1039      # J/kg/K

class cond:
    c_p = 4200      # J/kg/K
    rho = 1000      # kg/m^3
    H = lambda T: fixed.R/vapor.M*(7235 - 8.2*T + 5.71e-3*T**2)

# Nozzle properties
Nozzle = namedtuple('Nozzle', [ 'L', 'A', 'dAdx'])
phi = args.phi*pi/180
H, W = args.height, args.height*args.width
R = args.radius
nozzle = {
    'conic': Nozzle(
        (sqrt(args.Aratio) - 1)*R/phi,
        lambda x: pi*(R + phi*np.abs(x))**2,
        lambda x: 2*pi*(R + phi*np.abs(x))*phi*np.sign(x)),
    'trapezoid': Nozzle(
        (args.Aratio - 1)*H/phi,
        lambda x: W*(H + phi*np.abs(x)),
        lambda x: W*phi*np.sign(x))
}[args.nozzle]

# Material properties
_g = lambda y3: np.minimum(4*pi/3*cond.rho*y3, args.w0)
_Mmean = lambda y3: 1/((1-args.w0)/inert.M + (args.w0-_g(y3))/vapor.M)
_c_p = lambda y3: (1-args.w0)*inert.c_p + (args.w0-_g(y3))*vapor.c_p + _g(y3)*cond.c_p
_rho = lambda g, P, T, y3: g*P/(g-1)/_c_p(y3)/T
_Pvap = lambda y3, P: (args.w0 - _g(y3))*P*_Mmean(y3)/vapor.M
_gamma = lambda y3: 1/(1 - fixed.R/_c_p(y3)/_Mmean(y3))
_P_eq_hat = lambda T, r: np.minimum(vapor.P_eq(T)*exp(2*vapor.sigma(T)*vapor.M/cond.rho/fixed.R/T/r), fixed.Plarge)
_dotr = lambda Pvap, T, r: 5*pi/16*args.alpha/cond.rho*sqrt(vapor.M/2/pi/fixed.R/T)*(Pvap - _P_eq_hat(T, r))
_J = np.vectorize(lambda T, S: 0. if S <= 1 else vapor.J0*exp(-16*pi/3*(vapor.Omega*(vapor.T_crit/T - 1))**3/log(S)**2))
_r_crit = np.vectorize(lambda T, S: np.infty if S <= 1 else 2*vapor.sigma(T)*vapor.M/cond.rho/fixed.R/T/log(S))
_r_mean = np.vectorize(lambda y0, y2: np.nan if y0 <= 0 else sqrt(np.maximum(y2/y0, 0)))

if args.plots:
    T = np.linspace(fixed.Tmin, fixed.Tmax)
    X = nozzle.L*np.linspace(args.xmin, 1)
    fig, axs = plt.subplots(2, 2, figsize=(15, 8))
    axs[0, 0].plot(T, cond.H(T))
    axs[0, 0].set_title('H(T), J/kg')
    axs[0, 1].set_yscale('log')
    axs[0, 1].plot(T, vapor.P_eq(T), label='P_eq')
    axs[0, 1].plot(T, 0*T+args.pressure*fixed.P, label='P0')
    axs[0, 1].legend(loc='lower center')
    axs[0, 1].set_title('Pressure(T), Pa')
    axs[1, 0].plot(T, vapor.sigma(T))
    axs[1, 0].set_title('sigma(T), N/m')
    axs[1, 1].plot(X/nozzle.L, np.vectorize(nozzle.A)(X)/nozzle.A(0))
    axs[1, 1].set_title('A(x)/A*')
    plt.show()

def calc_all(t, y):
    A = nozzle.A(t)                             # nozzle area, m^2
    dAdt = nozzle.dAdx(t)                       # nozzle expansion rate, 1/m
    yc = (dotm/A)**2
    ym, ye, y0, y1, y2, y3 = y
    g = _gamma(y3)                              # heat capacity ratio
    c_p = _c_p(y3)                              # heat capacity
    a, b, c = ye*yc - ym**2/2, 2*ye*yc - g*ym**2/(g-1), ye*yc
    sgn = np.sign(dAdt)
    Ma = sqrt(np.maximum((-b + sgn*sqrt(np.maximum(b**2 - 4*a*c, 0)))/2/a/g, 1e-10))  # Mach number
    P = ym/(1 + g*Ma**2)                        # pressure, Pa
    Pvap = _Pvap(y3, P)                         # vapor pressure, Pa
    T = (g*P*Ma)**2/(g-1)/yc/c_p                # temperature, K
    rho = _rho(g, P, T, y3)                     # density, kg/m^3
    r_mean = _r_mean (y0, y2)                   # mean particle radius, m
    S = Pvap/vapor.P_eq(T)                      # oversaturation
    J = _J(T, S)                                # nucleation rate, 1/m^3/s
    dotr = _dotr(Pvap, T, r_mean)               # growth rate, m/s
    r0 = _r_crit(T, S)                          # critical radius, m
    mu_k = np.array([y0, y1, y2, y3])*rho
    return g, Ma, T, P, Pvap, rho, _g(y3), S, r0, r_mean, J, dotr, mu_k

def stop(t, y):
    gamma, Ma, T, P, Pvap, rho, g, S, r0, r_mean, J, dotr, mu_k = calc_all(t, y)
    return S - args.maxS

def func(t, y):
    gamma, Ma, T, P, Pvap, rho, g, S, r0, r_mean, J, dotr, mu_k = calc_all(t, y)
    A = nozzle.A(t)                             # nozzle area, m^2
    dAdt = nozzle.dAdx(t)                       # nozzle expansion rate, 1/m
    rhoL = cond.rho                             # condensate density, kg/m^3
    H = cond.H(T)                               # latent heat, J/kg
    r0J1 = r0*J if r0 < np.infty else 0
    r0J2 = r0**2*J if r0 < np.infty else 0
    r0J3 = r0**3*J if r0 < np.infty else 0
    dotq = 4/3*pi*rhoL*(r0J3 + 3*dotr*mu_k[2])*H

    dydt = np.empty_like(y)
    dydt[0] = g*P*Ma**2*dAdt/A                  # P(1+g*Ma^2)
    dydt[1] = A*dotq/dotm                       # c_pT(1+(g-1)Ma^2/2)
    dydt[2] = J*A/dotm                          # mu_0/rho
    dydt[3] = (r0J1 + dotr*mu_k[0])*A/dotm      # mu_1/rho
    dydt[4] = (r0J2 + 2*dotr*mu_k[1])*A/dotm    # mu_2/rho
    dydt[5] = (r0J3 + 3*dotr*mu_k[2])*A/dotm    # mu_3/rho

    return dydt

g, c_p = _gamma(0), _c_p(0)
xmin, L = args.xmin*nozzle.L, nozzle.L
A, Amax = nozzle.A(xmin), nozzle.A(L)
Ma = A/nozzle.A(0)
if nozzle.dAdx(xmin) < 0:
    Ma = 1/Ma
P0, T = args.pressure*fixed.P, args.temp
dotm = g*P0*Ma*A/sqrt((g-1)*c_p*T)
yk = args.n0/A**2*args.r0**np.arange(4)/_rho(g, P0, T, 0)
y0 = [P0*(1+g*Ma**2), c_p*T*(1+(g-1)*Ma**2/2), yk[0], yk[1], yk[2], yk[3]]

print(f'Initial values: T = {T} K, P = {args.pressure} atm, Ma = {Ma:.3g}, w0 = {args.w0:.3g}')
print(f'Geometry: A = {A:.3g} -> {Amax:.3g} m^2, L = {L:.3g} m, diverging angle = {args.phi:.2g}°')
print(f'Vapor mass flow rate (g/min) = {dotm*args.w0*1e3*60:.3g}')

stop.terminal = True
sol = solve_ivp(func, [xmin, L], y0, atol=0, rtol=args.tol, events=stop, method='Radau')
gamma, Ma, T, P, Pvap, rho, g, S, r0, r_mean, J, dotr, mu_k = calc_all(sol.t, sol.y)

rhoL = cond.rho
Nucl = 4/3*pi*rhoL*r0**3*J*A
Grow = 4*pi*rhoL*dotr*mu_k[2]*A
P_eq_hat = _P_eq_hat(T, r_mean)
X = sol.t/L

fig, axs = plt.subplots(2, 5, figsize=(15, 8))

axs[0, 0].set_title('P, atm')
axs[0, 0].plot(X, P/fixed.P)

axs[0, 1].set_title('T, K')
axs[0, 1].plot(X, T)

axs[0, 2].set_title('Ma')
axs[0, 2].plot(X, Ma)

axs[0, 3].set_title('Condensate mass fraction')
axs[0, 3].plot(X, g)
axs[0, 3].plot(X, 0*X + args.w0, '--', label='maximum')
axs[0, 3].legend(loc='lower center')

axs[0, 4].set_title('Number of particles, 1/m')
axs[0, 4].plot(X, mu_k[0]*A)
axs[0, 4].set_yscale('log')

axs[1, 0].set_title('Mean particle radius, m')
axs[1, 0].plot(X, r_mean)
axs[1, 0].plot(X, r0, '--', label='critical')
axs[1, 0].set_yscale('log')
axs[1, 0].legend(loc='lower center')

axs[1, 1].set_title('Oversaturation')
axs[1, 1].plot(X, S)
axs[1, 1].plot(X, X*0 + 1, '--')
axs[1, 1].set_yscale('log')

axs[1, 2].set_title('Nucleation, kg/m/s')
axs[1, 2].plot(X, Nucl)

axs[1, 3].set_title('Growth, kg/m/s')
axs[1, 3].plot(X, Grow)

#axs[1, 4].set_title('P, atm')
#axs[1, 4].plot(X, Pvap/fixed.P, label='vapor')
#axs[1, 4].plot(X, P_eq_hat/fixed.P, '--', label='equil')
#axs[1, 4].legend(loc='lower center')

axs[1, 4].set_title(f'A(x)/A*, phi={args.phi:.2g}°')
axs[1, 4].plot(X, np.vectorize(nozzle.A)(sol.t)/nozzle.A(0), '.-')

fig.tight_layout()
if args.pdf:
    plt.savefig('profiles.pdf')
else:
    plt.show()
