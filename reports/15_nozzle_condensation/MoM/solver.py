#!/usr/bin/env python3

import sys, argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from scipy.integrate import solve_ivp
from numpy import pi, sqrt, log, exp, tanh

parser = argparse.ArgumentParser(description='Solver for the 1D nozzle using the MoM model')
parser.add_argument('-m', '--method', type=str, default='RK45', help='integration method')
parser.add_argument('-n', '--nozzle', type=str, default='conic', help='type of nozzle (conic|trapezoid)')
parser.add_argument('-d', '--diffuser', action='store_true', help='add a diffuser after the nozzle')
parser.add_argument('-P', '--pressure', type=float, default=0.6, help='initial pressure (atm)')
parser.add_argument('-T', '--temp', type=float, default=300, help='initial temperature (K)')
parser.add_argument('-a', '--alpha', type=float, default=0.6, help='condensation coefficient')
parser.add_argument('-R', '--radius', type=float, default=2.5e-3, help='nozzle throat radius (m)')
parser.add_argument('-H', '--height', type=float, default=5e-3, help='nozzle throat height (m)')
parser.add_argument('-W', '--width', type=float, default=0.4, help='nozzle throat width/height')
parser.add_argument('--xmin', type=float, default=0, help='initial distance from throat/length')
parser.add_argument('--Amax', type=float, default=1.25, help='max A/throat A')
parser.add_argument('--Amin', type=float, default=1.1, help='min A/throat A')
parser.add_argument('--phi', type=float, default=1.5, help='nozzle divergence angle (degree)')
parser.add_argument('--phi2', type=float, default=0.5, help='diffuser convergence angle (degree)')
parser.add_argument('-w0', type=float, default=0.1, help='initial vapor mass fraction')
parser.add_argument('-n0', type=float, default=1e-18, help='initial particle concentration (1/m)')
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

class vapor:        # Water
    gamma = 4/3
    M = 0.018       # kg/mol
    T_crit = 647    # K
    P_crit = 22e6   # Pa
    J0 = 1e32       # Hz/m^3
    Omega = 1.44
    sigma = lambda T: 1e-3*(93.6 + 9.13e-3*T - 2.75e-4*T**2)
    P_eq = lambda T: exp(77.34 - 7235/T - 8.2*log(T) + 5.71e-3*T)

class inert:        # Nitrogen
    gamma = 7/5
    M = 0.028       # kg/mol

class cond:
    c_p = 4200      # J/kg/K
    rho = 1000      # kg/m^3
    H = lambda T: fixed.R/vapor.M*(7235 - 8.2*T + 5.71e-3*T**2)

# Nozzle properties
Nozzle = namedtuple('Nozzle', [ 'L', 'A', 'dAdx'])
phi = args.phi*pi/180
phi2 = -args.phi2*pi/180
H, W = args.height, args.height*args.width
R = args.radius
profile0 = lambda x: phi*np.abs(x)*np.heaviside(L1-x, 1) + (phi2*(x-L1)+phi*L1)*np.heaviside(x-L1, 0)
profile1 = lambda x: phi*np.sign(x)*np.heaviside(L1-x, 1) + phi2*np.heaviside(x-L1, 0)
nozzle = {
    'conic': Nozzle(
        lambda A1, A2, phi: (sqrt(A2) - sqrt(A1))*R/phi,
        lambda x: pi*(R + profile0(x))**2,
        lambda x: 2*pi*(R + profile0(x))*profile1(x)),
    'trapezoid': Nozzle(
        lambda A1, A2, phi: (A2 - A1)*H/phi,
        lambda x: W*(H + profile0(x)),
        lambda x: W*profile1(x))
}[args.nozzle]
L1, L2 = nozzle.L(1, args.Amax, phi), nozzle.L(args.Amax, args.Amin, phi2)
L = L1 + args.diffuser*L2
xmin = args.xmin*L1 + 1e-8

# Material properties
_c_p_ = lambda mat: fixed.R*mat.gamma/(mat.gamma-1)/mat.M
_c_p = lambda g: (1-args.w0)*_c_p_(inert) + (args.w0-g)*_c_p_(vapor) + g*cond.c_p
_Mmean = lambda g: 1/((1-args.w0)/inert.M + (args.w0-g)/vapor.M)
_gamma = lambda g: 1/(1 - fixed.R/_c_p(g)/_Mmean(g))
_rho = lambda g, P, T: _gamma(g)*P/(_gamma(g)-1)/_c_p(g)/T
_Pvap = lambda g, P: (args.w0 - g)*P*_Mmean(g)/vapor.M
_kelvin = lambda T: 2*vapor.sigma(T)*vapor.M/cond.rho/fixed.R/T
_dotr = lambda Pvap, T, r: 5*pi/16*args.alpha/cond.rho*sqrt(vapor.M/2/pi/fixed.R/T)*(Pvap - vapor.P_eq(T)*exp(_kelvin(T)/r))
_J = np.vectorize(lambda T, S: 0. if S <= 1 else vapor.J0*exp(-16*pi/3*(vapor.Omega*(vapor.T_crit/T - 1))**3/log(S)**2))
_r_crit = np.vectorize(lambda T, S: np.infty if S <= 1 else _kelvin(T)/log(S))
_r_mean = np.vectorize(lambda y0, y2: np.nan if y0 <= 0 else sqrt(np.maximum(y2/y0, 0)))

def calc_all(t, y):
    A = nozzle.A(t)                             # nozzle area, m^2
    P, gMM, c_pT, y0, y1, y2, y3 = y            # pressure, Pa & others
    g = np.minimum(4*pi/3*cond.rho*y3, args.w0) # condensate mass fraction
    gamma = _gamma(g)                           # heat capacity ratio
    c_p = _c_p(g)                               # heat capacity
    Ma = sqrt(gMM/gamma)                        # Mach number
    T = c_pT/c_p                                # temperature, K
    Pvap = _Pvap(g, P)                          # vapor pressure, Pa
    rho = _rho(g, P, T)                         # density, kg/m^3
    r_mean = _r_mean(y0, y2)                    # mean particle radius, m
    S = Pvap/vapor.P_eq(T)                      # oversaturation
    J = _J(T, S)                                # nucleation rate, 1/m^3/s
    dotr = _dotr(Pvap, T, r_mean)               # growth rate, m/s
    r0 = _r_crit(T, S)                          # critical radius, m
    mu_k = np.array([y0, y1, y2, y3])*rho       # moments of PSD, m^{k-3}
    return gamma, Ma, T, P, Pvap, rho, g, S, r0, r_mean, J, dotr, mu_k

def stop(t, y):
    gamma, Ma, T, P, Pvap, rho, g, S, r0, r_mean, J, dotr, mu_k = calc_all(t, y)
    return S - args.maxS

def func(t, y):
    gamma, Ma, T, P, Pvap, rho, g, S, r0, r_mean, J, dotr, mu_k = calc_all(t, y)
    A = nozzle.A(t)                             # nozzle area, m^2
    dA = nozzle.dAdx(t)                         # nozzle expansion rate, 1/m
    rhoL = cond.rho                             # condensate density, kg/m^3
    H = cond.H(T)                               # latent heat, J/kg
    r0J1 = r0*J if r0 < np.infty else 0
    r0J2 = r0**2*J if r0 < np.infty else 0
    r0J3 = r0**3*J if r0 < np.infty else 0
    dotq = 4/3*pi*rhoL*(r0J3 + 3*dotr*mu_k[2])*H

    dydt = np.empty_like(y)
    dydt[3] = J*A/dotm                          # mu_0/rho
    dydt[4] = (r0J1 + dotr*mu_k[0])*A/dotm      # mu_1/rho
    dydt[5] = (r0J2 + 2*dotr*mu_k[1])*A/dotm    # mu_2/rho
    dydt[6] = (r0J3 + 3*dotr*mu_k[2])*A/dotm    # mu_3/rho

    dg = 4/3*pi*rhoL*dydt[6]
    Q = A*dotq/c_p/T/dotm
    dlngamma = (1-gamma)*((cond.c_p-_c_p_(vapor))/c_p + _Mmean(g)/vapor.M)*dg
    C = gamma/(gamma-1)

    dydt[0] = P*y[1]/(Ma**2-1)*(dlngamma/(gamma-1) - (3+2*y[1]/C)*dA/A + Q) # P
    dydt[1] = y[1]*dA/A - (1+y[1])*dydt[0]/P                                # gamma*M^2
    dydt[2] = y[2]/(1+y[1]/C/2)*(Q - dydt[1]/C/2 - dlngamma*Ma**2/2)        # c_p*T

    return dydt


if args.plots:
    T = np.linspace(fixed.Tmin, fixed.Tmax)
    X = np.linspace(xmin, L)
    G = np.linspace(0, args.w0)
    fig, axs = plt.subplots(2, 3, figsize=(12, 8))

    axs[0, 0].set_title('H(T), J/kg')
    axs[0, 0].plot(T, cond.H(T))

    axs[0, 1].set_title('Pressure(T), Pa')
    axs[0, 1].set_yscale('log')
    axs[0, 1].plot(T, vapor.P_eq(T), label='P_eq')
    axs[0, 1].plot(T, 0*T + args.pressure*fixed.P, label='P0')
    axs[0, 1].legend(loc='lower center')

    axs[0, 2].set_title('Molar mass (g), g/mol')
    axs[0, 2].plot(G, 1e3*_Mmean(G), label='mean')
    axs[0, 2].plot(G, 0*G + 1e3*vapor.M, '--', label='vapor')
    axs[0, 2].plot(G, 0*G + 1e3*inert.M, '--', label='inert')
    axs[0, 2].legend(loc='lower center')

    axs[1, 0].set_title('sigma(T), N/m')
    axs[1, 0].plot(T, vapor.sigma(T))

    axs[1, 1].set_title('gamma(g)')
    axs[1, 1].plot(G, _gamma(G))

    axs[1, 2].set_title('A(x)/A*')
    axs[1, 2].plot(X/L1, nozzle.A(X)/nozzle.A(0))

    fig.tight_layout()
    plt.show()

gamma, c_p = _gamma(0), _c_p(0)
A, Amax, Amin = nozzle.A(xmin), nozzle.A(L1), nozzle.A(L1+L2)
Ma = A/nozzle.A(0)
if nozzle.dAdx(xmin) < 0:
    Ma = 1/Ma
P0, T0 = args.pressure*fixed.P, args.temp
dotm = gamma*P0*Ma*A/sqrt((gamma-1)*c_p*T0)
U0 = Ma*sqrt((gamma-1)*c_p*T0)
S = _Pvap(0, P0)/vapor.P_eq(T0)
n0, r0 = np.maximum(args.n0/A, _J(T0, S)*sqrt(A)/U0), _r_crit(T0, S)
yk = n0*r0**np.arange(4)/_rho(0, P0, T0)
y0 = [P0, gamma*Ma**2, c_p*T0, *yk]

print(f'Initial values: T = {T0} K, P = {args.pressure} atm, Ma = {Ma:.3g}, w0 = {args.w0:.3g}')
print(f'Vapor mass flow rate (g/min) = {dotm*args.w0*1e3*60:.3g}')

if args.verbose:
    line = f'Geometry: A[m^2] = {A:.3g} -({L1:.3g} m, {args.phi:.2g}°)-> {Amax:.3g}'
    if args.diffuser:
        line += f' -({L2:.3g} m, {args.phi2:.2g}°)-> {Amin:.3g}'
    print(line)

    rho0 = _Mmean(0)*P0/fixed.R/T0
    v0, L = sqrt((_gamma(0)-1)*_c_p(0)*T0), sqrt(nozzle.A(0))
    n0 = 3*rho0*args.w0/4/pi/_kelvin(T0)**3/cond.rho
    nucl = vapor.J0*L/v0/n0
    grow = sqrt(vapor.M/_gamma(0)/_Mmean(0))*L*rho0/_kelvin(T0)/cond.rho

    print(f'Reference values:')
    print(f' -- rho0 = {_Mmean(0)*P0/fixed.R/T0:.3g} kg/m^3')
    print(f' -- t0 = {L/v0:.3g} s')
    print(f' -- n0 = {n0:.3g} 1/m^3')
    print(f' -- lambda_K = {_kelvin(T0):.3g} m')
    print(f'Dimensionless quantities:')
    print(f' -- T_crit/T_0 = {vapor.T_crit/T0:.3g}')
    print(f' -- latent heat/enthalpy = {cond.H(T0)/_c_p(0)/T0:.3g}')
    print(f' -- L/lambda_K = {L/_kelvin(T0):.3g}')
    print(f' -- rho_L/rho_0 = {cond.rho/rho0:.3g}')
    print(f' -- nucleation = {nucl:.3g}')
    print(f' -- growth = {grow:.3g}')

stop.terminal = True
# Without `first_step` constrain the first step can be too large when J is too small
sol = solve_ivp(func, [xmin, L], y0, atol=0, rtol=args.tol, events=stop,
    method=args.method, first_step=sqrt(A)/100)
print(f'Number of points = {sol.t.size}')
gamma, Ma, T, P, Pvap, rho, g, S, r0, r_mean, J, dotr, mu_k = calc_all(sol.t, sol.y)

Nucl = 4/3*pi*cond.rho*np.nan_to_num(r0, posinf=0)**3*J*A
Grow = 4*pi*cond.rho*dotr*mu_k[2]*A
X = sol.t/L1

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
axs[1, 1].plot(X, X*0 + 1, '--', label='saturation')
axs[1, 1].set_yscale('log')
axs[1, 1].legend(loc='upper center')

axs[1, 2].set_title('Nucleation, kg/m/s')
axs[1, 2].plot(X, Nucl)

axs[1, 3].set_title('Growth, kg/m/s')
axs[1, 3].plot(X, Grow)

#axs[1, 4].set_title('P, atm')
#axs[1, 4].plot(X, Pvap/fixed.P, label='vapor')
#axs[1, 4].plot(X, P_eq_hat/fixed.P, '--', label='equil')
#axs[1, 4].legend(loc='lower center')

axs[1, 4].set_title(f'A(x)/A*, phi={args.phi:.2g}°')
axs[1, 4].plot(X, nozzle.A(sol.t)/nozzle.A(0), '.-')

fig.tight_layout()
if args.pdf:
    plt.savefig('profiles.pdf')
else:
    plt.show()
