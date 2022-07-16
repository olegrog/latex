#!/usr/bin/env python3

import sys, argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
from numpy import pi, sqrt, log, log10, exp, tanh

parser = argparse.ArgumentParser(description='Solver for the 1D nozzle using the MoM model')
parser.add_argument('-c', '--cond', type=str, default='H2O', help='condensable substance')
parser.add_argument('-m', '--method', type=str, default='RK45', help='integration method')
parser.add_argument('-n', '--nozzle', type=str, default='conic', help='type of nozzle (conic|trapezoid)')
parser.add_argument('-d', '--diffuser', action='store_true', help='add a diffuser after the nozzle')
parser.add_argument('-P', '--pressure', type=float, default=0.6, help='initial pressure (atm)')
parser.add_argument('-T', '--temp', type=float, default=300, help='initial temperature (K)')
parser.add_argument('-a', '--alpha', type=float, default=0.6, help='condensation coefficient')
parser.add_argument('-R', '--radius', type=float, default=2.5e-3, help='nozzle throat radius (m)')
parser.add_argument('-H', '--height', type=float, default=5e-3, help='nozzle throat height (m)')
parser.add_argument('-W', '--width', type=float, default=0.4, help='width/height at the nozzle throat')
parser.add_argument('--xmin', type=float, default=1e-6, help='initial distance from throat/length')
parser.add_argument('--Amax', type=float, default=1.25, help='max A/throat A')
parser.add_argument('--Amin', type=float, default=1.1, help='min A/throat A')
parser.add_argument('--phi', type=float, default=1.5, help='nozzle divergence angle (degree)')
parser.add_argument('--phi2', type=float, default=0.5, help='diffuser convergence angle (degree)')
parser.add_argument('-w0', type=float, default=0.1, help='initial vapor mass fraction')
parser.add_argument('-g0', type=float, default=1e-12, help='initial condensate mass fraction')
parser.add_argument('-r0', type=float, default=None, help='initial mean particle radius (m)')
parser.add_argument('-n0', type=float, default=None, help='initial particle concentration (1/m^3)')
parser.add_argument('-t', '--tol', type=float, default=1e-6, help='solution tolerance')
parser.add_argument('--Smax', type=float, default=1e8, help='maximum oversaturation')
parser.add_argument('--Nmin', type=int, default=50, help='minumum number of integration points')
parser.add_argument('--pdf', action='store_true', help='save PDF file instead')
parser.add_argument('--algebraic', action='store_true', help='find Ma using the algebraic equation')
parser.add_argument('--dry', action='store_true', help='plot solution without vapor as well')
parser.add_argument('-p', '--phase-diagram', action='store_true', help='draw trajectory on the phase diagram')
parser.add_argument('--plots', action='store_true', help='draw additional plots')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
args = parser.parse_args()

### Constants
class fixed:
    P = 101325              # Pa
    R = 8.3145              # J/mol/K
    N_A = 6.022e23          # 1/mol
    Tmin = 100              # K
    Tmax = 500              # K

class inert:        # Nitrogen
    gamma = 7/5
    M = 0.028       # kg/mol

### Condensable substance

class H2O:
    gamma = 1.327   # 4/3
    M = 0.018       # kg/mol
    T_crit = 647    # K
    P_crit = 22e6   # Pa
    T_trip = 216.5  # K
    J0 = 1e32       # Hz/m^3
    Omega = 1.44    # [Sinha, Wyslouzil, Wilemski 2009]
    sigma = lambda T: 1e-3*(93.6 + 9.13e-3*T - 2.75e-4*T**2) # [Sinha, Wyslouzil, Wilemski 2009]
    P_eq = lambda T: exp(77.34 - 7235/T - 8.2*log(T) + 5.71e-3*T)
    class cond:
        H = lambda T: fixed.R/vapor.M*(7235 - 8.2*T + 5.71e-3*T**2)
        c_p = 4200  # J/kg/K
        rho = 930   # kg/m^3 (at T < 200K)

class CO2:
    gamma = 1.289   # 7/5
    M = 0.044       # kg/mol
    T_crit = 304.26 # K
    P_crit = 7.38e6 # Pa
    T_trip = 216.55 # K
    J0 = 1e18       # Hz/m^3
    Omega = 1.70    # [Dingilian et al. 2020]
    sigma = lambda T: 1e-3*(0.05902*(304.26 - T)**1.25) # [Jasper 1972]
    #P_eq = lambda T: 101325*10**(-1353/T - 8.143*log10(T) + 0.006259*T + 24.619) # [Dingilian et al. 2020, for T>T_trip]
    P_eq = lambda T: 1e6*10**np.minimum(3.33084 - 692.853/(T - 25.007), 6.97568 - 1789.33/(T + 29.838)) # [Fernandez-Fassnacht & Del Rio 1984]
    class cond:
        H = lambda T: 1e3/vapor.M*(28.4363 - 0.021174*T + 1.19722e-4*T**2 - 4.64865e-7*T**3) # [Dingilian et al. 2020]
        c_p = 850   # J/kg/K
        rho = 1600  # kg/m^3 (dry ice)

vapor = globals()[args.cond]
cond = vapor.cond

# Nozzle properties
Nozzle = namedtuple('Nozzle', [ 'L', 'A', 'dAdx', 'var', 'value' ])
phi = args.phi*pi/180
phi2 = -args.phi2*pi/180
H_, W_ = args.height, args.height*args.width
R_ = args.radius
profile0 = lambda x: phi*np.abs(x)*np.heaviside(L1-x, 1) + (phi2*(x-L1)+phi*L1)*np.heaviside(x-L1, 0)
profile1 = lambda x: phi*np.sign(x)*np.heaviside(L1-x, 1) + phi2*np.heaviside(x-L1, 0)
nozzle = {
    'conic': Nozzle(
        lambda A1, A2, phi: (sqrt(A2) - sqrt(A1))*R_/phi,
        lambda x: pi*(R_ + profile0(x))**2,
        lambda x: 2*pi*(R_ + profile0(x))*profile1(x),
        'R', lambda x: R_ + profile0(x)
    ),
    'trapezoid': Nozzle(
        lambda A1, A2, phi: (A2 - A1)*H_/phi,
        lambda x: W_*(H_ + profile0(x)),
        lambda x: W_*profile1(x),
        'H', lambda x: H_ + profile0(x)
    )
}[args.nozzle]
L1, L2 = nozzle.L(1, args.Amax, phi), nozzle.L(args.Amax, args.Amin, phi2)
L = L1 + args.diffuser*L2
xmin = args.xmin*L1

# Material properties
_c_p_ = lambda mat: fixed.R*mat.gamma/(mat.gamma-1)/mat.M
_c_p = lambda g: (1-args.w0)*_c_p_(inert) + (args.w0-g)*_c_p_(vapor) + g*cond.c_p
_Mmean = lambda g: 1/((1-args.w0)/inert.M + (args.w0-g)/vapor.M)
_gamma = lambda g: 1/(1 - fixed.R/_c_p(g)/_Mmean(g))
_rho = lambda g, P, T: _gamma(g)*P/(_gamma(g)-1)/_c_p(g)/T
_rho_eq = lambda rho, rhoc, P, T: rho - P/fixed.R/T*_Mmean(rhoc/rho)
_Pvap = lambda g, P: (args.w0-g)*P*_Mmean(g)/vapor.M
_kelvin = np.vectorize(lambda T: 0 if T > vapor.T_crit else 2*vapor.sigma(T)*vapor.M/cond.rho/fixed.R/T)
_dotr = np.vectorize(lambda Pvap, T, r: 0 if T > vapor.T_crit else
    5*pi/16*args.alpha/cond.rho*sqrt(vapor.M/2/pi/fixed.R/T)*(Pvap - vapor.P_eq(T)*exp(_kelvin(T)/r)))
_J = np.vectorize(lambda T, S: 0. if S <= 1 or T > vapor.T_crit else
    vapor.J0*exp(-16*pi/3*(vapor.Omega*(vapor.T_crit/T - 1))**3/log(S)**2))
_r_crit = np.vectorize(lambda T, S: np.infty if S <= 1 else _kelvin(T)/log(S))
_r_mean = np.vectorize(lambda y0, y2: np.nan if y0 <= 0 else sqrt(np.maximum(y2/y0, 0)))

def calc_gas_dynamics(t, y, g):
    gamma = _gamma(g)                           # heat capacity ratio
    c_p = _c_p(g)                               # heat capacity, J/kg/K

    if args.algebraic:
        yc = (dotm/nozzle.A(t))**2
        ym, ye = y[:2]
        a, b, c = ye*yc - ym**2/2, 2*ye*yc - gamma*ym**2/(gamma-1), ye*yc
        Ma = sqrt(np.maximum((-b + sqrt(np.maximum(b**2 - 4*a*c, 0)))/2/a/gamma, 0))
        P = ym/(1 + gamma*Ma**2)
        T = (gamma*P*Ma)**2/(gamma-1)/yc/c_p
    else:
        P, gMM, T = y[:3]
        Ma = sqrt(np.maximum(gMM/gamma, 0))

    return gamma, Ma, T, P

def calc_dot_moments(mu, A, J, r_crit, dotr):
    dmu_rhodt = np.empty_like(mu)
    dmu_rhodt[0] = np.where(r_crit < np.infty, J*A/dotm, 0) # mu_0/rho
    for k in range(1, 4):
        with np.errstate(invalid='ignore'):
            r_critJ = np.where(r_crit < np.infty, r_crit**k*J, 0)
        dmu_rhodt[k] = (r_critJ + k*dotr*mu[k-1])*A/dotm    # mu/rho
    dg = 4/3*pi*cond.rho*dmu_rhodt[-1]
    return dg, dmu_rhodt

def calc_all(t, y):
    y0, y1, y2, y3 = y[-4:]
    g = 4*pi/3*cond.rho*y3                      # condensate mass fraction
    g = np.minimum(g, args.w0)
    gamma, Ma, T, P = calc_gas_dynamics(t, y, g)

    Pvap = _Pvap(g, P)                          # vapor pressure, Pa
    rho = _rho(g, P, T)                         # density, kg/m^3
    r_mean = _r_mean(y0, y2)                    # mean particle radius, m
    S = Pvap/vapor.P_eq(T)                      # oversaturation
    J = _J(T, S)                                # nucleation rate, 1/m^3/s
    dotr = _dotr(Pvap, T, r_mean)               # growth rate, m/s
    r_crit = _r_crit(T, S)                      # critical radius, m
    mu = np.array(y[-4:])*rho                   # moments of PSD, m^{k-3}

    return gamma, Ma, T, P, Pvap, rho, g, S, r_crit, r_mean, J, dotr, mu

def overSaturation(t, y):
    gamma, Ma, T, P, Pvap, rho, g, S, r_crit, r_mean, J, dotr, mu = calc_all(t, y)
    return S - args.Smax

def overHeating(t, y):
    gamma, Ma, T, P, Pvap, rho, g, S, r_crit, r_mean, J, dotr, mu = calc_all(t, y)
    return T - vapor.T_crit

def overSlowing(t, y):
    gamma, Ma, T, P, Pvap, rho, g, S, r_crit, r_mean, J, dotr, mu = calc_all(t, y)
    return Ma - 1

def overVaporization(t, y):
    gamma, Ma, T, P, Pvap, rho, g, S, r_crit, r_mean, J, dotr, mu = calc_all(t, y)
    return g

def func(t, y):
    if wet:
        gamma, Ma, T, P, Pvap, rho, g, S, r_crit, r_mean, J, dotr, mu = calc_all(t, y)
    else:
        g = g0
        gamma, Ma, T, P = calc_gas_dynamics(t, y, g)

    A = nozzle.A(t)                             # nozzle area, m^2
    dA = nozzle.dAdx(t)                         # nozzle expansion rate, m
    H = cond.H(T)                               # latent heat, J/kg
    dc_p = cond.c_p - _c_p_(vapor)              # dc_p/dg, J/kg/K

    dydt = np.empty_like(y)

    if wet:
        dg, dydt[-4:] = calc_dot_moments(mu, A, J, r_crit, dotr)
    else:
        dg = 0

    if args.algebraic:
        dydt[0] = -gamma*P*Ma**2*dA/A           # P*(1+gamma*Ma^2)
        dydt[1] = (H + dc_p*T)*dg               # c_p*T*(1+(gamma-1)*Ma^2/2)
    else:
        dM = _Mmean(g)/vapor.M*dg
        dQ = H/_c_p(g)/T*dg

        dydt[0] = P*y[1]/(Ma**2-1)*(-dA/A + dQ - dM)        # P
        dydt[1] = -y[1]*dA/A - (1+y[1])*dydt[0]/P           # gamma*M^2
        dydt[2] = T*(dA/A + (y[1]-1)/y[1]*dydt[0]/P + dM)   # T

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

### 1. Set the initial conditions

A0, Amax, Amin = nozzle.A(xmin), nozzle.A(L1), nozzle.A(L1+L2)
Ma0 = A0/nozzle.A(0)
if nozzle.dAdx(xmin) < 0:
    Ma0 = 1/Ma0
P0, T0 = args.pressure*fixed.P, args.temp
r0 = args.r0 if args.r0 else _kelvin(T0)
if args.n0:
    n0 = args.n0
    rhoc0 = n0*4/3*pi*r0**3*cond.rho
    rho0 = root_scalar(_rho_eq, args=(rhoc0, P0, T0), bracket=[_rho(0, P0, T0), _rho(args.w0, P0, T0)]).root
    g0 = rhoc0/rho0
else:
    g0 = args.g0
    rho0 = _rho(g0, P0, T0)
    rhoc0 = g0*rho0
    n0 = 3*rho0*g0/4/pi/cond.rho/r0**3

gamma0, c_p0 = _gamma(g0), _c_p(g0)
dotm = gamma0*P0*Ma0*A0/sqrt((gamma0-1)*c_p0*T0)
U0 = Ma0*sqrt((gamma0-1)*c_p0*T0)
yk = n0*r0**np.arange(4)/rho0

if args.algebraic:
    y0 = [P0*(1+gamma0*Ma0**2), c_p0*T0*(1+(gamma0-1)*Ma0**2/2), *yk]
else:
    y0 = [P0, gamma0*Ma0**2, T0, *yk]

### 2. Check the initial conditions

gamma, Ma, T, P, Pvap, rho, g, S, r_crit, r_mean, J, dotr, mu = calc_all(xmin, y0)
hatJ = J/vapor.J0
hatdotr = dotr*cond.rho/P*sqrt(fixed.R*T/vapor.M)
c0 = sqrt(gamma*fixed.R*T/_Mmean(g))
T_stag, P_stag = T0*(gamma+1)/2, P0*((gamma+1)/2)**(gamma/(gamma-1))
print('Initial values:')
print(f' -- T_o = {T_stag:.3g} K, P_o = {P_stag/fixed.P:.3g} atm (stagnation values)')
print(f' -- T_* = {T0:.3g} K, P_* = {args.pressure:.3g} atm, Ma_* = {Ma0:.3g}, w0 = {args.w0:.3g}')
print(f' -- S_* = {S:.3g}, c_* = {c0:.3g} m/s, g_* = {g:.3g}')
print(f' -- vapor mass flow rate = {dotm*args.w0*1e3*60:.3g} g/min')
print(f' -- condensate mass flow rate = {dotm*g*1e3*60:.3g} g/min')
print(f' -- particle number flow rate = {n0*c0*nozzle.A(xmin):.3g} per second')
if args.verbose:
    print(f' -- J_* = {J:.3g} m^-3 s^-1, hatJ_* = {hatJ:.3g}')
    print(f' -- dotr_* = {dotr:.3g} m/s, hatdotr_* = {hatdotr:.3g}')

if T != T0 or P != P0 or Ma != Ma0:
    raise RuntimeError(f'T = {T}, P = {P}, Ma = {Ma}')

if S > args.Smax:
    raise ValueError('Initial oversaturation is too high!')

if args.verbose:
    print(f'Geometry:')
    linePhi = f' -- phi = {args.phi:.2g}°'
    lineL = f' -- L[m] = {L1:.3g}'
    lineR = f' -- {nozzle.var}[m] = {nozzle.value(xmin):.3g} -> {nozzle.value(L1):.3g}'
    lineA = f' -- A[m^2] = {A0:.3g} -> {Amax:.3g}'
    linedA = f' -- dA/dx/A*[m^-1] = {nozzle.dAdx(xmin)/A0:.3g} -> {nozzle.dAdx(L1)/A0:.3g}'
    if args.diffuser:
        linePhi += f' -> {args.phi2:.2g}°'
        lineL += f' + {L2:.3g}'
        lineR += f' -> {nozzle.value(L):.3g}'
        lineA += f' -> {Amin:.3g}'
        linedA += f' -> {nozzle.dAdx(L)/A0:.3g}'
    print(linePhi); print(lineL); print(lineR); print(lineA); print(linedA);

    rho0 = _Mmean(0)*P0/fixed.R/T0
    v0 = sqrt(_gamma(0)*P0/rho0)
    L0 = sqrt(nozzle.A(0))
    n0 = 3*rho0/4/pi/_kelvin(T0)**3/cond.rho
    nucl = vapor.J0*L0/v0/n0 if n0 > 0 else 0
    grow = sqrt(vapor.M/_gamma(0)/_Mmean(0))*L0*rho0/_kelvin(T0)/cond.rho

    print(f'Reference values:')
    print(f' -- rho0 = {rho0:.3g} kg/m^3')
    print(f' -- t0 = {L0/v0:.3g} s')
    print(f' -- n0 = {n0:.3g} m^-3')
    print(f' -- lambda_K = {_kelvin(T0):.3g} m')
    print(f'Dimensionless quantities:')
    print(f' -- T_crit/T_0 = {vapor.T_crit/T0:.3g}')
    print(f' -- latent heat/enthalpy = {cond.H(T0)/_c_p(0)/T0:.3g}')
    print(f' -- L/lambda_K = {L0/_kelvin(T0):.3g}')
    print(f' -- rho_L/rho_0 = {cond.rho/rho0:.3g}')
    print(f' -- c_p(w0)/c_p(0) = {_c_p(args.w0)/_c_p(0):.3g}')
    print(f' -- nucleation = {nucl:.3g}')
    print(f' -- growth = {grow:.3g}')

### 3. Solve the initial-value problem

overSaturation.terminal = True
overHeating.terminal = True
overSlowing.terminal = True
overVaporization.terminal = True
overHeating.direction = 1
events = [ overSaturation, overHeating, overSlowing, overVaporization ]
wet = True
# Without `first_step` constrain the first step can be too large when J is too small
solver_kwargs = {
    'atol': 0,
    'rtol': args.tol,
    'method': args.method,
    'first_step': sqrt(A0)*args.tol,
    'max_step': L1/args.Nmin
}
sol = solve_ivp(func, [xmin, L], y0, events=events, **solver_kwargs)
gamma, Ma, T, P, Pvap, rho, g, S, r_crit, r_mean, J, dotr, mu = calc_all(sol.t, sol.y)

print(f'Number of points = {sol.t.size}')
if not sol.success:
    print('Terminated with the reason:', sol.message)
elif args.verbose:
    print(sol.message)

### 4. Plot the results

Nucl = 4/3*pi*cond.rho*np.nan_to_num(r_crit, posinf=0)**3*J*A0
Grow = 4*pi*cond.rho*dotr*mu[2]*A0
X = sol.t/L1

fig, axs = plt.subplots(2, 5, figsize=(15, 8))

axs[0, 0].set_title('P, atm')
axs[0, 0].plot(X, P/fixed.P)

axs[0, 1].set_title('T, K')
axs[0, 1].plot(X, T)
if np.min(T) < vapor.T_crit and vapor.T_crit < np.max(T):
    axs[0, 1].plot(X, T*0 + vapor.T_crit, ':', color='black', lw=1, label='critical')
    axs[0, 1].legend(loc='upper center')
if np.min(T) < vapor.T_trip and vapor.T_trip < np.max(T):
    axs[0, 1].plot(X, T*0 + vapor.T_trip, '--', color='black', lw=1, label='triple')
    axs[0, 1].legend(loc='upper center')
if args.verbose:
    h = _c_p(g)*T*(1 + (gamma-1)/2*Ma**2)
    ax = axs[0, 1].twinx()
    ax.plot(X, h, color='red', label='total enthalpy')
    ax.legend(loc='lower center')

axs[0, 2].set_title('Ma')
axs[0, 2].plot(X, Ma)

axs[0, 3].set_title('Condensate mass fraction')
axs[0, 3].plot(X, g)
axs[0, 3].plot(X, 0*X + args.w0, '--', label='maximum')
axs[0, 3].legend(loc='lower center')

axs[0, 4].set_title('Number of particles, 1/m')
axs[0, 4].plot(X, mu[0]*nozzle.A(sol.t))
axs[0, 4].set_yscale('log')

axs[1, 0].set_title('Mean particle radius, m')
axs[1, 0].plot(X, r_mean)
axs[1, 0].plot(X, r_crit, '--', label='critical')
if args.verbose:
    axs[1, 0].plot(X, _kelvin(T), ':', label='kelvin')
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

axs[1, 4].set_title(f'A(x)/A*, phi={args.phi:.2g}°')
axs[1, 4].plot(X, nozzle.A(sol.t)/nozzle.A(0), '.-')

Twet = np.copy(T)
Pwet = np.copy(Pvap)

if args.dry:
    wet = False
    print('Dry solution')

    sol = solve_ivp(func, [xmin, L], y0[:-4], **solver_kwargs)
    print(f'Number of points = {sol.t.size}')
    if not sol.success:
        print('Terminated with the reason:', sol.message)
    gamma, Ma, T, P = calc_gas_dynamics(sol.t, sol.y, g0)
    X = sol.t/L1

    opts = ['g:']
    axs[0, 0].plot(X, P/fixed.P, *opts)
    axs[0, 1].plot(X, T, *opts)
    axs[0, 2].plot(X, Ma, *opts)

    # Area Mach number relation for isentropic flow
    def isentropic(Ma, Aratio):
        gp, gm = gamma+1, gamma-1
        return Aratio**2 - (2/gp*(1+gm/2*Ma**2))**(gp/gm)/Ma**2

    gamma, c_p = _gamma(g0), _c_p(g0)
    Aratio = nozzle.A(X*L1)/nozzle.A(0)
    Ma = Ma0*np.array([ root_scalar(isentropic, args=x, x0=x, bracket=[1,1e2]).root for x in Aratio ])
    T = T0*(2 + (gamma-1)*Ma0**2)/(2 + (gamma-1)*Ma**2)
    P = dotm*sqrt((gamma-1)*c_p*T)/gamma/Ma/nozzle.A(X*L1)
    kwargs = { 'color': 'black', 'linewidth': 1, 'label': 'dry' }
    axs[0, 0].plot(X, P/fixed.P, **kwargs); axs[0, 0].legend(loc='upper center')
    axs[0, 1].plot(X, T, **kwargs); axs[0, 1].legend(loc='upper center')
    axs[0, 2].plot(X, Ma, **kwargs); axs[0, 2].legend(loc='upper center')

fig.tight_layout()
if args.pdf:
    plt.savefig('profiles.pdf')
else:
    plt.show()

if args.phase_diagram:
    plt.plot(Twet, Pwet, label='wet')
    Tb = np.linspace(vapor.T_trip, vapor.T_crit, args.Nmin)
    plt.plot(Tb, vapor.P_eq(Tb), label='saturation line')
    Tb = np.linspace(np.min(Twet), vapor.T_trip, args.Nmin)
    plt.plot(Tb, vapor.P_eq(Tb), label='sublimation line')
    kwargs = { 'color': 'black', 'marker': '.', 'linestyle': 'None', 'zorder': 10 }
    plt.plot([vapor.T_trip, vapor.T_crit], [vapor.P_eq(vapor.T_trip), vapor.P_eq(vapor.T_crit)], **kwargs)
    plt.semilogy()
    plt.legend()
    plt.show()
