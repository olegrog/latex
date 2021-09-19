#!/usr/bin/env python

import numpy as np
import shapely.geometry as geom
import sys, os, vtk
from functools import partial
from collections import namedtuple
from scipy.interpolate import interp1d
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

_, infile, outfile, kn = sys.argv
k = np.sqrt(np.pi)*float(kn)/2
eta_max = 20

Curve = namedtuple('Curve', [ 'x', 'y', 'nx', 'ny', 't', 'tmin', 'tmax' ])
Fields = namedtuple('Fields', [ 'T0', 'dT0n', 'd2T0nt', 'dU1nt', 'dTt', 'dTn', 'curv', 'dUnt' ])

def get_point_data(out, name):
    vtk_array = out.GetPointData().GetArray(name)
    return vtk_to_numpy(vtk_array)

def get_points(out):
    vtk_array = out.GetPoints().GetData()
    numpy_array = vtk_to_numpy(vtk_array)
    return numpy_array[:,0], numpy_array[:,1], numpy_array[:,2]

def get_cell_data(out, name):
    vtk_array = out.GetCellData().GetArray(name)
    return vtk_to_numpy(vtk_array)

def get_cell_centers(out):
    centers = vtk.vtkCellCenters()
    centers.AddInputData(out)
    centers.Update()
    return get_points(centers.GetOutput())

def log_interp(eta, F):
#    func = interp1d(eta, np.log(F), kind='linear')
#    return lambda x: np.exp(func(x))
    return interp1d(eta, F, kind='cubic')

def add_array(data, F, name):
    vtk_array = numpy_to_vtk(F, array_type=vtk.VTK_FLOAT)
    vtk_array.SetName(name)
    data().AddArray(vtk_array)

def ellipse(a, b, sign):
    s = lambda t: np.sqrt(np.square(b*np.cos(t)) + np.square(a*np.sin(t)))
    return Curve(
        lambda t: a*np.cos(t),
        lambda t: b*np.sin(t),
        lambda t: sign*b/s(t)*np.cos(t),
        lambda t: sign*a/s(t)*np.sin(t),
        lambda x, y: np.arctan2(a*y, b*x),
        0, np.pi/2
    )

def norm_marg(f):
    f[0], f[-1] = 2*f[0], 2*f[-1]
    return f

def norm_marg2(f):
    f[0], f[-1] = f[1], f[-2]
    return f

class Boundary:
    def __init__(self, curve, X, Y, Z, fields):
        T = np.linspace(curve.tmin, curve.tmax, 1e3)
        self.__curve = curve
        self.__line = geom.LineString(list(zip(curve.x(T), curve.y(T))))
        mask = np.vectorize(self.__is_on_curve)(X, Y, Z)
        self.__create_interps(X, Y, fields, mask)

    def __is_on_curve(self, x, y, z):
        return self.__line.distance(geom.Point(x, y)) < 1e-4 and z == 0

    def __create_interps(self, X, Y, fields, m):
        print("Points on the boundary:", X[m].size)
        s = np.array(list(map(lambda x, y: self.project(geom.Point(x, y)), X[m], Y[m])))
        t = np.array(list(map(lambda x, y: self.tangent(geom.Point(x, y)), X[m], Y[m])))
        arg = np.argsort(s)
        kind = 'linear'
        # NB: we should use 'bounds_error=False' since argument can exceed s.max() by 1e-9 (looks like a bug)
        interp_wall = lambda f: interp1d(s[arg], norm_marg(f[arg]), kind, bounds_error=False, fill_value=2*f[arg][-1])
        interp_int = lambda f: interp1d(s[arg], norm_marg2(f[arg]), kind, bounds_error=False, fill_value=f[arg][-2])
        self.fields = Fields(
            interp_int(fields.T0[m]),
            interp_wall(fields.dT0n[m]),
            interp_wall(np.sum(fields.d2T0nt[m]*t, axis=1)),
            interp_wall(np.sum(fields.dU1nt[m]*t, axis=1)),
            interp_wall(np.sum(fields.dTt[m]*t, axis=1)),
            interp_wall(fields.dTn[m]),
            interp_wall(fields.curv[m]),
            interp_wall(np.sum(fields.dUnt[m]*t, axis=1)))
        '''
        import pylab as py
        S = np.array(map(lambda x, y: self.project(geom.Point(x, y)), X, Y))
        print S.min(), S.max()
        print s.min(), s.max()
        x = np.linspace(S.min(), S.max(), 1000)
        for name in self.fields._fields:
            print name, getattr(self.fields, name)(x).min(), getattr(self.fields, name)(x).max()
            py.plot(x, getattr(self.fields, name)(x), label=name)
            py.show()
        '''

    def normal(self, point):
        project = self.__line.interpolate(self.__line.project(point))
        t = self.__curve.t(project.x, project.y)
        return [ self.__curve.nx(t), self.__curve.ny(t), 0 ]

    def tangent(self, point):
        project = self.__line.interpolate(self.__line.project(point))
        t = self.__curve.t(project.x, project.y)
        return [ self.__curve.ny(t), -self.__curve.nx(t), 0 ]

    def project(self, point):
        return self.__line.project(point)

    def distance(self, point):
        return self.__line.distance(point)

def get_kn_layer_corr(boundary, x, y):
    point = geom.Point(x, y)
    eta = boundary.distance(point)/k
    if eta >= eta_max:
        return np.zeros(4)
    s = boundary.project(point)
    T0, dT0n, d2T0nt, dU1nt, dTt, dTn, curv, dUnt = boundary.fields
    T0, dT0n, d2T0nt, dU1nt, dTt, dTn, curv, dUnt = T0(s), dT0n(s), d2T0nt(s), dU1nt(s), dTt(s), dTn(s), curv(s), dUnt(s)
    rho1K = 1./T0 * dT0n * Omega_1(eta)
    T1K = T0/p0 * dTn * Theta_1(eta)
    T2K = (T0/p0)**2 * dT0n * (dT0n*Theta_3(eta)/T0 + (Theta_5(eta) - 2*Theta_3(eta))*curv)/2
    U2Kt = T0/p0*( np.sqrt(T0)*(d2T0nt*Y_a4(eta) - curv*dTt*(Y_a6(eta)-Y_a5(eta)/2)) + dUnt*Y_0(eta) )
    return rho1K, T1K, U2Kt, T2K

def add_real_data(x, y, T0, U1, T_old, U_old, data):
    rho, T, U = p0/T_old, np.copy(T_old), U_old*k/p0
    for boundary in boundaries:
        tang, norm = np.zeros((len(X), 3)), np.zeros((len(X), 3))
        for i in range(len(X)):
            tang[i,:] = boundary.tangent(geom.Point(X[i], Y[i]))
            norm[i,:] = boundary.normal(geom.Point(X[i], Y[i]))
        rho1K, T1K, U2Kt, T2K = np.transpose(list(map(partial(get_kn_layer_corr, boundary), x, y)))
        rho += k*rho1K
        T += k*T1K + k**2*T2K
        U += k**2*tang*np.tile(U2Kt, (3, 1)).T/p0
    names = [ 'rho', 'T', 'U', 'T0', 'U1', 'diff(T0)', 'diff(U1)' ]
    fields = [ rho, T, U, T0, U1/p0, T-T_old, U/k-U_old/p0 ]
    for name, field in zip(names, fields):
        add_array(data, field, name)
    # NB: we should return fields to prevent memory freeing,
    #     because VTK arrays refer to numpy arrays
    return fields


### Read the Knudsen-layer functions
pwd = os.path.dirname(os.path.realpath(__file__))
eta, Y_0, Y_1, Omega_1, Theta_1, Y_a4, Y_a5, Y_a6, Omega_3, Theta_3, Omega_5, Theta_5, intY_1 \
    = np.loadtxt(os.path.join(pwd, '../tables/kn-layer-hs.txt')).T
Y_0 = log_interp(eta, Y_0)
Y_1 = log_interp(eta, Y_1)
Y_a4 = log_interp(eta, Y_a4)
Y_a5 = log_interp(eta, Y_a5)        # NB: Y_a5 = -Y_a5(eta)
Y_a6 = log_interp(eta, Y_a6)
Theta_1 = log_interp(eta, Theta_1)  # NB: Theta_1 = -Theta_1(eta)
Omega_1 = log_interp(eta, Omega_1)
Theta_3 = log_interp(eta, Theta_3)  # NB: Theta_3 = -Theta_3(eta)
Omega_3 = log_interp(eta, Omega_3)
Theta_5 = log_interp(eta, Theta_5)  # NB: Theta_5 = -Theta_5(eta)
Omega_5 = log_interp(eta, Omega_5)
intY_1 = log_interp(eta, intY_1)    # NB: intY_1 = int_eta^infty Y_1(z)dz/2

### Read a VTK-file
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(infile)
reader.Update()
out = reader.GetOutput()

### Add points and cells to the new grid
grid = vtk.vtkUnstructuredGrid()
grid.SetPoints(out.GetPoints())
grid.SetCells(out.GetCellType(0), out.GetCells())

### Read point data
X, Y, Z = get_points(out)
T = get_point_data(out, 'T')
U = get_point_data(out, 'U')
T0 = get_point_data(out, 'T0')
U1 = get_point_data(out, 'U1')
p0 = get_point_data(out, 'p0')[0]
fields = Fields(T0,
    get_point_data(out, 'wallGradT0n'),
    get_point_data(out, 'wallGradGradT0nt'),
    get_point_data(out, 'wallGradU1nt'),
    get_point_data(out, 'wallGradTt'),
    get_point_data(out, 'wallGradTn'),
    get_point_data(out, 'curvature'),
    get_point_data(out, 'wallGradU1nt'),
)

### Create boundaries
inner = Boundary(ellipse(0.3, 0.7, 1), X, Y, Z, fields)
outer = Boundary(ellipse(1.5, 1, -1), X, Y, Z, fields)
boundaries = [ inner, outer ]

### Add point data
pd = add_real_data(X, Y, T0, U1, T, U, lambda: grid.GetPointData())

### Read cell data
X, Y, Z = get_cell_centers(out)
T0 = get_cell_data(out, 'T0')
U1 = get_cell_data(out, 'U1')
T = get_cell_data(out, 'T')
U = get_cell_data(out, 'U')

### Add cell data
cd = add_real_data(X, Y, T0, U1, T, U, lambda: grid.GetCellData())

### Write a VTK-file
writer = vtk.vtkUnstructuredGridWriter()
writer.SetInputData(grid)
writer.SetFileName(outfile)
writer.Update()

