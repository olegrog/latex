#!/usr/bin/env python

import numpy as np
import sys, os, vtk
from scipy.interpolate import griddata, interp1d
from vtk.util.numpy_support import vtk_to_numpy
from pyevtk.hl import pointsToVTK

_, infile, outfile, kn = sys.argv
k = np.sqrt(np.pi)*float(kn)/2

def get_point_data(reader, name):
    vtk_array = reader.GetOutput().GetPointData().GetArray(name)
    return vtk_to_numpy(vtk_array)

def get_points(reader):
    vtk_array = reader.GetOutput().GetPoints().GetData()
    numpy_array = vtk_to_numpy(vtk_array)
    return numpy_array[:,0], numpy_array[:,1], numpy_array[:,2]

def interp(eta, F):
    log_interp = interp1d(eta, np.log(F), kind='linear')
    eta_max = 20
    f = lambda x: np.where(x < eta_max, x, eta_max)
    return lambda x: np.where(x < eta_max, np.exp(log_interp(f(x))), 0)

def vectorize(X, Y):
    values = dict(zip(X, Y))
    return np.vectorize(lambda x: values[x])

vec2tuple = lambda v: (v[:,0], v[:,1], v[:,2])

eta, Y_0, Y_1, Omega_1, Theta_1, H_A, H_B, Y_a4 = np.loadtxt('../tables/kn-layer-hs.txt').T
Y_1 = interp(eta, Y_1)
Theta_1 = interp(eta, Theta_1)
Omega_1 = interp(eta, Omega_1)

reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(infile)
reader.Update()

x, y, z = get_points(reader)
T0 = get_point_data(reader, 'T0')
U1 = get_point_data(reader, 'U1')
p0 = get_point_data(reader, 'p0')[0]
gradT0 = get_point_data(reader, 'grad(T0)')

wall = np.where((y==0) & (z==0))
X = x[wall]
wallT0 = vectorize(X, T0[wall])
wallDTt = vectorize(X, gradT0[wall][:,0])
wallDTn = vectorize(X, gradT0[wall][:,1])
wallU1 = vectorize(X, U1[wall][:,0])

X = np.sort(X)
#print X, wallT0(X), wallDTt(X), wallDTn(X)

eta = p0/wallT0(x)*y/k
T = T0 - k*wallT0(x)/p0*wallDTn(x)*Theta_1(eta)
rho = p0/T0 + k/wallT0(x)*wallDTn(x)*Omega_1(eta)
U = k*U1/p0
U[:,0] += -k*np.sqrt(wallT0(x))/p0*wallDTt(x)*Y_1(eta)


pointsToVTK(os.path.splitext(outfile)[0], x, y, z, data = {
    'T': T, 'U': vec2tuple(U), 'rho': rho
})

