#!/usr/bin/env python

import numpy as np
import sys, os, vtk
from scipy.interpolate import interp1d
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

_, infile, outfile, kn = sys.argv
k = np.sqrt(np.pi)*float(kn)/2

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

def interp(eta, F):
    log_interp = interp1d(eta, np.log(F), kind='linear')
    eta_max = 20
    f = lambda x: np.where(x < eta_max, x, eta_max)
    return lambda x: np.where(x < eta_max, np.exp(log_interp(f(x))), 0)

def vectorize(X, Y):
    return interp1d(X, Y, kind='cubic')

def add_array(data, F, name):
    vtk_array = numpy_to_vtk(F, array_type=vtk.VTK_FLOAT)
    vtk_array.SetName(name)
    data.AddArray(vtk_array)

def add_correction(x, y, T0, p0, U1, T, U):
    if k > 0:
        eta = p0/wallT0(x)*y/k
    else:
        eta = np.where(y==0, 0, 100)
    T = T + k*wallT0(x)/p0*wallDT0n(x)*Theta_1(eta)
    rho = p0/T0 + k/wallT0(x)*wallDT0n(x)*Omega_1(eta)
    V1 = U/p0
    V1[:,0] += -np.sqrt(wallT0(x))/p0*wallDT0t(x)*Y_1(eta) + wallT0(x)*k/p0**2*(
        np.sqrt(wallT0(x))*wallDDT0nt(x)*Y_a4(eta) + wallDU1nt(x)*Y_0(eta) )
    V = k*V1
    return rho, T, V, V1

def integrate_patch(x, y, patch):
    m = np.argsort(x[patch])
    return np.trapz(y[patch][m], x[patch][m])

### Read the Knudsen-layer functions
pwd = os.path.dirname(os.path.realpath(__file__))
eta, Y_0, Y_1, Omega_1, Theta_1, H_A, H_B, Y_a4 = np.loadtxt(os.path.join(pwd, '../tables/kn-layer-hs.txt')).T
Y_0 = interp(eta, Y_0)
Y_1 = interp(eta, Y_1)
Y_a4 = interp(eta, Y_a4)
Theta_1 = interp(eta, Theta_1)  # NB: Theta_1 = -Theta_1(eta)
Omega_1 = interp(eta, Omega_1)

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
x, y, z = get_points(out)
T = get_point_data(out, 'T')
U = get_point_data(out, 'U')
T0 = get_point_data(out, 'T0')
U1 = get_point_data(out, 'U1')
p0 = get_point_data(out, 'p0')[0]
gradT0 = get_point_data(out, 'wallGradT0n')
grad2T0 = get_point_data(out, 'wallGradGradT0nt')
gradU1 = get_point_data(out, 'wallGradU1nt')

### Find wall functions
wall = np.where((y==0) & (z==0))
X = x[wall]
alpha = 0.5
wallT0 = lambda x: 1 - alpha*np.cos(2*np.pi*x);
wallDT0t = lambda x: 2*np.pi*alpha*np.sin(2*np.pi*x);
wallDT0n = vectorize(X, gradT0[wall])
wallDDT0nt = vectorize(X, grad2T0[wall][:,0])
wallDU1nt = vectorize(X, gradU1[wall][:,0])

### Debug info
#X = np.sort(X)
#print X, wallT0(X), wallDT0t(X), wallDT0n(X), wallDDT0nt(X), wallDU1nt(X)

### Add point data
rho, T, V, V1 = add_correction(x, y, T0, p0, U1, T, U)
data = grid.GetPointData()
add_array(data, T, 'T')
add_array(data, V, 'U')
add_array(data, rho, 'rho')
add_array(data, V1, 'U/k')

### Integrate patches and print
bottom = np.where((y==0) & (z==0))
top = np.where((y==0.5) & (z==0))
print '%.6e %.5e %.5e %.5e %.5e' % (float(kn),
    integrate_patch(x, T, top),
    integrate_patch(x, V1[:,0], top),
    integrate_patch(x, T, bottom),
    integrate_patch(x, V1[:,0], bottom))

### Read cell data
x, y, z = get_cell_centers(out)
T0 = get_cell_data(out, 'T0')
U1 = get_cell_data(out, 'U1')
T = get_cell_data(out, 'T')
U = get_cell_data(out, 'U')

### Add cell data
rho, T, V, V1 = add_correction(x, y, T0, p0, U1, T, U)
data = grid.GetCellData()
add_array(data, T, 'T')
add_array(data, V, 'U')
add_array(data, rho, 'rho')
add_array(data, V1, 'U/k')

### Write a VTK-file
writer = vtk.vtkUnstructuredGridWriter()
writer.SetInputData(grid)
writer.SetFileName(outfile)
writer.Update()

