#!/usr/bin/env python

import numpy as np
import sys, os, vtk
from scipy.interpolate import griddata, interp1d
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

def add_correction(x, y, T0, p0, U1):
    eta = p0/wallT0(x)*y/k
    T = T0 - k*wallT0(x)/p0*wallDTn(x)*Theta_1(eta)
    rho = p0/T0 + k/wallT0(x)*wallDTn(x)*Omega_1(eta)
    U = k*U1/p0
    U[:,0] += -k*np.sqrt(wallT0(x))/p0*wallDTt(x)*Y_1(eta)
    return rho, T, U

### Read the Knudsen-layer functions
eta, Y_0, Y_1, Omega_1, Theta_1, H_A, H_B, Y_a4 = np.loadtxt('../tables/kn-layer-hs.txt').T
Y_1 = interp(eta, Y_1)
Theta_1 = interp(eta, Theta_1)
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
T0 = get_point_data(out, 'T0')
U1 = get_point_data(out, 'U1')
p0 = get_point_data(out, 'p0')[0]
gradT0 = get_point_data(out, 'grad(T0)')

### Find wall functions
wall = np.where((y==0) & (z==0))
X = x[wall]
wallT0 = vectorize(X, T0[wall])
wallDTt = vectorize(X, gradT0[wall][:,0])
wallDTn = vectorize(X, gradT0[wall][:,1])

### Debug info
#X = np.sort(X)
#print X, wallT0(X), wallDTt(X), wallDTn(X)

### Add point data
rho, T, U = add_correction(x, y, T0, p0, U1)
data = grid.GetPointData()
add_array(data, T, 'T')
add_array(data, U, 'U')
add_array(data, rho, 'rho')

### Read point data
x, y, z = get_cell_centers(out)
T0 = get_cell_data(out, 'T0')
U1 = get_cell_data(out, 'U1')

### Add cell data
rho, T, U = add_correction(x, y, T0, p0, U1)
data = grid.GetCellData()
add_array(data, T, 'T')
add_array(data, U, 'U')
add_array(data, rho, 'rho')

### Write a VTK-file
writer = vtk.vtkUnstructuredGridWriter()
writer.SetInputData(grid)
writer.SetFileName(outfile)
writer.Update()

