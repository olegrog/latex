#!/usr/bin/env python

import numpy as np
import sys, os, vtk, argparse
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

parser = argparse.ArgumentParser(description='extract patch data from vtk-file')
parser.add_argument('vtkfile', help='an input VTK file')
parser.add_argument('--output', metavar='filename', default=sys.stdout, help='output txt file')
parser.add_argument('--fields', type=lambda s: s.split(','), default='T,U', metavar='list', help='list of field names')
parser.add_argument('--equation', type=lambda s: lambda x,y: eval(s), metavar='expr', help='equation of the patch')
parser.add_argument('--phi', type=lambda s: lambda x,y: eval(s), metavar='expr', help='parameter of the patch')
parser.add_argument('--tau', type=lambda s: lambda x,y: eval(s), metavar='list', help='direction of the tangent vector')
parser.add_argument('--corr', type=lambda s: lambda x: eval(s), metavar='expr', help='correction of the vector field')
parser.add_argument('--exkind', metavar='type', help='extrapolation kind: linear | quad')
parser.add_argument('--kn', type=float, default=1, metavar='value', help='a Knudsen number for transform (x,y)-coords')
args = parser.parse_args()

def get_point_data(out, name):
    vtk_array = out.GetOutput().GetPointData().GetArray(name)
    return vtk_to_numpy(vtk_array)

def get_points(out):
    vtk_array = out.GetOutput().GetPoints().GetData()
    numpy_array = vtk_to_numpy(vtk_array) * args.kn
    return numpy_array[:,0], numpy_array[:,1], numpy_array[:,2]

def extrapolate(X, F):
    func = { 'xlog': lambda x: x*np.log(np.abs(x)), 'quad': lambda x: x**2 }[args.exkind]
    n = len(F)
    return np.linalg.solve(np.vstack((np.ones(n), X[1:1+n]-X[0], func(X[1:1+n]-X[0]))).T, F)[0]

def correct_marginal(X, F):
    n = 3
    if args.exkind:
        F[0] = extrapolate(X[0:1+n], F[1:1+n])
        F[-1] = extrapolate(X[-1:-2-n:-1], F[-2:-2-n:-1])
    return F

### Read a VTK-file
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(args.vtkfile)
reader.Update()
plane = vtk.vtkPlane()
plane.SetOrigin(reader.GetOutput().GetCenter())
plane.SetNormal(0, 0, 1)
planeCut = vtk.vtkCutter()
planeCut.SetInputConnection(reader.GetOutputPort())
planeCut.SetCutFunction(plane)
planeCut.GenerateTrianglesOff()
planeCut.Update()

### Read point data
x, y, z = get_points(planeCut)
patch = abs(args.equation(x,y)) < 1e-4
x, y = x[patch], y[patch]
phi, tau = args.phi(x,y), np.array(args.tau(x,y), dtype=float)
if len(np.array(tau).shape) == 1:
    tau = np.tile(tau.T, (len(x),1)).T
reorder = np.argsort(phi)
phi, tau = phi[reorder], tau[:,reorder]
tau /= np.linalg.norm(tau, axis=0)

names, fields = [ 'phi' ], [ phi ]
for name in args.fields:
    field = get_point_data(planeCut, name)[patch][reorder]
    if len(field.shape) == 1:   # scalar
        names += [ name ]
        fields += [ correct_marginal(phi, field) ]
    else:                       # vector
        names += [ '%s_t' % name ]
        Fx, Fy, Fz = field.T
        Fi = [Fx, Fy]
        #fields += [ correct_marginal(phi, args.corr(np.einsum('ij->j', tau*Fi))) ]
        fields += [ args.corr(np.einsum('ij->j', tau*Fi)) ]

np.savetxt(args.output, np.transpose(fields), fmt='%.5g', header='%8s'*len(names) % tuple(names))

