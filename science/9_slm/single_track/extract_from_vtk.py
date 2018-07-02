#!/usr/bin/env python

import numpy as np
import sys, os, vtk, argparse
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

parser = argparse.ArgumentParser(description='extract data from vtk-file')
parser.add_argument('vtkfile', help='an input VTK file')
parser.add_argument('-o', '--output', metavar='filename', default=sys.stdout, help='output txt file')
parser.add_argument('-v', '--verbose', action='store_true', help='verbose messages')
args = parser.parse_args()

def get_point_data(out, name):
    vtk_array = out.GetOutput().GetPointData().GetArray(name)
    return vtk_to_numpy(vtk_array) if vtk_array else 0

def get_coords(out):
    vtk_array = out.GetOutput().GetPoints().GetData()
    return vtk_to_numpy(vtk_array)

### Read a VTK-file
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(args.vtkfile)
reader.ReadAllScalarsOn()
reader.Update()
plane = vtk.vtkPlane()
bounds = reader.GetOutput().GetBounds()
dim = 2 if bounds[4] == bounds[5] else 3
origin, normal = np.zeros(3), np.zeros(3)
origin[dim-1] = bounds[2*dim-1]
normal[dim-1] = 1
if args.verbose:
    print 'origin:', origin, 'normal:', normal
plane.SetOrigin(origin)
plane.SetNormal(normal)
planeCut = vtk.vtkCutter()
planeCut.SetInputConnection(reader.GetOutputPort())
planeCut.SetCutFunction(plane)
planeCut.GenerateTrianglesOff()
planeCut.Update()

### Read the coordinates
x = get_coords(planeCut)[:,0]
reorder = np.argsort(x)

### Read the field data
names, fields = [ 'x' ], [ x[reorder] ]
for i in xrange(reader.GetNumberOfScalarsInFile()):
    name = reader.GetScalarsNameInFile(i)
    field = get_point_data(planeCut, name)[reorder]
    names.append(name)
    fields.append(field)

np.savetxt(args.output, np.transpose(fields), fmt='%+.5e', header='%11s '*len(names) % tuple(names))

