#!/usr/bin/env python

import argparse, my_tecplot
import numpy as np
from matplotlib import tri

parser = argparse.ArgumentParser(description='2D refiner for Tecplot data')
parser.add_argument('src', help='source dat-file')
parser.add_argument('dest', help='destination dat-file')
parser.add_argument('-n', '--refine', default=1, type=int, help='divide each element into 4^n')
args = parser.parse_args()

fields = my_tecplot.read_header(args.src)
my_tecplot.write_header(args.dest, fields, 'Fields refined with the cubic interpolation')

for zone, data in my_tecplot.read_ij_zones(args.src, fields):
    I, J = zone['i'],  zone['j']
    base = np.tile(np.delete(np.arange(I*(J-1)), I-1 + I*np.arange(J-1)), 3).reshape(3,-1).T
    triangles = np.vstack(( base + (0, 1, I), base + (1, I, I+1) ))
    triang = tri.Triangulation(data['x'], data['y'], triangles)
    refiner = tri.UniformTriRefiner(triang)
    refi_triang, refi_data = refiner.refine_triangulation(subdiv=args.refine), []
    for field in fields:
        field = field.lower()
        if field.lower() in ['x', 'y']:
            refi_data.append(getattr(refi_triang, field))
        else:
            tri_interp = tri.CubicTriInterpolator(triang, data[field], kind='geom')
            refi_triang, refi_field = refiner.refine_field(data[field], subdiv=args.refine, triinterpolator=tri_interp)
            refi_data.append(refi_field)
    refi_data = np.vstack(refi_data).T
    del zone['i'], zone['j']
    my_tecplot.write_fe_zone(args.dest, zone, refi_data, refi_triang.triangles)


