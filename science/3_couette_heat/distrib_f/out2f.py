#!/usr/bin/env python

import sys, math, array, struct, itertools
import json
import numpy as np

from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.interpolate import RectBivariateSpline, SmoothBivariateSpline

sys.path.append('/Users/olegrog/kesolver/tools')
from kepy.ximesh import read_ximesh
from out2 import readNodesCells, center

def sqr(x):
    return x*x

def init_plot_params(columnwidth, height):
    fig_width_pt = columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = height*inches_per_pt
    fig_size =  [fig_width,fig_height]
    params = {
                'backend': 'pdf',
                'axes.labelsize': 10,
                'font.size': 10,
                'legend.fontsize': 9,
                'xtick.labelsize': 7,
                'ytick.labelsize': 7,
                'text.usetex': False,
                'figure.figsize': fig_size,
                'font.family':'serif'
            }
    plt.rcParams.update(params)

width = 350
height = width * 0.75
init_plot_params(width, height)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
fig.subplots_adjust(left=-0.15, right=1.1, bottom=0, top=1.2)

# open .kei file
with open(sys.argv[1], 'rb') as fd:
    data = json.load(fd)

nodes, cells = readNodesCells(sys.argv[1])

symmetry, rad, circl, xyz, vol, r, d3v = read_ximesh(data)

with open(sys.argv[2], 'rb') as fd:
    while fd:
        data = fd.read(4+4)
        i, size = struct.unpack('=ii', data)

        print i, "size (file) = ", size

        print center([nodes[n] for n in cells[i].nodes])


        f = np.zeros_like(r)
        f.fill(np.NaN) #-1e-10)
        a = array.array('d')
        size = np.sum(circl)
        print "size (ximesh) = ", size
        a.fromfile(fd, size)

        if i == int(sys.argv[3]):
            f[circl] = np.array(a)
            f /= vol
            xp = xyz[0][:, :, rad]
            yp = xyz[1][:, :, rad]
            fp = f[:, :, rad]
            
            N = 200j
            stride = int(N.imag/50)
            xcut, ycut = np.max(xp), np.max(yp)
            cut, cut2 = 3.7, 4.0
            xcut, ycut = cut, cut
            offset = 1
            grid_x, grid_y = np.mgrid[-xcut:xcut:N, -ycut:ycut:N]
            base_x, base_y = np.linspace(-1, 1, 91), np.linspace(-1, 1, 121)
            grid_x, grid_y = np.meshgrid(base_x*cut, np.sign(base_y)*np.abs(base_y)**1.5*cut, sparse=False, indexing='ij')
            X = xyz[0][:,rad,rad]
            Y = xyz[1][rad,:,rad]
            fnew = np.nan_to_num(f)
            smooth = 0
            if len(sys.argv) > 5:
                smooth = sys.argv[5]
            spline = RectBivariateSpline(X, Y, fnew[:, :, rad], kx=3, ky=3, s=smooth)
            grid_z = np.zeros_like(grid_x)
            grid_z.fill(np.NaN)
            mask = grid_x**2 + grid_y**2 < cut**2
            grid_z[mask] = np.maximum(spline(grid_x[mask], grid_y[mask], grid=False), 0)
            zmin, zmax = -0.075, np.round(np.nanmax(grid_z)*1.01, 1-int(np.log(np.nanmax(grid_z))))
            print np.nanmin(grid_z), np.nanmax(grid_z), zmax
            ax.set_xlim(-cut2, cut2)
            ax.set_ylim(-cut2, cut2)
            ax.set_zlim(zmin, zmax)
            #surf = ax.plot_wireframe(xp, yp, fp, color='k', linewidth=0.5, rstride=1, cstride=1)
            #surf = ax.plot_wireframe(grid_x, grid_y, grid_z, color='k', linewidth=0.5, rstride=stride, cstride=stride)
            cset = ax.contourf(grid_x, grid_y, grid_z, zdir='z', offset=zmin, cmap=plt.cm.jet,
                levels=np.linspace(0, zmax, 20))
            #cset = ax.contourf(grid_x, grid_y, grid_z, zdir='x', offset=-offset-cut, cmap=plt.cm.jet)
            #cset = ax.contourf(grid_x, grid_y, grid_z, zdir='y', offset=offset+cut, cmap=plt.cm.jet)
            surf = ax.plot_surface(grid_x, grid_y, grid_z, cmap=plt.cm.jet, rstride=stride, cstride=stride,
                linewidth=0.25, vmin=0, vmax=zmax, alpha=0.8)

            cbaxes = fig.add_axes([0.8, 0.3, 0.03, 0.6]) 
            cbar = fig.colorbar(surf, cax=cbaxes, drawedges=False)
            cbar.solids.set_edgecolor("face")
            break

ax.set_xlabel(r'$\xi_x$')
ax.set_ylabel(r'$\xi_y$')

for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis:
    axis.pane.set_visible(False)
    axis.gridlines.set_visible(False)
    axis.set_rotate_label(False)
    axis.set_tick_params(which='both', direction='in', pad=-5)

ax.tick_params(which='both', direction='out', pad=-5) # this does not work

for elt in ax.w_zaxis.get_ticklines() + ax.w_zaxis.get_ticklabels():
    elt.set_visible(False)

ax.view_init(30, -20)
fig.savefig(sys.argv[4])

#plt.show()

