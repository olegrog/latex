#!/usr/bin/env python

import sys, os, json, glob
from scipy.interpolate import interp1d

name = glob.glob('*.geo')[0].split('.')[0].split('-')[0]
num = 4 if len(sys.argv) == 1 else int(sys.argv[1])

H = [ 30, 50, 70, 90, 105 ]
Kn = [ 2.48e-5, 4.4460e-4, 5.5130e-3, 0.13338, 1.8903 ]
Ma = [ 2.278, 3.190, 3.032, 2.403, 1.33 ]
Cut = [ -1, -1, 8.5, 7.5, 6.5 ]
idx = 2
H, Kn, Ma, Cut = [H[idx]], [Kn[idx]], [Ma[idx]], [Cut[idx]]

base_korob, base_endtime = 50e3, 10e3   # crude: 5e1, 1e4   fine: 2e2, 5e2
rad, macro = 8, 100                     # crude: 8, 50      fine: 20, 50

korob = lambda kn: base_korob #*rad**2
endtime = lambda kn: base_endtime

vec2str = lambda x: '(' + ' '.join( map(str, (x, 0., 0.)) ) + ')'

for h, kn, ma, cut in zip(H, Kn, Ma, Cut):
    dirname = '%dkm' % h
    filename = os.path.join(dirname, name)
    geo, msh, kep, kem, kei = [ filename + '.' + ext for ext in ['geo', 'msh', 'kep', 'kem', 'kei'] ]
    if os.path.isdir(dirname):
        if os.path.isfile(kep):
            continue
        else:
            print "Overwrite %s..." % dirname
    else:
        os.system('mkdir %s' % dirname)
    with open(geo, 'w') as f:
        print >> f, 'Kn = %.6e;' % kn
    os.system('cat ' + name + '.geo >> ' + geo)
    os.system('gmsh -3 -part %d -o %s %s > %s' % (num, msh, geo, os.path.join(dirname, 'log.gmsh')))

    T_B, U = 1 + ma**2/3, ma*(5./3)**.5
    with open(name + '.kep', 'r') as f:
        data = json.load(f)
    if not 'volume' in data['gas']: # for uniform grid
        cut = cut - 2
    data['printer']['dir']                      = 'result/'
    data['printer']['savemacro']                = int(endtime(kn)/macro)
    data['printer']['savefuncfreq']             = int(endtime(kn)/macro)*macro
    data['num_steps']                           = int(endtime(kn)/macro)*macro + 1
    data['gas']['rad']                          = rad
    data['gas']['cut']                          = cut*2**.5
    data['gas']['v']                            = vec2str(0*U/2)
    data['integral']['power']                   = int(korob(kn))
    data['boundary_conditions']['body']['T']    = str(T_B)
    data['boundary_conditions']['free']['u']    = vec2str(U)
    data['initial_conditions']['internal']['u'] = vec2str(U)

    print "Kn = %.2e, Korob = %.2e, endTime = %.2e" % (kn, korob(kn), endtime(kn))
    with open(kep, 'wb') as f: 
        json.dump(data, f, indent=2, sort_keys=True)
    os.system('../../tools/msh2kem.py %s %s' % (msh, kem))
    os.system('../../tools/attach_mesh.py %s %s %s' % (kep, kem, kei))

