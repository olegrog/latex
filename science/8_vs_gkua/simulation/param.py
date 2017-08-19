#!/usr/bin/env python

import sys, os, json, glob
import numpy as np
from scipy.interpolate import interp1d

name = glob.glob('*.geo')[0].split('.')[0]
num = 4 if len(sys.argv) == 1 else int(sys.argv[1])

Kn = [ 1.8903 ]
Ma = [ 1.33 ]
base_korob, base_endtime = 2e1, 10e3   # crude: 5e1, 1e4   fine: 2e2, 5e2
rad, macro = 12, 100                     # crude: 8, 50      fine: 20, 50

def interp(x, y):
    log_interp = interp1d(np.log(x), np.log(y), kind='linear')
    return lambda x: np.exp(log_interp(np.log(x)))

x = np.logspace(np.log10(1e-2), np.log10(1e+1), num=4)
korob = interp1d(x, np.array([1,2,10,2])*base_korob*rad**2)
endtime = interp1d(x, np.array([3,4,1,1])*base_endtime)

for kn, ma in zip(Kn, Ma):
    dirname = '%.5f' % kn
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

    with open(name + '.kep', 'r') as f:
        data = json.load(f)
    data['printer']['dir']          = 'result/'
    data['printer']['savemacro']    = int(endtime(kn)/macro)
    data['printer']['savefuncfreq'] = int(endtime(kn)/macro)*macro
    data['num_steps']               = int(endtime(kn)/macro)*macro + 1
    data['gas']['rad']              = rad
    data['integral']['power']       = int(korob(kn))
    print "Kn = %.2e, Korob = %.2e, endTime = %.2e" % (kn, korob(kn), endtime(kn))
    with open(kep, 'wb') as f: 
        json.dump(data, f, indent=2, sort_keys=True)
    os.system('sed -i.sed s/_U_/%.6e/g %s' % (ma*np.sqrt(5./3), kep))
    os.system('sed -i.sed s/_Tb_/%.6e/g %s' % (1+ma**2/3, kep))
    os.system('../../tools/msh2kem.py %s %s' % (msh, kem))
    os.system('../../tools/attach_mesh.py %s %s %s' % (kep, kem, kei))

