#!/usr/bin/env python

import os
import sys
import numpy as np

params = {
    'Pxy':      [ 'Pxy', 'center bottom', '-P_{xy}/\\Delta{v}', -0.3, 0.6, '-y', int(True) ],
    'vx':       [ 'vx', 'center top', 'v_x/\\Delta{v}', 0, 0.5, 'y', int(True) ],
    'qx':       [ 'qx', 'center top', '-q_x/(\\Delta{v})^2', 0, 0.4, '-y/v', int(True) ],
    'qy':       [ 'qy', 'center top', '-q_y/(\\Delta{v})^2', 0, 0.16, '-y/v', int(True) ],
    'Pxx':      [ 'Pxx', 'center bottom', '(P_{xx}-P_{yy})/(\\Delta{v})^2', -0.2, 0.5, 'y', int(False) ],
    'Pyy':      [ 'Pyy', 'center bottom', 'P_{yy}/(\\Delta{v})^2', -0.03, 0.07, 'y', int(False) ],
    'Pzz':      [ 'Pzz', 'center bottom', '(P_{zz}-P_{yy})/(\\Delta{v})^2', -0.02, 0.04, 'y', int(False) ],
    'tau':      [ 'tau', 'center bottom', '\\tau/(\\Delta{v})^2', 0, 0.18, 'y', int(False) ],
    'P':        [ 'P', 'center bottom', 'P/(\\Delta{v})^2', 0, 0.18, 'y', int(False) ],
    'omega':    [ 'omega', 'center top', '\\omega/(\\Delta{v})^2', -0.02, 0.04, 'y', int(False) ]
}
dummies = [ 'macro', 'place', 'ylabel', 'ymin', 'ymax', 'base', 'filter' ]

macro = sys.argv[1].split('.')[0]
data = params[macro]
    
print 'Creating %s.sh...' % macro
filedata = None
with open('template.sh', 'r') as file:
    filedata = file.read().replace('<name>', macro)

for i, dummy in enumerate(dummies):
    filedata = filedata.replace('<%s>' % dummy, str(data[i]))

output = '%s.sh' % macro
with open(output, 'w') as file:
    file.write(filedata)
os.chmod(output, 0755)
