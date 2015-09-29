#!/usr/bin/env python

import os
import sys
import numpy as np

params = {
    'top_temp':     [ '\left.T\\right|_{y=\\frac12}', 0.495, 0.533, 0.39, 2, '10*y', int(True) ],
    'top_flow':     [ '\left.v_x\\right|_{y=\\frac12}', -4.5e-3, 6e-4, 0.40, 3, '10*y/p0*kn*sqrt(pi)/2', int(False) ],
    'bottom_temp':  [ '\left.T\\right|_{y=0}', 0.485, 0.508, 0.38, 4, '10*y', int(True) ],
    'bottom_flow':  [ '\left.v_x\\right|_{y=0}', 0, 7e-3, 0.39, 5, '10*y/p0*kn*sqrt(pi)/2', int(False) ],
}
dummies = [ 'name', 'expr', 'ymin', 'ymax', 'xcoord', 'column', 'asym', 'is_heat' ]

name = sys.argv[1].split('.')[0]
data = [ name ] + params[name]
    
print 'Creating %s.sh...' % name
filedata = None
with open('template.sh', 'r') as file:
    filedata = file.read().replace('<name>', name)

for i, dummy in enumerate(dummies):
    filedata = filedata.replace('<%s>' % dummy, str(data[i]))

output = '%s.sh' % name
with open(output, 'w') as file:
    file.write(filedata)
os.chmod(output, 0755)
