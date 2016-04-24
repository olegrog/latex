#!/usr/bin/env python

import os
import sys
import numpy as np

params = {
    'top_temp':     [ '\\left.T\\right|_{y=\\frac12}',              0.495, 0.535, 0.39, 2, '5.15537e-01', '5.24818e-01', 1, 1               ],
    'top_flow':     [ '\\left.\\frac{v_x}{k}\\right|_{y=\\frac12}', -0.14, 0,     0.38, 3, '0',           '-1.2149e-01', 0, '2/sqrt(pi)/$1' ],
    'bottom_temp':  [ '\\left.T\\right|_{y=0}',                     0.485, 0.505, 0.38, 4, '0.5',         '5.00000e-01', 1, 1               ],
    'bottom_flow':  [ '\\left.\\frac{v_x}{k}\\right|_{y=0}',        0,     0.25,  0.37, 5, '0',           '1.94545e-01', 0, '2/sqrt(pi)/$1' ],
}
dummies = [ 'name', 'expr', 'ymin', 'ymax', 'xcoord', 'column', 'heat', 'snit', 'is_heat', 'corr' ]

name = sys.argv[1]
data = [ name ] + params[name]
    
print 'Creating %s.sh...' % name
filedata = None
with open('template.sh', 'r') as file:
    filedata = file.read().replace('<name>', name)

for i, dummy in enumerate(dummies):
    filedata = filedata.replace('<%s>' % dummy, str(data[i]))

output = '_%s.sh' % name
with open(output, 'w') as file:
    file.write(filedata)
os.chmod(output, 0755)
