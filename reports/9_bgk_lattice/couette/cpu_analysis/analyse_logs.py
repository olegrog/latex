#!/usr/bin/env python
    
import glob
import pandas as pd
import pylab as py

def read_log(filename):
    N, m1, m2 = filename.split('.')[0].split('-')
    result = { 'N': N, 'model1': m1, 'model2': m2 }
    with open(filename, 'r') as lines:
        for line in lines:
            if 'User time' in line:
                result['cpu'] = float(line.split()[-1])
            if 'System time' in line:
                result['sys'] = line.split()[-1]
            if 'Percent of CPU' in line:
                result['percent'] = line.split()[-1]
            if 'wall clock' in line:
                result['clock'] = line.split()[-1]
            if 'Major (requiring I/O) page faults' in line:
                result['major'] = line.split()[-1]
            if 'Minor (reclaiming a frame) page faults' in line:
                result['minor'] = line.split()[-1]
            if 'raise' in line:
                print filename, lines
    return result

data = {}
for filename in glob.glob('*.log'):
    result = read_log(filename)
    N = int(result['N'])
    if not N in data:
        data[N] = {}
    if result['model1'] != result['model2']:
        data[N]['hybrid'] = result['cpu']
    elif result['model1'] == 'dvm':
        data[N]['dvm'] = result['cpu']
    else:
        data[N]['lbm'] = result['cpu']
    
df = pd.DataFrame(data).transpose()
py.plot(df.dvm/df.hybrid, '*-')
py.loglog()
py.xlabel('Number of cells in the Knudsen layer')
py.ylabel('Acceleration (t_DVM/t_hybrid)')
py.grid()
py.xlim(2, 1e3)
py.ylim(1, 15)
py.show()

