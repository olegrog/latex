import re
import numpy as np

def read_header(datafile):
    with open(datafile, 'r') as lines:
        for line in lines:
            words = re.findall("[a-zA-Z]+", line)
            if words and words[0].lower() == 'variables':
                return _parse_fields(line)

def read_ij_zones(datafile, fields):
    with open(datafile, 'r') as lines:
        total = sum(1 for line in lines)
    with open(datafile, 'r') as lines:
        for n, line in enumerate(lines):
            words = re.findall("[a-zA-Z]+", line)
            if words and words[0].lower() == 'zone':
                zone, skip = _parse_zone(line), n+1
                points = zone['i'] * zone['j']
                yield zone, np.genfromtxt(datafile, skip_header=skip,
                    skip_footer=total-skip-points, names=map(lambda f: f.lower(), fields))

def write_header(datafile, fields, title=''):
    with open(datafile, 'w') as f:
        if len(title):
            f.write('Title="{}"\n'.format(title))
        f.write('Variables="{}"\n'.format('","'.join(fields)))

def write_fe_zone(datafile, zone, data, triangles):
    zone['datapacking'] = 'point'
    zone['zonetype'] = 'fetriangle'
    zone['nodes'] = data.shape[0]
    zone['elements'] = triangles.shape[0]
    with open(datafile, 'a') as f:
        f.write('Zone {}\n'.format(', '.join(map(lambda k: '{}={}'.format(k, zone[k]), zone.keys()))))
        np.savetxt(f, data, fmt='%.7f')
        np.savetxt(f, triangles + 1, fmt='%d')


_parse_fields = lambda line: re.findall("\w+", line)[1:]

def _parse_zone(line):
    zone = dict(re.split(' *= *', pair) for pair in re.findall("[a-z]+ *= *\w+", line.lower()))
    try:
        zone['i'], zone['j'] = map(int, (zone['i'], zone['j']))
    except KeyError:
        raise NotImplementedError('Only ij-ordered grids are supported')
    if zone.get('f'):
        zone['datapacking'] = zone['f']
        del zone['f']
    return zone
