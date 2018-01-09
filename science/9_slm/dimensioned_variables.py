#!/usr/bin/env python

import argparse, pint
import numpy as np
from collections import namedtuple

parser = argparse.ArgumentParser(description='Solver for the plane Couette-flow problem')
parser.add_argument('-a', '--alloy', type=str, default='ss_316l.txt', metavar='<file>', help='powder characteristics')
parser.add_argument('-b', '--bed', type=str, default='yadroitsev10.txt', metavar='<file>', help='test bed characteristics')
args = parser.parse_args()

dict2nt = lambda name, dict: namedtuple(name, dict.keys())(*dict.values())

def load_data(system, filename):
    result = {}
    print '{} parameters:'.format(system)
    with open(filename) as lines:
        for line in lines:
            words = line.split()
            name, quantity = words[0], ureg.Quantity(float(words[1]), words[2])
            result[name] = quantity
            print ' - {0:20s} = {1:.3e}'.format(name, quantity)
    return dict2nt(system, result)

def to_dimensionless(quant, basis, dimensional=False):
    dim = lambda q: pint.formatter(q.dimensionality.items())
    quantities = { '_': dim(quant) }
    for system, quantity in basis:
        quantities['{0}.{1}'.format(type(system).__name__, quantity)] = dim(system._asdict()[quantity])
    powers = pint.pi_theorem(quantities)[0]
    #print '  ###:', powers
    power = powers.pop('_')
    if dimensional:
        quant = ureg('')
        power = -power
    return eval(pint.formatter(map(lambda item: (item[0], item[1]/power), powers.items())))*quant

ureg = pint.UnitRegistry(auto_reduce_dimensions=True, autoconvert_offset_to_baseunit=True)
alloy = load_data('alloy', args.alloy)
bed = load_data('bed', args.bed)
derived1 = dict2nt('derived1', {
    'total_time': bed.track_length/bed.scanning_speed
})

basis1 = [  # for diffusive time scale
    ( alloy, 'density' ),
    ( alloy, 'heat_capacity' ),
    ( alloy, 'conductivity100' ),
    ( bed, 'beam_radius' ),
]
basis2 = [  # for absolute temperatures
    ( alloy, 'conductivity100' ),
    ( bed, 'beam_radius' ),
    ( bed, 'bulk_temperature' ),
]
dimensional = [
    ( 'time', 'second', basis1 ),
]
dimensionless = [
    ( bed, 'scanning_speed', basis1 ),
    ( bed, 'laser_power', basis2 ),
    ( bed, 'layer_thickness', basis2 ),
    ( bed, 'track_length', basis2 ),
    ( alloy, 'liquidus', basis2 ),
    ( alloy, 'solidus', basis2 ),
] + [ ( derived1, d, basis1 ) for d in derived1._asdict() ]

print 'Dimension units:'
for name, unit, basis in dimensional:
    print ' - [{}]: {:.3e}'.format(name, to_dimensionless(ureg(unit), basis, dimensional=True))

print 'Dimensionless variables:'
for system, quantity, basis in dimensionless:
    print ' {:>8s}) {:20s}: {:.3e}'.format(type(system).__name__, quantity, to_dimensionless(system._asdict()[quantity], basis))
