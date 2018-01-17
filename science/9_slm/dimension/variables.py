#!/usr/bin/env python

import argparse, pint
import numpy as np
from collections import namedtuple

parser = argparse.ArgumentParser(description='Solver for the plane Couette-flow problem')
parser.add_argument('-a', '--alloy', type=str, default='ss_316.txt', metavar='<file>', help='powder characteristics')
parser.add_argument('-b', '--bed', type=str, default='yadroitsev10.txt', metavar='<file>', help='test bed characteristics')
parser.add_argument('-m', '--material', type=str, default='solid', metavar='<solid|liquid>', help='reference material')
parser.add_argument('-t', '--temp', type=str, default='bed.bulk_temperature', metavar='<expression>', help='reference temperature')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
args = parser.parse_args()

dict2nt = lambda name, dict: namedtuple(name, dict.keys())(*dict.values())
der_name = lambda name: name + '_der'

def load_data(system, filename):
    quantities = {}
    print '{} parameters (temperature in kelvins):'.format(system)
    with open(filename) as lines:
        for line in lines:
            if line.startswith('#'):
                continue
            words = line.split()
            #print ' #words', words
            name, quantity = words[0], ureg.Quantity(float(words[1]), words[2])
            subnames = name.split('-')
            if len(subnames) > 1:
                name = subnames[1]
                if subnames[0] != args.material:
                    continue
            if len(words) > 3:
                if name in quantities:
                    quantities[der_name(name)] = 0    # first derivative
                else:
                    quantities[name] = {}       # dict with 2 points
                temperature = ureg(words[3]).magnitude
                quantities[name][temperature] = quantity
                name = name + ' ({})'.format(temperature)
            else:
                quantities[name] = quantity
            print ' - {0:30s} = {1:.3e}'.format(name, quantity)
    return dict2nt(system, quantities)

def calc_derivatives(quantities, temp):
    new, temp = {}, temp.to_base_units()
    for name, quant in quantities._asdict().items():
        if type(quant) == dict:
            if len(quant) > 1:
                coeffs = np.polyfit(quant.keys(), [ q.magnitude for q in quant.values() ], 1)
                units = quant.values()[0].units
                new[name] = ureg.Quantity(np.poly1d(coeffs)(temp), units)
                new[der_name(name)] = ureg.Quantity(coeffs[0]*units/temp.units)
            else:
                new[name] = quant.values()[0]
    return quantities._replace(**new)

def to_dimensionless(quant, basis, dimensional=False):
    dim = lambda q: pint.formatter(q.dimensionality.items())
    quantities = { '_': dim(quant) }
    for system, quantity in basis:
        quantities['{0}.{1}'.format(type(system).__name__, quantity)] = dim(system._asdict()[quantity])
    powers = pint.pi_theorem(quantities)[0]
    power = powers.pop('_')
    expression = lambda sgn: pint.formatter([ (item[0], sgn*item[1]/power) for item in powers.items() ])
    if args.verbose:
        print '\t\t{}'.format('=' if dimensional else '*'), expression(-1)
    return eval(expression(-1 if dimensional else 1))*(ureg('') if dimensional else quant)

ureg = pint.UnitRegistry(auto_reduce_dimensions=True, autoconvert_offset_to_baseunit=True)
alloy = load_data('alloy', args.alloy)
bed = load_data('bed', args.bed)
alloy = calc_derivatives(alloy, eval(args.temp))

basis1 = [  # for diffusive time scale
    ( alloy, 'density' ),
    ( alloy, 'heat_capacity' ),
    ( alloy, 'conductivity' ),
    ( bed, 'beam_radius' ),
]
basis2 = [  # for absolute temperatures
    ( alloy, 'heat_capacity' ),
    ( alloy, 'conductivity' ),
    ( bed, 'beam_radius' ),
    ( bed, 'bulk_temperature' ),
]
basis3 = [  # for derivatives
    ( alloy, 'density' ),
    ( alloy, 'heat_capacity' ),
    ( alloy, 'conductivity' ),
    ( bed, 'bulk_temperature' ),
]

derived = dict2nt('derived', {
    'total_time': bed.track_length/bed.scanning_speed
})
get_derivatives1 = lambda system: map(lambda q: (system, der_name(q), basis3),
        filter(lambda q: der_name(q) in system._asdict(), system._asdict().keys()))
dimensionless = [
    ( bed, 'scanning_speed', basis1 ),
    ( bed, 'laser_power', basis2 ),
    ( bed, 'layer_thickness', basis2 ),
    ( bed, 'track_length', basis2 ),
    ( alloy, 'liquidus', basis2 ),
    ( alloy, 'solidus', basis2 ),
    ( alloy, 'fusion_heat', basis2),
] + [ ( derived, d, basis1 ) for d in derived._asdict() ] + get_derivatives1(alloy)

get_derivatives2 = lambda system: map(lambda q: (q + '({})'.format(args.temp), system._asdict()[q], [(system, q)]),
        filter(lambda q: der_name(q) in system._asdict(), system._asdict().keys()))
dimensional = [
    ( 'time', 'second', basis1 ),
] + get_derivatives2(alloy)

print 'Dimension units:'
for name, unit, basis in dimensional:
    print ' - {:40s}: {:.3e}'.format('[{}]'.format(name), to_dimensionless(ureg.Quantity(unit), basis, dimensional=True))

print 'Dimensionless variables:'
for system, quantity, basis in dimensionless:
    print ' {:>8s}) {:20s}: {:.3e}'.format(type(system).__name__, quantity, to_dimensionless(system._asdict()[quantity], basis))
