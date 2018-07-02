#!/usr/bin/env python

import argparse, pint
import numpy as np
from collections import namedtuple  # needs for json syntax

parser = argparse.ArgumentParser(description='Solver for the plane Couette-flow problem')
parser.add_argument('-a', '--alloy', type=str, default='ss_316.txt', metavar='<file>', help='powder properties')
parser.add_argument('-b', '--bed', type=str, default='yadroitsev10.txt', metavar='<file>', help='test bed properties')
parser.add_argument('-c', '--const', type=str, default='constants.txt', metavar='<file>', help='physical constants')
parser.add_argument('-s', '--state', type=str, default='solid', metavar='<solid|liquid>', help='reference state')
parser.add_argument('-z', '--zero-temp', type=str, default='bed.bulk_temperature', metavar='<expression>', help='zero temperature')
parser.add_argument('-u', '--unit-temp', type=str, default='(alloy.solidus+alloy.liquidus)/2', metavar='<expression>', help='reference temperature')
parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
args = parser.parse_args()

dict2nt = lambda name, dict: namedtuple(name, dict.keys())(*dict.values())
der_name = lambda name: name + '_der'
jump_name = lambda name: name + '_jump'
state_name = lambda name, state: name + '_' + state if len(state) else name

def load_data(system, filename):
    quantities = {}
    print '{} parameters (temperature in kelvins) read from "{}":'.format(system, filename)
    with open(filename) as lines:
        for line in lines:
            if line.startswith('#'):
                continue
            words = line.split()
            #print ' #words', words
            name, quantity, state = words[0], ureg.Quantity(float(words[1]), words[2]).to_base_units(), ''
            subnames = name.split('-')
            if len(subnames) > 1:
                state, name = subnames
                sname = state_name(name, state)
            if len(words) > 3:
                temperature = ureg(words[3]).magnitude
                if name in quantities:
                    if not state:
                        quantities[der_name(name)] = 0      # register the first derivative of a state-independent quantity
                else:
                    quantities[name] = {}                   # register a new quantity
                if state:
                    if state in quantities[name]:
                        quantities[der_name(sname)] = 0     # register the first derivative for the specific state
                    else:
                        quantities[name][state] = {}        # register a new state of the quantity
                    if state != args.state:
                        quantities[jump_name(name)] = 0     # register the quantity jump at the phase transition
                    quantities[name][state][temperature] = quantity
                else:
                    quantities[name][temperature] = quantity
                name = '{} ({})'.format(name, temperature)
                if state:
                    name = '{} for {}'.format(name, state)
            else:
                quantities[name] = quantity
            print ' - {0:40s} = {1:.3e}'.format(name, quantity)
    return dict2nt(system, quantities)

def calc_derivatives(quantities, base_temp, jump_temp):
    def calc_der(quant, state=''):
        if len(quant) == 1:
            value = quant.values()[0]       # impossible to calculate derivative for 1 point
        else:
            coeffs = np.polyfit(quant.keys(), [ q.magnitude for q in quant.values() ], 1)
            units = quant.values()[0].units
            calc_value = lambda temp: ureg.Quantity(np.poly1d(coeffs)(temp), units)
            value = calc_value(base_temp)
            at_jump_temp[state] = calc_value(jump_temp)
            sname = state_name(name, state)
            new[der_name(sname)] = ureg.Quantity(coeffs[0]*units/base_temp.units)
        if state == args.state or not state:
            new[name] = value

    new, at_jump_temp = {}, {}
    for name, quant in quantities._asdict().items():
        if type(quant) != dict:
            continue
        if type(quant.keys()[0]) == str and len(quant) == 1:
            quant = quant.values()[0]               # remove state if only one is provided
        if type(quant.keys()[0]) != str:
            calc_der(quant)
        else:
            for state in quant.keys():
                calc_der(quant[state], state)
            for state in quant.keys():
                if state != args.state:
                    new[jump_name(name)] = at_jump_temp[state] - at_jump_temp[args.state]
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
const = load_data('const', args.const)

alloy = calc_derivatives(alloy, eval(args.zero_temp), eval(args.unit_temp))
derived = dict2nt('derived', {
    'unit_temperature': eval(args.unit_temp) - eval(args.zero_temp),
    'total_time': bed.track_length/bed.scanning_speed,
    'fusion_capacity': alloy.fusion_heat/(alloy.liquidus - alloy.solidus),
    'abs_temperature': eval(args.zero_temp),
    'solidus': alloy.solidus - eval(args.zero_temp),
    'liquidus': alloy.liquidus - eval(args.zero_temp),
    'fusion_delta': alloy.liquidus - alloy.solidus,
    'radiative_transfer': const.stefan_boltzmann,
})

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
    ( derived, 'unit_temperature' ),
]
basis3 = [  # for derivatives
    ( alloy, 'density' ),
    ( alloy, 'heat_capacity' ),
    ( alloy, 'conductivity' ),
    ( derived, 'unit_temperature' ),
]

get_derivatives = lambda system: map(lambda q: (system, q, basis3),
        filter(lambda q: q.split('_')[-1] == 'der' or q.split('_')[-1] == 'jump', system._asdict().keys()))
dimensionless = [
    ( bed, 'scanning_speed', basis1 ),
    ( bed, 'laser_power', basis2 ),
    ( bed, 'layer_thickness', basis2 ),
    ( bed, 'track_length', basis2 ),
    ( bed, 'convective_transfer', basis1),
    ( alloy, 'fusion_heat', basis2),
    ( derived, 'liquidus', basis2 ),
    ( derived, 'solidus', basis2 ),
    ( derived, 'fusion_delta', basis2 ),
    ( derived, 'fusion_capacity', basis1 ),
    ( derived, 'total_time', basis1 ),
    ( derived, 'radiative_transfer', basis2 ),
    ( derived, 'abs_temperature', basis2 ),
] + get_derivatives(alloy)

get_references = lambda system: map(lambda q: (q + '({})'.format(args.zero_temp), system._asdict()[q], [(system, q)]),
        filter(lambda q: der_name(q) in system._asdict(), system._asdict().keys()))
dimensional = [
    ( 'time', 'second', basis1 ),
] + get_references(alloy)

print 'Dimension units:'
for name, unit, basis in dimensional:
    print ' - {:38s}: {:.3e}'.format('[{}]'.format(name), to_dimensionless(ureg.Quantity(unit), basis, dimensional=True))

print 'Dimensionless variables:'
for system, quantity, basis in dimensionless:
    print ' {:>8s}) {:30s}: {:.3e}'.format(type(system).__name__, quantity, to_dimensionless(system._asdict()[quantity], basis))
