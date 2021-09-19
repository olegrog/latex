#!/usr/bin/env python

from sympy import *

z, r = symbols('z r')
Q1, Q2, Q11, Q12, Q21, Q22, Q3 = symbols('Q1 Q2 Q11 Q12 Q21 Q22 Q3')
zzzz, xxxx, xxzz, xzxz, xxyy, xyxy = symbols('zzzz xxxx xxzz xzxz xxyy xyxy')
r2, z2 = r**2, z**2
r4, z4 = r**4, z**4

system5 = Matrix((
    (1, 2, 2*z2, 4*z2, z4, zzzz),
    (9, 6, 6*r2, 4*r2, r4, 2*xxxx+zzzz + 4*xxzz+2*xxyy),
    (3, 12, 2*r2, 8*r2, r4, 2*xxxx+zzzz + 4*xzxz+2*xyxy),
    (3, 2, r2+3*z2, 4*z2, z2*r2, zzzz + 2*xxzz),
    (1, 4, 2*z2, r2+5*z2, z2*r2, zzzz + 2*xzxz)
))

system3 = Matrix((
    (3, 6*z2, z4, zzzz),
    (15, 10*r2, r4, 2*xxxx+zzzz + 4*xxzz+2*xxyy),
    (5, r2+7*z2, z2*r2, zzzz + 2*xxzz)
))

solution3 = solve_linear_system(system3, Q1, Q2, Q3)
for var in [Q2, Q3]:
    print(var, '=', solution3[var])
print('(7*Q2+r**2*Q3)*(3*z**2-r**2) =', simplify((7*solution3[Q2]+r2*solution3[Q3])*(3*z2-r2)))

print()

solution5 = solve_linear_system(system5, Q11, Q12, Q21, Q22, Q3)
for var in [Q22, Q3]:
    print(var, '=', solution5[var])
print('(7*Q22+r**2*Q3)*(3*z**2-r**2) =', simplify((7*solution5[Q22]+r2*solution5[Q3])*(3*z2-r2)))


