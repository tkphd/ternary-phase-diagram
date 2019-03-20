# -*- coding: utf-8 -*-

# Numerical libraries
# Thermodynamics and computer-algebra libraries
from sympy import Matrix, diff, pprint, symbols
from sympy.abc import x, y
from sympy.utilities.codegen import codegen
from sympy.solvers import solve_linear_system

# Define SymPy variables

x1, x2 = symbols('x1 x2')
xo, xb, xc, xd = symbols('xo xb xc xd')
yo, yb, yc, yd = symbols('yo yb yc yd')

# Define paraboloid origins

xe1A, xe2A = (0.20, 0.20)
xe1B, xe2B = (0.60, 0.20)
xe1C, xe2C = (0.20, 0.60)

# Free Energies

skew = 0.001
GA = (x1 - xe1A)**2 + 1 * (x1 - xe1A) * (x2 - xe2A) + (x2 - xe2A)**2
GB = (x1 - xe1B)**2 + (1+skew) * (x1 - xe1B) * (x2 - xe2B) + (x2 - xe2B)**2
GC = (x1 - xe1C)**2 + (1-skew) * (x1 - xe1C) * (x2 - xe2C) + (x2 - xe2C)**2

# First Derivatives

dGAdx1 = diff(GA, x1)
dGAdx2 = diff(GA, x2)

dGBdx1 = diff(GB, x1)
dGBdx2 = diff(GB, x2)

dGCdx1 = diff(GC, x1)
dGCdx2 = diff(GC, x2)

# Second Derivatives

d2GAdx11 = diff(GA, x1, x1)
d2GAdx12 = diff(GA, x1, x2)
d2GAdx22 = diff(GA, x2, x2)

d2GBdx11 = diff(GB, x1, x1)
d2GBdx12 = diff(GB, x1, x2)
d2GBdx22 = diff(GB, x2, x2)

d2GCdx11 = diff(GC, x1, x1)
d2GCdx12 = diff(GC, x1, x2)
d2GCdx22 = diff(GC, x2, x2)

# Define lever rule equations, cf. TKR4p160-161
#
# For a ternary coexistence triangle, this system takes two line
# segments---AB from a vertex to the opposite edge and CD between the
# other two vertices---and returns the (xA, xB) coordinates of the
# intersection. This is critical for evaluation of the phase fractions.

levers = solve_linear_system(Matrix( ((yo - yb, xb - xo, xb * yo - xo * yb),
                                      (yc - yd, xd - xc, xd * yc - xc * yd)) ), x, y)

# Generate numerically efficient C-code

codegen(
    [   # Single-phase minima
        ('xe1A', xe1A), ('xe2A', xe2A),
        ('xe1B', xe1B), ('xe2B', xe2B),
        ('xe1C', xe1C), ('xe2C', xe2C),
        # Free energies
        ('GA', GA),
        ('GB', GB),
        ('GC', GC),
        # First derivatives
        ('dGAdx1', dGAdx1), ('dGAdx2', dGAdx2),
        ('dGBdx1', dGBdx1), ('dGBdx2', dGBdx2),
        ('dGCdx1', dGCdx1), ('dGCdx2', dGCdx2),
        # Second derivatives
        ('d2GAdx11', d2GAdx11), ('d2GAdx12', d2GAdx12), ('d2GAdx22', d2GAdx22),
        ('d2GBdx11', d2GBdx11), ('d2GBdx12', d2GBdx12), ('d2GBdx22', d2GBdx22),
        ('d2GCdx11', d2GCdx11), ('d2GCdx12', d2GCdx12), ('d2GCdx22', d2GCdx22),
        # Lever rule compositions
        ('xAintersect', levers[x]),
        ('xBintersect', levers[y])
    ],
    language='C', prefix='paraboloids', project='paraboloids', to_files=True)
