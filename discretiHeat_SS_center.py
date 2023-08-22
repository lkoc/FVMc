# Create stencil and substitutions for the diffusion equation
'''
This is a spatially variable physic parameters of the heat equation:
 ∂T(x,y,t)/∂t = 
         (k(x,y)/ (rho(x,y)*C(x,y))*∇2T(x,y,t) + 
                  (1/(rho(x,y)*C(x,y))* ∇ k(x,y)∇T(x,y,t) + 
                           f(x,y,t) /((rho(x,y) * C(x,y))
'''

import time
from argparse import ArgumentParser
import numpy as np

from sympy import solve, Poly, Eq, Function, exp, simplify, expand, collect
from sympy import *
from sympy.abc import T, x, y, t, rho, F, C, K, h, s
T = Function('T')
F = Function('F')
rho = Function('rho')
C = Function('C')
K = Function('K')

dt = T(x, y, t).diff(t).as_finite_difference([t, t + s])

gradTx = T(x, y, t).diff(x).as_finite_difference([x - h, x, x + h])
gradTy = T(x, y, t).diff(y).as_finite_difference([y - h, y, y + h])
LapTx = T(x, y, t).diff(x, x).as_finite_difference([x - h, x, x + h])
LapTy = T(x, y, t).diff(y, y).as_finite_difference([y - h, y, y + h])
gradKx = K(x, y).diff(x).as_finite_difference([x - h/2, x, x + h/2,])
gradKy = K(x, y).diff(y).as_finite_difference([y - h/2, y, y + h/2])

# Ecuación
eqn = Eq(dt, K(x, y)/(rho(x, y)*C(x, y))*(LapTx + LapTy) +
         1/(rho(x, y)*C(x, y))*(gradKx * gradTx + gradKy * gradTy) +
         F(x, y)/(rho(x, y) * C(x, y)))

eqn = simplify(eqn)

stencil = solve(eqn, T(x, y, t + s))
stencil2 = poly(stencil[0], T(-h + x, y, t), T(h+x, y, t),
                T(x, -h + y, t), T(x, h + y, t), T(x, y, t))

# coeficietnes
# print(stencil2.all_coeffs)
fuente = stencil2.as_dict()[(0, 0, 0, 0, 0)]
Po = stencil2.as_dict()[(0, 0, 0, 0, 1)]
Pn = stencil2.as_dict()[(0, 0, 0, 1, 0)]
Ps = stencil2.as_dict()[(0, 0, 1, 0, 0)]
Pe = stencil2.as_dict()[(0, 1, 0, 0, 0)]
Pw = stencil2.as_dict()[(1, 0, 0, 0, 0)]

# Condicones de probabilidad
solve(Pw+Pe+Pn+Ps - 1, s/h**2)

# obtención de la ecuación original
limit(eqn.rhs, h, 0).doit()
simplify(limit(stencil[0], h, 0).doit())
