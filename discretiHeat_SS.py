# Create stencil and substitutions for the diffusion equation in steady state
'''
This is a spatially variable physic parameters of the heat equation:
(k(x,y)/ (rho(x,y)*C(x,y))*∇2T(x,y,t) + 
                    (1/(rho(x,y)*C(x,y))* ∇ k(x,y)∇T(x,y,t) + 
                                f(x,y,t) /((rho(x,y) * C(x,y)) = 0
'''

import time
from argparse import ArgumentParser
import numpy as np

from sympy import solve, Poly, Eq, Function, exp, simplify, expand, collect
from sympy import *
from sympy.abc import T, x, y, t, rho, F, C, K, h, s
T = Function('T')       # Temperatura
F = Function('F')       # Fuente
rho = Function('rho')   # Densidad
C = Function('C')       # Calor específico
K = Function('K')       # Conductividad térmica

# Ecuación de calor

gradTx = T(x, y).diff(x).as_finite_difference([x - h, x, x + h])
gradTy = T(x, y).diff(y).as_finite_difference([y - h, y, y + h])
LapTx = T(x, y).diff(x, x).as_finite_difference([x - h, x, x + h])
LapTy = T(x, y).diff(y, y).as_finite_difference([y - h, y, y + h])
gradKx = K(x, y).diff(x).as_finite_difference([x - h/2, x, x + h/2])
gradKy = K(x, y).diff(y).as_finite_difference([y - h/2, y ,y + h/2])

# Ecuación
eqn = Eq(K(x, y)/(rho(x, y)*C(x, y))*(LapTx + LapTy) +
         1/(rho(x, y)*C(x, y))*(gradKx * gradTx + gradKy * gradTy) +
         F(x, y)/(rho(x, y) * C(x, y)),0)

eqn = simplify(eqn)

stencil = solve(eqn, T(x, y))
stencil2 = poly(stencil[0], T(-h + x, y), T(h+x, y),
                T(x, -h + y), T(x, h + y))

# coeficients
# print(stencil2.all_coeffs)
fuente = stencil2.as_dict()[(0, 0, 0, 0)]
Ps = stencil2.as_dict()[(0, 0, 0, 1)]
Pn = stencil2.as_dict()[(0, 0, 1, 0)]
Pe = stencil2.as_dict()[(0, 1, 0, 0)]
Pw = stencil2.as_dict()[(1, 0, 0, 0)]

# Condicones de probabilidad
solve(Pw+Pe+Pn+Ps - 1, h**2)

# obtención de la ecuación original
limit(eqn.lhs, h, 0).doit()
simplify(limit(stencil[0], h, 0).doit())
