# -*- coding: utf-8 -*-
"""
Created on Thu May 12 19:16:56 2022
tomado de:
https://github.com/gregorv/opencl_heat/blob/master/waerme.py
Benchmar : ejemplo 1 del paper:
A meshless model for transient heat conduction in functionally graded materials

https://www.researchgate.net/publication/226961451_A_meshless_model_for_transient_heat_conduction_in_functionally_graded_materials

Wang, H., Qin, Q. H., & Kang, Y. (2006). A meshless model for transient heat 
conduction in functionally graded materials. Computational mechanics, 38(1), 51-60.

@author: QU1267
"""
import sys
import time
import pyopencl as cl

from numba import jit, prange, njit
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from matplotlib import cm
from tqdm import tqdm


# Set "nopython" mode for best performance,
@jit(nopython=True, parallel=True, fastmath=True)
def bordes(u):
    # Set the boundary conditions
    # Bottom Diritchlet BC
    #u[-1, :] = u_bottom
    # Bottom Neumann
    # discretización central en espacio de segundo orden con nodo fantasma
    u_1 = u[-2, 1:-1] - 2 * delta_x * (u_bottom)  # Nodo fantasma
    u[-1, 1:-1] += gamma * (u[-1, 2:] + u[-1, 0:-2] +
                            u[-2, 1:-1] + u_1 - 4*u[-1, 1:-1])

    # Top Dirichlet BC
    #u[0, :] = u_top
    # Top Robin BC
    # u_1= u[1, 1:-1] - A_r * u[0,1:-1]  + B_r # Nodo fantasma
    #u[0,1:-1] = gamma*(u[0,2:]+u[0,0:-2]+u[1,1:-1]+u_1-4*u[0,1:-1])+u[0,1:-1]
    # Top Neumann
    # discretización central en espacio de segundo orden con nodo fantasma
    u_1 = u[1, 1:-1] - 2 * delta_x * (u_top)  # Nodo fantasma
    u[0, 1:-1] += gamma * (u[0, 2:] + u[0, 0:-2] +
                           u[1, 1:-1] + u_1 - 4*u[0, 1:-1])

    # Right Diritchlet BC
    u[:, -1] = u_right
    # Right Robin BC
    # u_1= u[1:-1,-2] - A_r * u[:,1:-1,-1]  + B_r # Nodo fantasma
    # Neumman
    # u_1= u[1:-1 , -2] - 2* delta_y * (u_right) # Nodo fantasma
    #u[1:-1,-1] += gamma*(u[2:,-1]+u[0:-2,-1]+u[1:-1,-2]+u_1-4*u[1:-1,-1])

    # left Neumann
    # discretización central en espacio de segundo orden con nodo fantasma
    # u_1= u[1:-1 , 1] - 2* delta_y * (u_left) # Nodo fantasma
    #u[1:-1, 0] = gamma * (u[2:,0] + u[0:-2, 0] + u[1:-1,1] + u_1 - 4*u[1:-1,0]) + u[1:-1,0]
    #u[1:-1, 0] = (1./2.) * ( u[1:-1, 1] + u_1 )
    ##u[i, j] = gamma *   (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4*u[i][j]) + u[i][j]

    # discretización de un solo lado en el espacio de segundo orden (LeVeque p.31)
    #u[1:-1, 0] = -u_left * delta_y + (4./3.) * u[1:-1, 1] - (1./3.) * u[1:-1, 2]

    # left Dirichlet
    #u[:, 0]=u_left
    # Neumman
    u_1 = u[1:-1, 1] - 2 * delta_y * (u_left)  # Nodo fantasma
    u[1:-1, 0] += gamma*(u[2:, 0]+u[0:-2, 0]+u[1:-1, 1]+u_1-4*u[1:-1, 0])

    # Esquinas
    u[0, 0] = (1./2.) * (u[1, 0] + u[0, 1])  # Esquina sup-izq
    u[0, -1] = (1./2.) * (u[1, -1] + u[0, -2])  # Esquina sup-der
    u[-1, -1] = (1./2.) * (u[-1, -2] + u[-2, -1])  # Esquina inf-derecha
    u[-1, 0] = (1./2.) * (u[-1, 1] + u[-2, 0])  # Esquina inf-izq

    return u


plate_length = 1.
plate_width = 1.

ni = int(2**10+1)  # 50
nj = int(2**10+1)  # 50
width, height = [ni, nj]
delta_x = np.double(plate_length/(ni))
delta_y = delta_x

kr = np.double(1.0)   # W/(m . K)
rho = np.double(1.0)  # kg/m^3
Cp = np.double(1.0)   # J /( kg . K)
alpha = kr/(rho*Cp)
K = np.double(kr)  # Corregir por Cp y densidad
delta_t = 0.75 * 0.25*(delta_x**2/alpha)
delta_t = min(0.1, delta_t)
t_final = 0.2
max_iter = int(t_final/delta_t)
gamma = np.double((alpha * delta_t) / (delta_x * delta_x))

# Boundary conditions Neumman
u_top = np.double(0.0)
# Boundary conditions Diritchlet
u_right = np.double(1.00)
# Boundary conditions Neumman
u_left = np.double(0.0)  # Neumman
# Boundary conditions Neumman
u_bottom = np.double(0.0)
# Initial condition everywhere inside the grid
u_initial = np.double(0.0)

tin = np.full((ni, nj), u_initial).astype(np.float32)

# Usa solución anterior como valor inicial s iestá disponible
try:
    # pass
    tin = np.load("tin_BNCH_EJ1_1025200ms.npy",
                  allow_pickle=True).astype(np.float32)
except:
    tin = np.full((ni, nj), u_initial).astype(np.float32)
tout = np.zeros_like(tin)
# alpha= k(x,y) * delta_t/(Cp(x,y) * rho(x,y)) --> (x,y)) alpha es "a"
init_conductivity = np.full((ni, nj), alpha*delta_t).astype(np.float32)
# f = q * delta_t/(rho*Cp) = q* delta_t *a/K
init_f = np.full((ni, nj), 0 * delta_t/(rho*Cp)).astype(np.float32)

source = False  # False
if source:
    x = np.arange(0, ni)
    y = np.arange(0, nj)
    cx = 5
    cy = 1
    r = 0.025

    mask = (x[np.newaxis, :]*delta_x-cx)**2 + \
        (y[:, np.newaxis]*delta_y-cy)**2 < r**2
    init_f[mask] = np.float32(0.1)
    mask = (x[np.newaxis, :]*delta_x-cx-2*r-0.15)**2 + \
        (y[:, np.newaxis]*delta_y-cy)**2 < r**2
    init_f[mask] = np.float32(0.1)
    mask = (x[np.newaxis, :]*delta_x-cx+2*r+0.15)**2 + \
        (y[:, np.newaxis]*delta_y-cy)**2 < r**2
    init_f[mask] = np.float32(0.1)


# Boundary conditions Robin
T_inf = 25.  # Celcius
h = 1.0  # W/(m**2 . K)
A_r = np.double(2*delta_x*h/K)
B_r = np.double(2*delta_x*h/K*T_inf)
M_Br = np.full_like(tin, B_r).astype(np.float32)
tin = bordes(tin)

ctx = cl.create_some_context()
queue = cl.CommandQueue(ctx)
mf = cl.mem_flags
cl_tin = cl.Buffer(ctx, mf.READ_WRITE
                   | mf.COPY_HOST_PTR, hostbuf=tin)
cl_tout = cl.Buffer(ctx, mf.READ_WRITE
                    | mf.COPY_HOST_PTR, hostbuf=tout)
cl_K = cl.Buffer(ctx, mf.READ_ONLY
                 | mf.COPY_HOST_PTR, hostbuf=init_conductivity)
cl_C = cl.Buffer(ctx, mf.READ_ONLY
                 | mf.COPY_HOST_PTR, hostbuf=init_f)
cl_Br = cl.Buffer(ctx, mf.READ_ONLY
                  | mf.COPY_HOST_PTR, hostbuf=M_Br)


prg = cl.Program(ctx, open('heat2d_BNCH1.cl').read()).build()
knl = prg.heat2d
tt = 0
for i in tqdm(range(max_iter)):
    knl(queue, tin.shape, None, cl_tin, cl_tout, cl_Br, cl_K, cl_C,
        gamma.astype(np.float32), delta_x.astype(np.float32),
        delta_y.astype(np.float32), A_r.astype(np.float32))
    # knl(queue, tin.shape, None, cl_tin, cl_tout, cl_Br,
    # gamma.astype(np.float32), A_r.astype(np.float32))
    # Swap variables
    tt += delta_t
    if i % 1000 == 0:
        print("t="+str(tt))
    cl_tin, cl_tout = cl_tout, cl_tin

cl.enqueue_copy(queue, tout, cl_tout)


def grafica(u):
    fig, ax = plt.subplots(figsize=(8, 6))
    mappable = ax.imshow(u, interpolation=None, cmap=plt.cm.jet)
    fig.colorbar(mappable, label="Temperature (°C)")
    ax.set_xlabel("x samples")
    ax.set_xlabel("y samples")
    fig.suptitle("Temperature map at final time")
    fig.tight_layout()
    fig.show()


grafica(tout)
# This is another useful figure, a contour plot
xvec = np.linspace(0.0, plate_length, ni)
yvec = np.linspace(0.0, plate_width, nj)
xmat, ymat = np.meshgrid(xvec, yvec)

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
cs = ax4.contour(yvec, xvec, tout, colors='k')
ax4.contourf(yvec, xvec, tout)
ax4.clabel(cs, inline=1, fontsize=12)
ax4.invert_yaxis()
ax4.set_ylabel("Depth")
ax4.set_xlabel("Width in Meters")
ax4.set_title("Temperature field as a contour")
ax4.set_aspect('equal', adjustable='box')
fig4.show()
np.save("tin_BNCH_EJ1_"+str(nj)+"400ms.npy", tout)

# # yvec2 = np.linspace ( 0.45*plate_length, 0.55*plate_length, int(0.1*nj))
# # xvec2 =  np.linspace ( 0, 0.15*plate_length, int(0.15*ni))
# # xmat2, ymat2 = np.meshgrid ( xvec, yvec )

# # tout_zoom= tout[0:int(0.15*nj),int(0.45*ni):int(0.55*ni)]

# # fig5 = plt.figure()
# # ax5  = fig5.add_subplot(111)
# # cs = ax5.contour(yvec2,xvec2,tout_zoom,colors ='k')
# # ax5.contourf(yvec2,xvec2,tout_zoom)
# # ax5.clabel(cs,inline=1,fontsize=12)
# # ax5.invert_yaxis()
# # ax5.set_ylabel("Depth")
# # ax5.set_xlabel("Width in Meters")
# # ax5.set_title("Temperature field as a contour")
# # ax5.set_aspect('equal', adjustable='box')
# # fig5.show()

# yvec2 = np.linspace ( 0.49*plate_length, 0.51*plate_length, int(0.02*nj))
# xvec2 =  np.linspace (0.09*plate_length, 0.11*plate_length, int(0.02*ni))
# xmat2, ymat2 = np.meshgrid ( xvec2, yvec2 )

# tout_zoom= tout[int(0.09*nj+1):int(0.11*nj),int(0.49*ni):int(0.51*ni)]

# fig5 = plt.figure()
# ax5  = fig5.add_subplot(111)
# cs = ax5.contour(yvec2,xvec2,tout_zoom,colors ='k')
# ax5.contourf(yvec2,xvec2,tout_zoom)
# ax5.clabel(cs,inline=1,fontsize=12)
# ax5.invert_yaxis()
# ax5.set_ylabel("Depth")
# ax5.set_xlabel("Width in Meters")
# ax5.set_title("Temperature field as a contour")
# ax5.set_aspect('equal', adjustable='box')
# fig5.show()

# Solución teórica
# @jit(nopython=True, parallel = True , fastmath=True)


def u_teor(t, xvec):
    u_teor = 0
    for i in range(20000):
        mu = (2*i + 1)*np.pi/2
        u_teor += (-1)**(i)*4/((2*i+1)*np.pi) * \
            np.cos(mu*xvec) * np.exp((-mu**2)*t)
    return 1 - u_teor


#  Ploteo
# Assign variables to the y axis part of the curve
# Plotting both the curves simultaneously
fig5 = plt.figure()
plt.plot(xvec, tout[int(ni/2), :], color='r', label='Calculado')
plt.plot(xvec, u_teor(delta_t*max_iter, xvec), color='g', label='Teórico')

# Naming the x-axis, y-axis and the whole graph
plt.xlabel("x(m)")
plt.ylabel("T (°C)")
plt.title("Solución en t="+str(delta_t*max_iter)+" s")
# Adding legend, which helps us recognize the curve according to it's color
plt.legend()
# To load the display window
plt.show()

fig6 = plt.figure()
plt.plot(xvec, (1-tout[int(ni/2), :]/u_teor(delta_t*max_iter, xvec))*100,
         color='r', label='Calculado')
# Naming the x-axis, y-axis and the whole graph
plt.xlabel("x(m)")
plt.ylabel("error (%)")
plt.title("Solución en t="+str(delta_t*max_iter)+" s")
# Adding legend, which helps us recognize the curve according to it's color
plt.legend()
# To load the display window
plt.show()
