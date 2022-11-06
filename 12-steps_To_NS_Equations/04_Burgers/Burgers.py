import matplotlib
import matplotlib.pyplot as plt
import numpy
import matplotlib.font_manager
import time
import sys
import sympy
from sympy.utilities.lambdify import lambdify

x, nu, t = sympy.symbols('x nu t')
phi = (sympy.exp(-(x - 4 * t)**2 / (4 * nu * (t + 1))) +
       sympy.exp(-(x - 4 * t - 2 * sympy.pi)**2 / (4 * nu * (t + 1))))
phiprime = phi.diff(x)

print(phiprime)

u = -2 * nu * (phiprime / phi) + 4

ufunc = lambdify((t, x, nu), u)
print(ufunc(1, 4, 3))

matplotlib.rcParams['text.usetex'] = True

nx = 41
dx = 2*numpy.pi/(nx - 1)
nu = 0.07
sigma = 0.2
dt = dx*nu
nt = 20

x = numpy.linspace(0, 2 * numpy.pi, nx)
u = numpy.asarray([ufunc(0, x0, nu) for x0 in x])

u_ini = numpy.ones(nx)
u_ini = u.copy()
un = numpy.ones(nx)

with plt.style.context(['science', 'ieee', 'grid']):
    fig, ax = plt.subplots()
    ax.plot(x, u)

    # matplot formatting.
    ax.set_xlim([0, 2*numpy.pi])
    # ax.set_ylim([0.8, 2.2])
    ax.set_ylabel(r'$u$')
    ax.set_xlabel(r'$x$')
    fig.savefig('initial_state.png', dpi=300)

for n in range(nt):
    un = u.copy()
    for i in range(1, nx-1):
        u[i] = un[i] - dt * un[i]*(un[i] - un[i-1])/dx \
            + nu * dt/(dx**2)*(un[i+1]+un[i-1]-2*un[i])

    u[0] = un[0] - dt * un[0]*(un[0] - un[-2])/dx \
        + nu * dt/(dx**2)*(un[1]+un[-2]-2*un[0])
    u[-1] = u[0]

u_analytical = numpy.asarray([ufunc(nt*dt, x0, nu) for x0 in x])

with plt.style.context(['science', 'ieee', 'grid']):
    fig, ax = plt.subplots()
    ax.plot(x, u_analytical, label=r'analytical',
            color='black', marker='D', markersize=3, markevery=0.1, linestyle="-")
    ax.plot(x, u, label=r'numerical',  color='tab:red',
            marker='v', markersize=3, markevery=0.1, linestyle="--")

    # matplot formatting.
    ax.legend()
    # ax.legend(bbox_to_anchor=(1.0, 0.8, 0.3, 0.2), loc='upper left')
    ax.set_xlim([0, 2*numpy.pi])
    ax.set_ylim([0, 8])
    ax.set_ylabel(r'$u$')
    ax.set_xlabel(r'$x$')
    fig.suptitle('1D Burgers, step=20')
    fig.savefig('Burgers_step20.png', dpi=300)
