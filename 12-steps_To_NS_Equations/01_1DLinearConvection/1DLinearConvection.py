import sys
import time
import matplotlib.font_manager
import numpy
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True


nx = 41
dx = 2/(nx - 1)

nt = 25
dt = 0.025
c = 1

u = numpy.ones(nx)
u[int(0.5/dx):int(1/dx+1)] = 2

with plt.style.context(['science', 'ieee', 'grid']):
    fig, ax = plt.subplots()
    ax.plot(numpy.linspace(0, 2, nx), u)

    # matplot formatting.
    ax.set_xlim([0, 2])
    ax.set_ylim([0.8, 2.2])
    ax.set_ylabel(r'$u$')
    ax.set_xlabel(r'$x$')
    fig.savefig('initial_state.png', dpi=300)

u_ini = numpy.ones(nx)
u_ini = u.copy()
un = numpy.ones(nx)

for n in range(nt):
    un = u.copy()
    for i in range(1, nx):
        u[i] = (1 - c*dt/dx)*un[i] + c*dt/dx * un[i-1]

with plt.style.context(['science', 'ieee', 'grid']):
    fig, ax = plt.subplots()
    ax.plot(numpy.linspace(0, 2, nx), u_ini, label=r't=0',
            color='blue', marker='D', markersize=3, markevery=0.1, linestyle="-")
    ax.plot(numpy.linspace(0, 2, nx), un, label=r't=25',  color='tab:red',
            marker='v', markersize=3, markevery=0.1, linestyle="--")

    # matplot formatting.
    ax.legend()
    # ax.legend(bbox_to_anchor=(1.0, 0.8, 0.3, 0.2), loc='upper left')
    ax.set_xlim([0, 2])
    ax.set_ylim([0.8, 2.2])
    ax.set_ylabel(r'$u$')
    ax.set_xlabel(r'$x$')
    fig.suptitle('1D Linear Convection')
    fig.savefig('t25.png', dpi=300)
