import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as ticker
from numpy.typing import NDArray

def plot_EV(time:NDArray, E_atom:NDArray, V:NDArray, ax1):
    #ax1.set_aspect(1)

    color = 'tab:red'
    ax1.set_xlabel('time (ps)')
    ax1.set_ylabel('total energy (eV/atom)', color=color)
    ax1.plot(time, E_atom, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel(r'volume (nm$^3$)', color=color)  # we already handled the x-label with ax1
    ax2.plot(time, V/1000, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    ax1.xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(25))

    ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.01))
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.0025))
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.025))

# def plot_theta(time:NDArray, theta:NDArray, ax):
#     '''
#     theta is a list of tilt angle present in each timestep
#     '''
#     ax.plot(r_grid, gr)
#     ax.legend(leg, fontsize=5)
#     ax.set_ylabel(r'g(r) norm.')
#     ax.set_xlabel(r'r ($\AA{}$)')
#     ax.set_xlim([0,r_grid[-1] - 0.5])