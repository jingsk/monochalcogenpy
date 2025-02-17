#Name: Tina Mihm
#Date: 12/10/2024
#Description: Pulls the volume, temperature, timesteps, theta and structures and graphs volume vs the time

import numpy as np
import math as math
import matplotlib.pyplot as plt
from ase.io import read
import sys
from plotter import plot_EV, plot_theta
from tilt_angle import tilt_angle

def dashboard(traj_name, time_step, image_name = './dashboard.png'):
    traj = read(traj_name, ':')
    time = time_step * np.arange(len(traj)) #ps
    fig, axs = plt.subplots(
        1,2,
        figsize=[6,2.5], 
        width_ratios = [1,1],
        layout="constrained")

    E_atom = np.array([atoms.get_total_energy() for atoms in traj])/len(traj[0])
    V = np.array([atoms.get_volume() for atoms in traj])
    plot_EV(time, E_atom, V, axs[0])
    for p,c in zip(['bc', 'ac'], ['tab:blue','tab:orange']):
        theta = [np.concatenate(list(tilt_angle(atoms, p).values())) for atoms in traj]
        plot_theta(time,theta,c,axs[1])
    #fig.tight_layout()
    #plt.subplots_adjust(left=0.25, bottom=0.2, right=0.8, top=0.8, wspace=None, hspace=None)
    fig.savefig(image_name, dpi=300)

if __name__=='__main__':
    time_step = float(sys.argv[2])
    dashboard(sys.argv[1], time_step, image_name = './dashboard.png')