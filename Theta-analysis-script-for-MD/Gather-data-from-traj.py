#Name: Tina Mihm
#Date: 12/10/2024
#Description: Pulls the volume, temperature, timesteps, theta and structures and graphs volume vs the time

import csv

import numpy as np
import subprocess
import re
import pandas as pd
import math as math
from ase.geometry.analysis import Analysis
##import seaborn as sns
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.colors import *
from matplotlib.font_manager import fontManager, FontProperties
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from matplotlib.ticker import MaxNLocator
import os

#---------------------------
#rcParams['savefig.dpi'] = 100
#rcParams['figure.dpi'] = 200
###
### Size of the figure
###
ratio=(np.sqrt(5)-1)/2.0     # golden ratio
#ratio=1                     # square figure
plt.rcParams["figure.figsize"] = 3.37, (3.37)*ratio
matplotlib.rcParams.update({'errorbar.capsize': 2.5})
#rcParams[‘figure.figsize’] = 3.37, 3.37
fig = figure()
#--------------------------
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.io.vasp import write_vasp
from ase.io.xyz import write_xyz
from ase.io import read
from ase.io import write
from os import system


#########################################
##  System info  ##
########################################

traj = Trajectory('traj/md.traj')

System = "Original_MD_370K_1000ps"
psmax = 1000
Tmax = 300
start_ps = 1

#########################################
##  pull data from trajectory  ##
########################################

print(len(traj))
# print(traj[0])
atoms = traj[0]
angle = atoms.get_cell_lengths_and_angles()
print(angle)
# at = atoms.get_temperature()
# print(at)
# av = atoms.get_volume()
# print(av)

#write("Traj_info_0.xyz", traj[0], format='extxyz')
#for i in range(0, 250, 50):
#    write("Traj_info_"+str(i)+".xyz", traj[i], format='extxyz')
#    system("cat Traj_info_"+str(i)+".xyz >> Traj_info.xyz")
#for i in range(250, len(traj), 10):
#    write("Traj_info_"+str(i)+".xyz", traj[i], format='extxyz')
#    system("cat Traj_info_"+str(i)+".xyz >> Traj_info.xyz")

T = []
V = []
Theta_all = []
Theta_av_1 = []
Theta_av_3 = []

for i in range(0, len(traj)):
    atoms = traj[i]
    av = atoms.get_volume()
    V +=[av]

for i in range((len(traj) - start_ps), len(traj)):
    print(len(traj))
    print(traj[0])
    atoms = traj[i]
    at = atoms.get_temperature()
    angle = atoms.get_cell_lengths_and_angles()
    print(angle)
    pos = atoms.get_positions()
    print(pos)
    Cell = atoms.get_cell()
    ## get lattice parameters to scale coordinates 
    A = Cell[0][0]
    B = Cell[1][1]
    C = Cell[2][2]

    # write_vasp("POSCAR_test.vasp", atoms)

    L1 = []
    L2 = []
    L3 = []
    L4 = []

    L1_index = []
    L2_index = []
    L3_index = []
    L4_index = []

    lst = [list(x) for x in pos]

    for a in pos: 
        if a[2] < 5.0:
            print("Layer 1:")
            print(a)
            L1 +=[a]
            l = a.tolist()
            i = lst.index(l)
            if i < 288:
                L1_index +=[i]
        if a[2] > 5.0 and a[2] < 11.0: 
            print("Layer 2:")
            print(a)
            L2 +=[a]
            l = a.tolist()
            i = lst.index(l)
            if i < 288:
                L2_index +=[i]
        if a[2] > 11.0 and a[2] < 17.0: 
            print("Layer 3:")
            print(a)
            L3 +=[a]
            l = a.tolist()
            i = lst.index(l)
            if i < 288:
                L3_index +=[i]
        if a[2] > 17.0: 
            print("Layer 4:")
            print(a)
            L4 +=[a]
            l = a.tolist()
            i = lst.index(l)
            if i < 288:
                L4_index +=[i]

    print(len(L1), len(L2), len(L3), len(L4))
    print(len(L1_index),len(L2_index), len(L3_index), len(L4_index))
    print(L4_index)

    an = Analysis(atoms)
    bonds = an.get_bonds("Ge", "Se", unique=True)
    bonds2 = an.get_values(bonds)
    bonds3 = an.get_values(bonds, vector=True)
    print("This is bond info:")
    print(bonds)
    print(bonds2)
    print(len(bonds[0]))
    print(len(bonds2[0]))

    b_index = []
    Theta_1 = []
    Theta_2 = []
    Theta_3 = []
    Theta_4 = []
    t = 1

    for i in L1_index: 
        print("Layer 1")
        # print(i)
        for x in bonds[0]:
            # print(x)
            if i ==x[0]:
                By_diff = abs(pos[x[1]][2] - pos[i][2])
                print(By_diff)
                if By_diff > 2: 
                    print("for bond:", x )
                    # print(t)
                    t = t+1
                    ind = bonds[0].index(x)
                    b_index+=[ind]
                    # print("for bond:", x )
                    # print(pos[x[0]][2], pos[x[1]][2])
                    # print(bonds2[0][ind])
                    bond = bonds2[0][ind]
                    # Diff = pos[x[1]] - pos[i]
                    # print(Diff)
                    # bond = np.sqrt(Diff[0]**2 + Diff[1]**2 + Diff[2]**2 )
                    # print(bond)
                    B_Delta_b = pos[x[1]][1] - pos[i][1]
                    # print(B_Delta_b)
                    if B_Delta_b > 4: 
                        B_Delta_b = B_Delta_b - Cell[1][1]
                    if B_Delta_b < -4: 
                        B_Delta_b = B_Delta_b + Cell[1][1]
                    # print(B_Delta_b)
                    B_theta = math.asin(B_Delta_b/bond)
                    print("this is the theta for the B point in the double well:", round(B_theta, 2))
                    Theta_1 +=[B_theta]
                    Theta_all +=[B_theta]
    # print("theta vs bonds:", len(Theta), len(bonds2[0]))
    for i in L2_index: 
        print("Layer 2")
        # print(i)
        for x in bonds[0]:
            # print(x)
            if i ==x[0]:
                By_diff = abs(pos[x[1]][2] - pos[i][2])
                print(By_diff)
                if By_diff > 2: 
                    # print(t)
                    t = t+1
                    ind = bonds[0].index(x)
                    b_index+=[ind]
                    print("for bond:", x )
                    # print(bonds2[0][ind])
                    bond = bonds2[0][ind]
                    # Diff = pos[x[1]] - pos[i]
                    # print(Diff)
                    # bond = np.sqrt(Diff[0]**2 + Diff[1]**2 + Diff[2]**2 )
                    # print(bond)
                    B_Delta_b = pos[x[1]][1] - pos[i][1]
                    # print(B_Delta_b)
                    if B_Delta_b > 4: 
                        B_Delta_b = B_Delta_b - Cell[1][1]
                    if B_Delta_b < -4: 
                        B_Delta_b = B_Delta_b + Cell[1][1]
                    # print(B_Delta_b)
                    B_theta = math.asin(B_Delta_b/bond)
                    print("this is the theta for the B point in the double well:", round(B_theta, 2))
                    Theta_2 +=[B_theta]
                    Theta_all +=[B_theta]
    for i in L3_index: 
        print("Layer 3")
        # print(i)
        for x in bonds[0]:
            # print(x)
            if i ==x[0]:
                By_diff = abs(pos[x[1]][2] - pos[i][2])
                print(By_diff)
                if By_diff > 2: 
                    # print(t)
                    t = t+1
                    ind = bonds[0].index(x)
                    print("for bond:", x)
                    b_index+=[ind]
                    # print(bonds2[0][ind])
                    bond = bonds2[0][ind]
                    # Diff = pos[x[1]] - pos[i]
                    # print(Diff)
                    # bond = np.sqrt(Diff[0]**2 + Diff[1]**2 + Diff[2]**2 )
                    # print(bond)
                    B_Delta_b = pos[x[1]][1] - pos[i][1]
                    # print(B_Delta_b)
                    if B_Delta_b > 4: 
                        B_Delta_b = B_Delta_b - Cell[1][1]
                    if B_Delta_b < -4: 
                        B_Delta_b = B_Delta_b + Cell[1][1]
                    # print(B_Delta_b)
                    B_theta = math.asin(B_Delta_b/bond)
                    print("this is the theta for the B point in the double well:", round(B_theta, 2))
                    Theta_3 +=[B_theta]
                    Theta_all +=[B_theta]
    for i in L4_index: 
        print("Layer 4")
        # print(i)
        for x in bonds[0]:
            # print(x)
            if i ==x[0]:
                By_diff = abs(pos[x[1]][2] - pos[i][2])
                print(By_diff)
                if By_diff > 2: 
                    # print(t)
                    t = t+1
                    ind = bonds[0].index(x)
                    b_index+=[ind]
                    print("for bond:", x )
                    # print(pos[x[0]][2], pos[x[1]][2])
                    # print(bonds2[0][ind])
                    bond = bonds2[0][ind]
                    # Diff = pos[x[1]] - pos[i]
                    # print(Diff)
                    # bond = np.sqrt(Diff[0]**2 + Diff[1]**2 + Diff[2]**2 )
                    # print(bond)
                    B_Delta_b = pos[x[1]][1] - pos[i][1]
                    # print(B_Delta_b)
                    if B_Delta_b > 4: 
                        B_Delta_b = B_Delta_b - Cell[1][1]
                    if B_Delta_b < -4: 
                        B_Delta_b = B_Delta_b + Cell[1][1]
                    # print(B_Delta_b)
                    B_theta = math.asin(B_Delta_b/bond)
                    print("this is the theta for the B point in the double well:", round(B_theta, 2))
                    Theta_4 +=[B_theta]
                    Theta_all +=[B_theta]
    #print(Theta)
    # print(len(Theta), len(bonds2[0]))
    # print(mean(Theta))
    Theta_av_1 += [mean(Theta_1)]
    Theta_av_1 += [mean(Theta_2)]
    Theta_av_1 += [mean(Theta_3)]
    Theta_av_1 += [mean(Theta_4)]
    Theta_av_3 += [mean(Theta_all)]
    print(len(Theta_all))
    T +=[at]
    #V +=[av]

print("Average Theta for Layer 1:", mean(Theta_1))
print("Average Theta for Layer 2:", mean(Theta_2))
print("Average Theta for Layer 3:", mean(Theta_3))
print("Average Theta for Layer 4:", mean(Theta_4))

print(len(Theta_av_1))
print(len(Theta_av_3))
print(len(T))

for i in Theta_av_1: 
    print(i)

Theta_av = mean(Theta_all)
Theta_av_2 = mean(Theta_av_1)
print("This is the average theta for the last x ps", Theta_av)
print("This is the average theta for the last x ps", Theta_av_2)

#for i in range(0, 250, 50):
#    print(i)
#    atoms = traj[i]
#    a = Atoms(atoms)
#    write_vasp("POSCAR_traj_"+str(i), a)
#
#for i in range(250, len(traj), 10):
#    print(i)
#    atoms = traj[i]
#    a = Atoms(atoms)
#    write_vasp("POSCAR_traj_"+str(i), a)

cwd = os.getcwd()
dir_path = cwd + r"/log"

Time = []
Etot = []
Epot = []
Ekin = []
Temp = []

for fl in os.listdir(dir_path):
    file = os.path.join(dir_path, fl)
    with open(file, 'r') as fp:
        for l_no, line in enumerate(fp):
            if l_no == 0:
                continue
            #print(l_no)
            #print(line)
            str = line.split()
            #print(str)
            Time += [float(str[0])]
            Etot += [float(str[1])]
            Epot += [float(str[2])]
            Ekin += [float(str[3])]
            Temp += [float(str[4])]
#########################################
##  Save data to a .csv  ##
########################################

start = len(traj) - start_ps

#df2 = {"Time[ps]": Time, "Volume (A^3)": Vol_sort, "Pressure (GPa)": P_GPa_sort, "T (K)": Temp, "E_tot (eV)": Etot, "E_pot (eV)": Epot, "E_k (eV)": Ekin }
df2 = {"Time[ps]": Time[start:], "T (K)": T, "V (A^3)":V[start:], "Av Theta (rad)": Theta_av_3, "E_tot (eV)": Etot[start:], "E_pot (eV)": Epot[start:], "E_k (eV)": Ekin[start:] }
df2 = pd.DataFrame(df2)
df2.to_csv(System+'.csv', index=False)
#print(df2)


#########################################
##  Graph data  ##
########################################

### Graph of V vs Time ###

plt.figure(1)
plt.plot(Time, V, linestyle = "", marker = "o", markersize =2,  color = "#00ff00")
plt.xlim(0, psmax)
##plt.ylim(-1.5, 2.5)
plt.ylabel(r"Volume [$\AA^3$]")
plt.xlabel(r"Time [ps]")
plt.savefig(System+"-Volumes.png", bbox_inches='tight')
#plt.show()

#---------------------------

plt.figure(2)
plt.plot(Time, Temp, linestyle = "", marker = "o", markersize =2,  color = "#00ff00")
#plt.xlim(987, psmax)
#plt.ylim(Tmax - 15, Tmax + 15)
plt.ylabel(r"Temp [K]")
plt.xlabel(r"Time [ps]")
plt.savefig(System+"-Temp.png", bbox_inches='tight')

#-----------------------------

plt.figure(3)
plt.plot(T, Theta_av_3, linestyle = "-", marker = "o", markersize =2,  color = "#00ff00")
#plt.xlim(0, Tmax)
##plt.ylim(-1.5, 2.5)
plt.ylabel(r"Average Theta [rad]")
plt.xlabel(r"Temp [K]")
plt.savefig(System+"-Theta_vs_Temp.png", bbox_inches='tight')

#-----------------------------

plt.figure(4)
plt.plot("Time[ps]", "Av Theta (rad)",data = df2, linestyle = "-", marker = "o", markersize =2,  color = "#00ff00")
#plt.xlim(0, psmax)
##plt.ylim(-1.5, 2.5)
plt.ylabel(r"Theta [rad]")
plt.xlabel(r"Time [ps]")
plt.savefig(System+"-Theta_vs_Time.png", bbox_inches='tight')
