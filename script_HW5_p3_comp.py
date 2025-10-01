# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 13:52:57 2025

@author: coled
"""

# ICF HW6 Problem 3 

# Part a: Numerically solve for length until M = 1

import numpy as np
import matplotlib.pyplot as plt


# define function dm/dx

def dm_dx(M, f, x, gamma=1.4):
    c1 = (1 + (gamma-1)*0.5 * M **2) / ( (M**2 - 1))
    c2 = (gamma*M**2) * f / ( np.exp(1/(x+1)))
    c3 = -2 / ((x+1)**2)
    return M * c1 * (c3-c2)
    

dx = 0.001 #dx in mm
xmin = 0
xmax = 100
n_samples = int((xmax-xmin)/dx)


x = np.linspace(0, 100, n_samples)
M = np.zeros_like(x)
M[0] = 6

f = 0.005

tol = 1e-2

def find_choke(M=M, x = x, f = f):
    for i in range(len(x) -1):
        M[i+1] = M[i] + dx * dm_dx(M[i], f, x[i])
        if (abs(M[i] - 1) < tol) or (M[i]<1):
            i_break = i
            break
    return (x[i_break], M[:i_break])


L1  = find_choke()[0]
M_array = find_choke()[1]



# part b: residence time 

T0 = 600 # stag temp in K
p0 = 200 * 10**3 # stag pressure in pa
R = 287 # gas constant for air, j/kg*K
n2 = int(L1/dx)
x2 = np.linspace(0, L1, n2)

#initialize all relevant flow properties

u = np.zeros_like(x2)
A = np.zeros_like(x2)
T = np.zeros_like(x2)
rho = np.zeros_like(x2)
p = np.zeros_like(x2)
p0_array = np.zeros_like(x2)

T0_array = np.zeros_like(x2)
time = 0

def pressure_isen(M, gamma = 1.4):
   return (1 + 0.2*M**2) ** (-1 * (gamma/(gamma-1)))

def flow_property_functions(M, T0= T0, p0 = p0, R= R, gamma = 1.4):
    c1 = 1 + (gamma-1) * 0.5 * (M **2.0)
    T = 1/c1 * T0 # T in K
    u = M * np.sqrt(gamma*R*T) # u here in m/s
    props = {
        "temp": T,
        "u" : u}
    return props

# intialize the density constant for later on
props = flow_property_functions(M_array[0])
A_init = np.pi*np.exp(2)
u_init = props["u"]
p_init = p0 * pressure_isen(M_array[0])
T_init = props["temp"]
rho_init = p_init /( R * T_init)
rho_const = rho_init * u_init * A_init


# for loop to solve for properties and can use it to add up residence times
for i in range(len(M_array)):
    props = flow_property_functions(M_array[i])
    T[i] = props["temp"]
    T0_array[i] = T0 # adiabatic assumption
    u[i] = props["u"]
    A[i] = np.pi * np.exp(2/(x2[i]+1))
    rho[i] = rho_const / (u[i] * A[i])
    p[i] = rho[i] * R * T[i]
    p0_array[i] = p[i] * 1/pressure_isen(M_array[i])
    
    #res time
    time += dx / u[i]


fig1, axs = plt.subplots(2,3)

Mplot = axs[0,0].plot(x2, M_array)
axs[0,0].set_title('Mach Number')
p_plot = axs[0,1].plot(x2, p)
axs[0,1].set_title('Pressure, Pa')
p0_plot = axs[0,2].plot(x2, p0_array)
axs[0,2].set_title('Stagnation Pressure, Pa')
T0_plot = axs[1,0].plot(x2, T0_array)
axs[1,0].set_title('Stagnation T, K')
T_plot = axs[1, 1].plot(x2, T)
axs[1,1].set_title('T, K')

plt.show()


