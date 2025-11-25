# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 15:49:36 2025

@author: coled
"""

import numpy as np
# part d: variable friction


# initial functions needed to define friction 
def sutherland(T, T0 = 273.15, mu0 = 1.716e-5):
    mu = ((T/T0) ** 1.5 ) * ((T0 + 110.4) / (T + 110.4))
    return mu

def Re(mu, u, rho, D): 
    return rho*u*D / mu

def friction_func(D, rho, u, mu, epsilon = 0.001):
    den = np.log10((epsilon / (3.7 * D)) + (5.74/(Re(mu, u, rho, D)**0.9)) )
    return 0.0625/den

def pressure_isen(M, gamma = 1.4):
   return (1 + 0.2*M**2) ** (-1 * (gamma/(gamma-1)))

def temp_isen(M, gamma = 1.4):
    return 1 + (gamma-1) * 0.5 * (M **2.0)

def u_func(M, T, gamma=1.4, R = 287):
    u = M * (gamma*R*T)**0.5
    return u

def dm_dx(M, f, x, gamma=1.4):
    c1 = (1 + (gamma-1)*0.5 * M **2) / ( (M**2 - 1))
    c2 = (gamma*M**2) * f / ( np.exp(1/(x+1)))
    c3 = -2 / ((x+1)**2)
    return M * c1 * (c3-c2)

def rho_func(mdot, u, A):
    return mdot/(u*A)


dx = 0.001 #dx in mm
xmin = 0
xmax = 100
n_samples = int((xmax-xmin)/dx)


x = np.linspace(0, 100, n_samples)
M = np.zeros_like(x)
D = np.zeros_like(x)
u = np.zeros_like(x)
A = np.zeros_like(x)
T = np.zeros_like(x)
rho = np.zeros_like(x)
p = np.zeros_like(x)
p0_array = np.zeros_like(x)
T0_array = np.zeros_like(x)
f = np.zeros_like(x)
Re_array = np.zeros_like(x)
mu = np.zeros_like(x)

# fill in diameter and area arrays

for i in range(len(x)):
    D[i] = 2 * np.exp(1/(x[i]+1)) 
    A[i] = np.pi * (D[i] *0.5)**2

# find the values at the start of the problem

# givens

T0 = 600
p0 = 200 * 10**3
R = 287

#initialize
M[0] = 6
T[0] = T0 * temp_isen(M[0])
A[0] = np.pi * (D[0] *0.5)**2
u[0] = u_func(M[0], T[0])
p[0] = p0 * pressure_isen(M[0])
rho[0] = p[0] / (R*T[0])
mdot = rho[0]*u[0]*A[0]
mu[0] = sutherland(T[0])
Re_array[0] = Re(mu[0], u[0], rho[0], D[0])
f[0] = friction_func(D[0], rho[0], u[0], mu[0])


tol = 1e-3
# same find choke function as before, but need to update f on each iteration


def find_choke(M=M, x = x, f = f, mdot = mdot):
    for i in range(len(x) -1):
        M[i+1] = M[i] + dx * dm_dx(M[i], f[i], x[i])
        # temp is purely function of mach number
        T[i+1] = T0 * temp_isen(M[i+1])
        # now can find velocity
        u[i+1] = u_func(M[i+1], T[i+1])
        # now we can find density 
        rho[i+1] = rho_func(mdot, u[i+1], A[i+1])
        # now find mu
        mu[i+1] = sutherland(T[i+1])
        # now find Re
        Re_array[i+1] = Re(mu[i+1], u[i+1], rho[i+1], D[i+1])
        # now, finally, update
        f[i+1] = friction_func(D[i+1], rho[i+1], u[i+1], mu[i+1])
        if (abs(M[i] - 1) < tol) or (M[i]<1):
            i_break = i
            break
    return (x[i_break], M[:i_break], i_break)

L1 = find_choke()[0]
check_M = find_choke()[1]
i_break = find_choke()[2]
avg_f = sum(f[:i_break])/len(f[:i_break])