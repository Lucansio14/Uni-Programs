# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 11:22:24 2025
Solver of the 3D Lorenz equations using simple Forward Time Centered Space (FTCS) algorithm with initial conditions (x0,y0,z0).
More information about the equation and FTCS method in these websites: https://en.wikipedia.org/wiki/Lorenz_system
                                                                       https://en.wikipedia.org/wiki/FTCS_scheme
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
#Definition of general constants, variables and arrays
t0=0#Start time
T=100#End time
dt=0.01#Time-step size
N=int((T-t0)/dt)#Time discretitzation
t_array=np.linspace(t0,T,N)#Values of time
alpha=10#System parameter (Prandtl number or ratio between momentum diffusivity and thermal diffusivity)
rho=28#System parameter (Rayleigh number or the ratio between the time scale for diffusive thermal transport and the time scale for convective thermal transport at a speed v.
beta=8/3#System parameter (relates to the physical dimensions of the fluid layer itself)
x_array=np.zeros(N+1)#Values of x (proportional to the intensity of the convection or the rate of fluid flow)
y_array=np.zeros(N+1)#Values of y (proportional to the temperature difference between the rising and falling air currents)
z_array=np.zeros(N+1)#Values of time (proportional to the distortion of the vertical temperature profile from a linear one)
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Initial conditions or ICs (slight modifications of these can result in significantly different trajectories)
x0=0
y0=1
z0=1.05
x_array[0]=x0
y_array[0]=y0
z_array[0]=z0
#Solution evolution
for n in range(N):#FTCS
    x_array[n+1]=x_array[n] + dt*alpha*(y_array[n] - x_array[n])
    y_array[n+1]=y_array[n] + dt*(- x_array[n]*z_array[n] + rho*x_array[n] - y_array[n])
    z_array[n+1]=z_array[n] + dt*(x_array[n]*y_array[n] - beta*z_array[n])
#Plots
#x vs y vs z (Trayectory)
fig=plt.figure(figsize=(11,7))
ax=plt.axes(projection="3d")
ax.plot(x_array,y_array,z_array,"blue",lw=0.5)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title(f'3D Lorenz system trayectory for ICs ($x_0$,$y_0$,$z_0$) = (${x0}$,${y0}$,${z0}$)')
ax.grid()
fig.set_tight_layout(True)
plt.show()
print("Program execution time:",time.process_time()-start_time_program,"seconds.")
