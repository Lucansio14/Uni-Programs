# -*- coding: utf-8 -*-
"""
Created on Mon Aug 4 12:50:18 2025
Solver of the 1D Diffusion equation (\partial_{t}u = D\partial^{2}_{x^{2}}u) with the FTCS (Forward Time-Centered Space) explicit method and analytical solution.
More information about the equation in this website: https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/heat.pdf
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
#Definition of general constants, variables, arrays and functions
D=1#Collective diffusion coefficient (constant)
A_n=1#Arbitrary constant (due to adimensionalization)
x0=0#Start of the space domain
xf=1#End of the space domain
L=abs(xf-x0)#Length of the space domain
tn=0#Start time
tf=0.1#End time
Dx=0.05#Space-step size
Dt=0.5*((Dx**2)/(D))#Time-step size. For stability purposes, we use the condition (7.10) from the website to obtain this value.
Jx=int((xf-x0)/Dx)#Number of spatial points
nt=int((tf-tn)/Dt)#Time discretitzation
x_array=np.linspace(x0,xf,Jx)#Values of x
t_array=np.linspace(tn,tf,nt)#Values of t
x_g_array,t_g_array=np.meshgrid(x_array,t_array)#Values of the x-t grid
A_n=sp.integrate.quad(lambda x: (2/L)*np.sin((np.pi/L)*x)*np.sin((np.pi/L)*x),0,L)#Constant of the analytical solution
def IC(x):#Initial condition, in this case, a sinusoidal function
    return np.sin((np.pi/L)*x)
def AnalyticalSol(t,x):#Analytical solution
    return A_n[0]*np.sin((np.pi/L)*x)*np.exp(-D*((np.pi/L)**2)*t)
def InitialCondition(x_array):#Preparation for the implementation of the initial condition
    u_array=np.zeros((nt,Jx))
    for jx in range(0,Jx):
        u_array[0,jx]=IC(x_array[jx])
    return u_array
def DirichletBC(u_array):#Dirichlet boundary conditions
    u_array[:,0]=0
    u_array[:,-1]=0
    return u_array
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Solution evolution
u_n_array=InitialCondition(x_array)#Numerical solution array
u_n_array=DirichletBC(u_n_array)
u_a_array=InitialCondition(x_array)#Analytical solution array
u_a_array=DirichletBC(u_a_array)
for n in range(1,nt):
    for i in range(1,Jx-1):
        u_a_array[n,i]=AnalyticalSol(t_array[n],x_array[i])
        u_n_array[n,i]=u_n_array[n-1,i]+D*(Dt/(Dx**2))*(u_n_array[n-1,i+1]+u_n_array[n-1,i-1]-2*u_n_array[n-1,i])
        u_n_array=DirichletBC(u_n_array)
Error_array=abs(u_a_array-u_n_array)
#Plots
#Numerical solution
figure,axes=plt.subplots()
graph=axes.contourf(x_g_array,t_g_array,u_n_array,levels=100,cmap="jet")
axes.set(title=r"Numerical solution of the 1D Diffusion equation",xlabel="$x$",ylabel="$t$")
plt.colorbar(graph)
plt.grid()
plt.tight_layout()
plt.show()
#Analytical solution
figure,axes=plt.subplots()
graph=axes.contourf(x_g_array,t_g_array,u_a_array,levels=100,cmap="jet")
axes.set(title=r"Analytical solution of the 1D Diffusion equation",xlabel="$x$",ylabel="$t$")
plt.colorbar(graph)
plt.grid()
plt.tight_layout()
plt.show()
print("Program execution time:",time.process_time()-start_time_program,"seconds.")
#Error between solutions
figure,axes=plt.subplots()
graph=axes.contourf(x_g_array,t_g_array,Error_array,levels=100,cmap="inferno")
axes.set(title=r"Error between solutions of the 1D Diffusion equation",xlabel="$x$",ylabel="$t$")
plt.colorbar(graph)
plt.grid()
plt.tight_layout()
plt.show()
print("Program execution time:",time.process_time()-start_time_program,"seconds.")