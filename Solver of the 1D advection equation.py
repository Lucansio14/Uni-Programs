# -*- coding: utf-8 -*-
"""
Created on Tue Mar 4 17:32:14 2025
Solver of the 1D advection equation (\partial_{t}u = -a\partial_{x}u) with RK4 (time part) and FD4 (spatial part).
More information about the equation in this website: https://farside.ph.utexas.edu/teaching/329/lectures/node90.html
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
#Definition of general constants, variables, arrays and functions
a=1#Flow speed (Units = cm/s)
dx=0.02#Space-step size (Units = cm) (100 points = 0.02, 200 points = 0.01, 400 points = 0.005, for CFL coefficient of 0.8)
dt=0.8*(dx/a)#Time-step size (Units = s). For stability purposes, we use the CFL criteria (lambda=((a*dt)/dx)<=1) to obtain this value.
x0=-1#Start of the space domain (Units = cm)
xf=1#End of the space domain (Units = cm)
t0=0#Start time (Units = s)
tf=100#End time (Units = s)
dPlot=125#For the plots (specifically, the number of iterations made between plots) (100 points = 125, 200 points = 250 and 400 points = 500, for CFL coefficient of 0.8)
J=int((xf-x0)/dx)#Number of spatial points
N=int((tf-t0)/dt)#Time discretitzation
x_array=np.linspace(x0,xf,J)#Values of x
dudr_array=np.zeros(J)#Finite differences spatial derivative scheme array
animation_frames_list=[]#Animation frames/artists for the results
Num_frame=0#Number of solution frame for the animation
def IC(x):#Initial condition, in this case, a Gaussian function
    return np.exp(-(1/2)*((x-0)/(0.1))**2)
def InitialCondition(x_array):#Preparation for the implementation of the initial condition
    j=0
    u_array=np.zeros(J)
    for x in x_array:
        u_array[j]=IC(x)
        j+=1
    return u_array
def Finite_Difference_Derivative_Order4(u_array):#Spatial discretization: Finite/discrete central differences (4rth order)
    for j in range(2,J-2):
        dudr_array[j]=(1/(12*dx))*(-u_array[j+2]+8*u_array[j+1]-8*u_array[j-1]+u_array[j-2])
    return dudr_array
def RHS(u_array,a):#Right Hand Side of the advection equation
    return -a*Finite_Difference_Derivative_Order4(u_array)
def RK4(u_array,a):#Temporal discretization: Runge-Kutta of fourth order (RK4)
    k_1=dt*RHS(u_array,a)
    u_i_array=u_array+0.5*k_1
    u_i_array=Periodic_BC(u_i_array)
    k_2=dt*RHS(u_i_array,a)
    u_i_array=u_array+0.5*k_2
    u_i_array=Periodic_BC(u_i_array)
    k_3=dt*RHS(u_i_array,a)
    u_i_array=u_array+k_3
    u_i_array=Periodic_BC(u_i_array)
    k_4=dt*RHS(u_i_array,a)
    dudt_array=u_array+((1)/(6))*(k_1+2*k_2+2*k_3+k_4)
    return dudt_array
def Periodic_BC(u_array):#(circular) Periodic boundary conditions
    u_array[0]=u_array[-4]
    u_array[1]=u_array[-3]
    u_array[-2]=u_array[2]
    u_array[-1]=u_array[3]
    return u_array
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Solution evolution and plotting of the results
u_array=InitialCondition(x_array)
fig,ax=plt.subplots(figsize=(9,5))
plt.xlim(x0,xf)
plt.title("Solution of the 1D Advection equation")
plt.xlabel(r"$x$")
plt.ylabel(r"$u(x,t)$")
plt.grid()
plt.tight_layout()
for n in range(N+1):
    if n % dPlot==0:
        frame=ax.plot(x_array,u_array,color="blue")
        print("Frame ",Num_frame," computed at t = ",str(n*dt)," s.",sep="")
        animation_frames_list.append(frame)
        Num_frame+=1
    u_new_array=RK4(u_array,a)#Time evolution
    u_new_array=Periodic_BC(u_new_array)#Imposing BCs
    u_array=u_new_array#Update for the next time step
print("All frames computed, creating animation...")
animation=animation.ArtistAnimation(fig=fig,artists=animation_frames_list,interval=200)#Animation of the results
animation.save(filename="1D Advection Solution.mp4",writer="ffmpeg")#Animation save
print("Animation saved at ",os.getcwd(),".",sep="")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")