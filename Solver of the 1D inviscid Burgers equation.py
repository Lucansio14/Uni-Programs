# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 10:47:56 2025
Solver of the 1D inviscid Burgers' equation: The most basic non linear equation (\partial_{t}u = -u\partial_{x}u) with RK4 (time part), FD4 (spatial part) and Kreiss-Oliger dissipation filter.
More information about the equation in this website: https://en.wikipedia.org/wiki/Burgers%27_equation
More information about CFL criteria for non linear equations in this website: https://scicomp.stackexchange.com/questions/31456/cfl-equation-for-non-linear-equation

Context about the computed solution:
Without the use of more sophisticated methods (like shock capturing), the solution is "corrupted" in some way by a spike that increases considerably in time, destroying it entirely in little time.
The reason of this issue is completely numerical. It can be mitigated somewhat by relying heavily in the the Kreiss-Oliger dissipation filter, small space-step sizes and restrictive CFL conditions.
However, this is not enough for a solution of a reasonable quality. Furthermore, this also produces problems for the CFL condition, and the "solutions/graphs" timers.
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
#Definition of general constants, variables, arrays and functions
dx=0.01#Space-step size (Units = cm)
x0=-1#Start of the space domain (Units = cm)
xf=1#End of the space domain (Units = cm)
t0=0#Start time (Units = s)
tf=1#End time (Units = s)
dPlot=1#For the plots (specifically, the number of iterations made between plots)
sigma=1.5#Control parameter of the Kreiss-Oliger dissipation filter term (Qd) (sigma >= 0 necessary!)
Amp_u=1#Amplitude of the initial gaussian perturbation
J=int((xf-x0)/dx)#Number of spatial points
x_array=np.linspace(x0,xf,J)#Values of x
Qd_array=np.zeros(J)#Kreiss-Oliger dissipation filter term array
RHS_array=np.zeros(J)#Right Hand Side array
animation_frames_list=[]#Animation frames/artists for the results
Num_frame=0#Number of solution frame for the animation
def IC(x):#Initial condition, in this case, a Gaussian function
    return Amp_u*np.exp(-(1/2)*((x-0)/(0.1))**2)
def InitialCondition(x_array):#Preparation for the implementation of the initial condition
    j=0
    u_array=np.zeros(J)
    for x in x_array:
        u_array[j]=IC(x)
        j+=1
    return u_array
def Finite_Difference_Derivative_Order4(u_array,j):#Spatial discretization: Finite/discrete central differences (4rth order)
    dudr=(1/(12*dx))*(-u_array[j+2]+8*u_array[j+1]-8*u_array[j-1]+u_array[j-2])
    return dudr
def RHS(u_array):#Right Hand Side of the inviscid Burgers' equation:
    for j in range(2,J-2):
        if j==2 or j==J-3:
            Qd_array[j]=0
        else:
            Qd_array[j]=(1/(64*dx))*(u_array[j+3]-6*u_array[j+2]+15*u_array[j+1]-20*u_array[j]+15*u_array[j-1]-6*u_array[j-2]+u_array[j-3])
        RHS_array[j]=-u_array[j]*Finite_Difference_Derivative_Order4(u_array,j)+sigma*Qd_array[j]
    return RHS_array
def RK4(u_array):#Temporal discretization: Runge-Kutta of fourth order (RK4)
    k_1=dt*RHS(u_array)
    u_i_array=u_array+0.5*k_1
    u_i_array=Periodic_BC(u_i_array)#Imposing periodic BCs
    k_2=dt*RHS(u_i_array)
    u_i_array=u_array+0.5*k_2
    u_i_array=Periodic_BC(u_i_array)#Imposing periodic BCs
    k_3=dt*RHS(u_i_array)
    u_i_array=u_array+k_3
    u_i_array=Periodic_BC(u_i_array)#Imposing periodic BCs
    k_4=dt*RHS(u_i_array)
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
dt=0.8*(dx/max(u_array))#Variable time-step size (Units = s). For stability purposes, we use the CFL criteria (lambda=((max(u)*dt)/dx)<=1) to obtain this value.
N=int((tf-t0)/dt)#Time discretitzation
fig,ax=plt.subplots(figsize=(9,5))
plt.xlim(x0,xf)
plt.ylim(0-0.2*Amp_u,Amp_u+0.6*Amp_u)
plt.title("Solution of the 1D inviscid Burgers equation")
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
    u_new_array=RK4(u_array)#Time evolution
    u_new_array=Periodic_BC(u_new_array)#Imposing BCs
    u_array=u_new_array#Update for the next time step
    dt=0.8*(dx/abs(max(u_array)))#Variable time-step size (Units = s). For stability purposes, we use the CFL criteria (lambda=((max(u)*dt)/dx)<=1) to obtain this value.
    N=int((tf-t0)/dt)#Time discretitzation
print("All frames computed, creating animation...")
animation=animation.ArtistAnimation(fig=fig,artists=animation_frames_list,interval=200)#Animation of the results
animation.save(filename="1D Inviscid Burgers Solution.mp4",writer="ffmpeg")#Animation save
print("Animation saved at ",os.getcwd(),".",sep="")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")