# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 11:02:55 2021
Solver of the 1D parabolic diffusion Parcial Derivative Equation (PDE) or the 1D heat equation (with diffusivity constant equal to 1): u_t = u_xx using the theta
method and Thomas algorithm in the region [x0,xf] in the x-axis with boundary conditions (BCs) u(t,x0) = u(t,xf) = 0 and initial condition u(t0,x) = u0(x).
More information on these websites: https://en.wikipedia.org/wiki/Heat_equation
                                    https://leifh.folk.ntnu.no/teaching/tkt4140/._main065.html
                                    https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants, functions, arrays and lists
Num_Spatial_points=41#Number of numerical spatial grid points
Time_Iter_list=[0,1,10,100,1000,2000,3000,4000]#Each time step (dt) iteration represented in the plots
x0=0
xf=1
x_array=np.linspace(x0,xf,Num_Spatial_points)#Spatial grid
dx=xf/len(x_array)
t0=0
tf=0.5
nu=0.2#Ratio between spatial and time grid discretizations (dx and dt, respectively)
dt=nu*dx**2
t_array=np.arange(t0,tf+dt,dt)#Time grid
theta=0.75#Parameter that performs a weighted average of two other numerical schemes (the FTCS (theta = 0) and BTCS (theta = 1) methods). An special case is obtained when theta = 0.5 (Crank Nicholson scheme).
a_array=np.zeros(Num_Spatial_points)#Vector of the Thomas algorithm
b_array=np.zeros(Num_Spatial_points)#Vector of the Thomas algorithm
c_array=np.zeros(Num_Spatial_points)#Vector of the Thomas algorithm
d_array=np.zeros(Num_Spatial_points)#Vector of the Thomas algorithm
for i in range(Num_Spatial_points):
    a_array[i]=-theta*nu
    b_array[i]=1+2*theta*nu
    c_array[i]=-theta*nu
def InitialCondition(x_array):#Preparation for the implementation of the initial condition
    f_array=np.zeros(len(x_array))
    n=0
    for i in x_array:
        if (np.abs(i-1/2))<=1/4:
            f_array[n]=1
        else:
            f_array[n]=0
        n+=1
    return f_array
def analytical_sol(t_array):#Specific analytical solution for the heat equation for these BCs and initial condition, for comparison purposes.
    x=0
    Num_Analytical_Spatial_points=1000#Number of analytical spatial grid points
    dx=1/Num_Analytical_Spatial_points
    x_list=[]
    u_list=[]
    while x<=xf:
        m=1
        sum_u=0
        while m<=Num_Analytical_Spatial_points:
            u=2/(np.pi*m)*(np.cos(m*np.pi/4)-np.cos(m*3/4*np.pi))*np.exp(-(m*np.pi)**2*t_array)*np.sin(m*np.pi*x)
            sum_u+=u
            m+=1
        u_list.append(sum_u)
        x_list.append(x)
        x+=dx
    return (x_list,u_list)
def Thomas(a_array,b_array,c_array,d_array):#Thomas or tridiagonal_matrix algorithm
    e_array=np.zeros(Num_Spatial_points)
    f_array=np.zeros(Num_Spatial_points)
    e_array[0]=0#Due to the specific BCs
    f_array[0]=0#Due to the specific BCs
    for i in range(1,Num_Spatial_points-1):
        e_array[i]=(-c_array[i]/(b_array[i]+a_array[i]*e_array[i-1]))
        f_array[i]=(d_array[i]-a_array[i]*f_array[i-1])/(b_array[i]+a_array[i]*e_array[i-1])
    x_array=np.zeros(Num_Spatial_points)
    x_array[-1]=f_array[-1]
    for j in range(len(x_array)-2,-1,-1):
        x_array[j]=f_array[j]+e_array[j]*x_array[j+1]
    return x_array
#Solution evolution
Theta_method=lambda u_array,i:nu*(1-theta)*u_array[i+1]+(1+2*nu*(theta-1))*u_array[i]+nu*(1-theta)*u_array[i-1]#Theta method for the heat equation
u_array=[0]*(len(t_array))
u0_array=InitialCondition(x_array)
u_array[0]=u0_array
m=1
while t0<(tf-dt):
    for i in range(1,len(u0_array)-1):
        d_array[i]=Theta_method(u0_array,i)
    u0_array=Thomas(a_array,b_array,c_array,d_array)
    u_array[m]=u0_array
    t0+=dt
    m+=1
#Graphs for each time step iterations
for time_iter in Time_Iter_list:
    plt.figure(figsize=(9,5))
    plt.xlim(x0-0.025,xf+0.025)
    plt.ylim(min(u_array[0])-0.1,max(u_array[0])+0.1)
    plt.plot(analytical_sol(time_iter*dt)[0],analytical_sol(time_iter*dt)[1],label='Analytical')
    plt.plot(x_array,u_array[time_iter],label='Numerical')
    plt.xlabel("$x$")
    plt.ylabel("$u(x,t)$")
    plt.title(f'Comparison between analytical and numerical solutions for $t={time_iter}dt$')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

print("Program execution time:",time.process_time()-start_time_program,"seconds.")
