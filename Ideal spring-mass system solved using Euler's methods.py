# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 20:38:11 2021
Context: An ideal spring-mass system of angular frecuency omega and phase phi, consisting of a mass m suspended from an 
ideal massless spring with a proporcionality constant K. This system will be solved (aka obtain the position and velocity 
with respect to the time, x(t) and v(t), respectively) using three different Euler's methods with the value of x(0) = x0 
= amplitude (A), v(0) = v0 and dt as chosen fixed constants.
More information of this system in this website: https://phys.libretexts.org/Courses/Berea_College/Introductory_Physics%3A_Berea_College/13%3A_Simple_Harmonic_Motion/13.01%3A_The_motion_of_a_spring-mass_system
@author: Lucas Romero Fernández
"""
import time
import matplotlib.pyplot as plt
import numpy as np
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants, functions and lists
m=1#Units = kg
K=1#Units = N/m
omega=np.sqrt(K/m)
phi=0#Note: It is 0 only if it is considered that x(0) = A
t0=0#Units = s
tf=10#End time (Units = s)
x0=1#Units = m
A=x0#Units = m
v0=0#Units = m/s
E0=1/2*m*v0**2+1/2*K*x0**2#Initial energy of the system (Units = J)
dt=0.1#Time step (Units = s)
t_list=[t0]
x1_list,x2_list,x3_list=[x0],[x0],[x0]
v1_list,v2_list,v3_list,v1c_list,v2c_list,v3c_list=[v0],[v0],[v0],[v0],[v0],[v0]
E1_list,E2_list,E3_list=[E0],[E0],[E0]
analytical_sol_list_x=[]
analytical_sol_list_v=[]
analytical_sol_list_E=[]
analytical_sol_list_vc=[]
def analytical_sol1(t):#x(t)
    return A*np.cos(omega*t+phi)
def analytical_sol2(t):#v(t)
    return -A*omega*np.sin(omega*t+phi)
def f1(t,x,v):#Movement equation 1 (dx/dt)
    return v
def f2(t,x,v):#Movement equation 2 (dv/dt)
    return (-(K/m)*x)
#Simple Euler method
while t0<tf:
    x1=x0+dt*f1(t0,x0,v0)
    v1=v0+dt*f2(t0,x0,v0)
    E1=1/2*m*v1**2+1/2*K*x1**2
    v1c=np.sqrt(np.abs((K*(A**2-x0**2))/(m)))
    x1_list.append(x1)
    v1_list.append(v1)
    E1_list.append(E1)
    v1c_list.append(v1c)
    t0=np.float64(t0+dt)
    t_list.append(t0)
    x0=x1
    v0=v1
#Modified Euler method
t0=0
x0=1
v0=0
t_list=[t0]
while t0<tf:
    t_mid=t0+dt/2
    x_mid=x0+(dt/2)*f1(t0,x0,v0)
    v_mid=v0+(dt/2)*f2(t0,x0,v0)
    x2=x0+dt*f1(t_mid,x_mid,v_mid)
    v2=v0+dt*f2(t_mid,x_mid,v_mid)
    E2=1/2*m*v2**2+1/2*K*x2**2
    v2c=np.sqrt(np.abs((K*(A**2-x0**2))/(m)))
    x2_list.append(x2)
    v2_list.append(v2)
    E2_list.append(E2)
    v2c_list.append(v2c)
    t0=np.float64(t0+dt)
    t_list.append(t0)
    x0=x2
    v0=v2
#Improved Euler method
t0=0
x0=1
v0=0
t_list=[t0]
while t0<tf:
    x3=x0+(dt/2)*(f1(t0,x0,v0)+f1(t0+dt,x0+dt*f1(t0,x0,v0),v0+dt*f1(t0,x0,v0)))
    v3=v0+(dt/2)*(f2(t0,x0,v0)+f2(t0+dt,x0+dt*f2(t0,x0,v0),v0+dt*f2(t0,x0,v0)))
    E3=1/2*m*v3**2+1/2*K*x3**2
    v3c=np.sqrt(np.abs((K*(A**2-x0**2))/(m)))
    x3_list.append(x3)
    v3_list.append(v3)
    E3_list.append(E3)
    v3c_list.append(v3c)
    t0=np.float64(t0+dt)
    t_list.append(t0)
    x0=x3
    v0=v3
#Analytical method
t_list=np.array(t_list)
for i in range(0,len(t_list)):
    analytical_sol_list_x.append(analytical_sol1(t_list[i]))
    analytical_sol_list_v.append(analytical_sol2(t_list[i]))
    analytical_sol_list_E.append(1/2*m*analytical_sol_list_v[i]**2+1/2*K*analytical_sol_list_x[i]**2)
    analytical_sol_list_vc.append(np.sqrt(np.abs((K*(A**2-analytical_sol_list_x[i]**2))/(m))))
#Graphs
#x(t) vs t
plt.figure(figsize=(9,5))
plt.plot(t_list,x1_list,c="red",label="Simple Euler method")
plt.plot(t_list,x2_list,c="green",label="Modified Euler method",linestyle="dashed")
plt.plot(t_list,x3_list,c="purple",label="Improved Euler method",linestyle="dotted")
plt.plot(t_list,analytical_sol_list_x,c="black",label="Analytical solution",linestyle="dashdot")
plt.xlabel("t [s]")
plt.ylabel("x(t) [m]")
plt.title("Position comparison of the methods")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#v(t) vs t (obtained through equations)
plt.figure(figsize=(9,5))
plt.plot(t_list,v1_list,c="red",label="Simple Euler method")
plt.plot(t_list,v2_list,c="green",label="Modified Euler method",linestyle="dashed")
plt.plot(t_list,v3_list,c="purple",label="Improved Euler method",linestyle="dotted")
plt.plot(t_list,analytical_sol_list_v,c="black",label="Analytical solution",linestyle="dashdot")
plt.xlabel("t [s]")
plt.ylabel("v(t) [m/s]")
plt.title("Velocity comparison of the methods (obtained through equations)")
plt.legend()
plt.show()
#E vs t
plt.figure(figsize=(9,5))
plt.plot(t_list,E1_list,c="red",label="Simple Euler method")
plt.plot(t_list,E2_list,c="green",label="Modified Euler method",linestyle="dashed")
plt.plot(t_list,E3_list,c="purple",label="Improved Euler method",linestyle="dotted")
plt.plot(t_list,analytical_sol_list_E,c="black",label="Analytical solution",linestyle="dashdot")
plt.xlabel("t [s]")
plt.ylabel("E [J]")
plt.title("Energy comparison of the methods")
plt.legend()
plt.show()
#v(t) vs t (obtained through energy conservation)
plt.figure(figsize=(9,5))
plt.plot(t_list,v1c_list,c="red",label="Simple Euler method")
plt.plot(t_list,v2c_list,c="green",label="Modified Euler method",linestyle="dashed")
plt.plot(t_list,v3c_list,c="purple",label="Improved Euler method",linestyle="dotted")
plt.plot(t_list,analytical_sol_list_vc,c="black",label="Analytical solution",linestyle="dashdot")
plt.xlabel("t [s]")
plt.ylabel("abs(v(t)) [m/s]")#Absolute value because of square root and negative values of v in the equation with energy conservation 
plt.title("Velocity comparison of the methods (obtained through energy conservation)")
plt.legend()
plt.show()
print("Period of the motion (T):",2*np.pi/omega,"seconds.")
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")