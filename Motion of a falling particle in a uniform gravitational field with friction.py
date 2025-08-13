# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 11:42:43 2021
Representation of the motion (NURM or Non Uniform Rectilinear Motion) of a falling particle with mass m in a uniform gravitational field g
with friction force fr = -kv^2 (k corresponds to the coefficient of friction, an arbitrary parameter in this case), and at an initial height x0.
More information on this website: https://en.wikipedia.org/wiki/Linear_motion 
@author: Lucas Romero Fernández
"""
import time
import matplotlib.pyplot as plt
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants, and lists
m=0.00005#Units: kg. Has to be a positive value.
g=9.8#Units: m/s^2. Note: Positive value because the minus sign is already considered in the expression.
k=0.01*(10)**(-6)#Units: kg/m. Has to be a positive value.
x0=1500#Units: m. Has to be a positive value.
dt=0.01#Time-step size. Units: s. Has to be a positive value. Decrease its value to gain accuracy.
x_list=[]
v_list=[]
t_list=[]
i=0#Counter for the amount of time-steps.
#Main process of the obtainment of the trajectory using the linear motion expressions.
while True:
    if i==0:
        x=x0
        v=0
        x_list.append(x)
        v_list.append(v)
        t_list.append(0)
        i+=1
    else:
        x=x_list[-1]+v_list[-1]*dt+0.5*(-g+(k/m)*(v_list[-1])**2)*dt**2
        v=v_list[-1]+(-g+(k/m)*(v_list[-1])**2)*dt
        if x_list[-1]>=0:
            x_list.append(x)
            v_list.append(v)
            t_list.append(i*dt)
            i+=1
        else:
            break
#Graphs
#x vs t
plt.figure(figsize=(9,5))
plt.plot(t_list,x_list)
plt.axhline(0,color="black",linestyle="dashed",alpha=0.7)
plt.axhline(x0,color="black",linestyle="dashed",alpha=0.7)
plt.xlabel("t [s]")
plt.ylabel("x [m]")
plt.title("Displacement of the falling particle")
plt.grid()
plt.tight_layout()
plt.show()
#v vs t
plt.figure(figsize=(9,5))
plt.plot(t_list,v_list)
plt.xlabel("t [s]")
plt.ylabel("v [m/s]")
plt.title("Falling velocity of the particle")
plt.grid()
plt.tight_layout()
plt.show()
print("Approximated estimation of the particle falling time:",t_list[-1],"seconds.")
print("")
print("Maximum absolute falling speed achieved:",v_list[-1],"m/s.")
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")