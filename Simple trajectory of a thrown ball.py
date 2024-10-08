# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 10:21:15 2021
Representation of a simple/ideal trajectory of a ball/proyectile thrown with a certain initial speed v0 and an angle theta, with respect to the horizontal axis,
and at an initial height y0. The expression of this type of trayectory is: y(x) = f(x) = x * tan(theta) − (g * x^2 / (2 * v0^2 * (cos(theta))^2)) + y0
More information on this website:https://en.wikipedia.org/wiki/Projectile_motion 
@author: Lucas Romero Fernández
"""
import time
import matplotlib.pyplot as plt
import numpy as np
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants, functions and lists
v0=25/6#Units: m/s. Has to be a positive value to obtain a trajectory for positive values of x.
theta_deg=60#Units: degrees. Has to be a value in the interval [0,90] to obtain a trajectory for positive values of x.
y0=100#Units: m. Has to be a positive value.
g=9.8#Units: m/s. Note: Positive value because the minus sign is already considered in the expression.
n=1000#Number of points in the x-axis
theta_rad=np.radians(theta_deg)
x_max=(v0**2*(np.cos(theta_rad))**2*(np.tan(theta_rad)+np.sqrt((np.tan(theta_rad))**2+(2*g*y0)/(v0**2*(np.cos(theta_rad))**2))))/g#End of the trajectory (y = 0 for positive values of x)
x_list=np.linspace(0,x_max,n)
y_list=[]
#Main process of the obtainment of the trajectory using the expression already shown
for x in x_list:
    y=x*np.tan(theta_rad)-((g*(x)**(2))/(2*(v0)**(2)*(np.cos(theta_rad))**(2)))+y0
    y_list.append(y)
#Graph
plt.figure(figsize=(9,5))
plt.plot(x_list,y_list)
plt.axhline(0,color="black",linestyle="dashed",alpha=0.7)
plt.axvline(0,color="black",linestyle="dashed",alpha=0.7)
plt.xlabel("x [m]")
plt.ylabel("y(x) [m]")
plt.title("Ideal trajectory of the thrown ball")
plt.grid()
plt.tight_layout()
plt.show()
print("Distance traveled by the ball:",x_max,"meters.")
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")