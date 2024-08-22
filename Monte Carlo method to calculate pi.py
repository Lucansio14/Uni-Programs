# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 13:23:36 2019
Calculate an aproximation of the number pi by using the Monte Carlo method.
@author: Lucas Romero Fernández
"""
import time
from random import random
import numpy as np
import matplotlib.pyplot as plt
#More information about the method in this website: https://en.wikipedia.org/wiki/Monte_Carlo_method
#main_program
#Definition of lists and constants
start_time_program=time.process_time()#To calculate the program execution time
xin=[]#X-coordinate of points inside the circle
yin=[]#Y-coordinate of points inside the circle
xout=[]#X-coordinate of points outside the circle
yout=[]#Y-coordinate of points ouside the circle
N=50000#Number of points. Increase it massively to gain precision if necessary. Must be an integer and bigger than 1.
#Main process of the obtainment of pi with the Monte Carlo Method
def MonteCarlo(N):
    for i in range(N):
        x=2*random()-1
        y=2*random()-1
        if (x**2+y**2<=1):
            xin.append(x)
            yin.append(y)
        else:
            xout.append(x)
            yout.append(y)
    pi=4*len(xin)/N
    return pi
MC_pi=MonteCarlo(N)
#Main process of the obtainment of pi with NumPy
numpy_pi=np.pi
E_t=abs(numpy_pi-MC_pi)#Error between methods
#Graph
fig=plt.figure(figsize=(9,5))
ax=plt.gca()
plt.xlim(-1,1)
plt.ylim(-1,1)
circle=plt.Circle((0,0),1,color="black",fill=False)
ax.add_patch(circle)
plt.scatter(xin,yin,color="g",s=1,alpha=0.5)
plt.scatter(xout,yout,color="r",s=1,alpha=0.5)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("Monte Carlo method to calculate pi ($N$ = {})".format(N))
ax.set_aspect("equal",adjustable="box")
plt.tight_layout()
plt.show()
#Results
print("Value of pi with the Monte Carlo method (N = {})".format(N),"=",MC_pi)#Result of the obtainment of pi with the Monte Carlo Method
print("")
print("Value of pi with NumPy =",numpy_pi)#Result of the obtainment of pi with NumPy
print("")
print("Value of the error =",E_t)#Result of the error between methods
print("")
print("Value of the absolute relative error =",abs(((E_t)/(numpy_pi)))*100,"%")#Result of the absolute relative error between methods
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")