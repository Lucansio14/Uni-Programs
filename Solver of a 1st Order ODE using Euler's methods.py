# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 20:20:34 2021
Solver of an Ordinary Differential Equation (ODE) of first order f(x): y'(x) = y^2 + 1 using Euler's different methods with different step sizes "h" in the region [x0,xf] in the x-axis with initial condition y0=y(0).
More information on this websites: https://en.wikipedia.org/wiki/Euler_method#:~:text=The%20Euler%20method%20is%20a,proportional%20to%20the%20step%20size.
                                   https://www.maths.tcd.ie/~ryan/teaching/11404/March18_20_HigherOrderMethods.pdf
                                   https://math.libretexts.org/Courses/Monroe_Community_College/MTH_225_Differential_Equations/03%3A_Numerical_Methods/3.02%3A_The_Improved_Euler_Method_and_Related_Methods
@author: Lucas Romero Fernández
"""
import time
import matplotlib.pyplot as plt
import numpy as np
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants, functions and lists
h_list=[0.05,0.10,0.15,0.20]#Different step sizes used for the Euler's methods
n=0#Number to "count" all the cases with different step sizes
for h in h_list:
    x0=0
    y0=0
    xf=1
    x_list=[x0]
    y1_list=[y0]
    y2_list=[y0]
    y3_list=[y0]
    y4_list=[y0]
    analytical_sol_list=[]
    def analytical_sol(x):
        x0=0
        y0=0
        return (np.tan(np.arctan(y0)+x-x0))
    def f(x,y):
        return ((y)**2+1)
    #Simple Euler method of first order
    while x0<xf:
        y1=y0+h*f(x0,y0)
        y1_list.append(y1)
        x0=np.float32(x0)+h
        x_list.append(x0)
        y0=y1
    #Simple Euler method with Taylor development up to 2nd order
    x0=0
    y0=0
    x_list=[x0]
    while x0<xf:
        y2=y0+h*f(x0,y0)+((h)**2/(2))*f(x0,y0)*2*y0
        y2_list.append(y2)
        x0=np.float32(x0)+h
        x_list.append(x0)
        y0=y2
    #Modified Euler method
    x0=0
    y0=0
    x_list=[x0]
    while x0<xf:
        x_mid=x0+h/2
        y_mid=y0+(h/2)*f(x0,y0)
        y3=y0+h*f(x_mid,y_mid)
        y3_list.append(y3)
        x0=np.float32(x0)+h
        x_list.append(x0)
        y0=y3
    #Improved Euler method
    x0=0
    y0=0
    x_list=[x0]
    while x0<xf:
        y4=y0+(h/2)*(f(x0,y0)+f(x0+h,y0+h*f(x0,y0)))
        y4_list.append(y4)
        x0=np.float32(x0)+h
        x_list.append(x0)
        y0=y4
    #Analytical solution
    x_list=np.array(x_list)
    for i in range(0,len(x_list)):
        analytical_sol_list.append(analytical_sol(x_list[i]))
    #Graphs
    plt.figure(figsize=(9,5))
    plt.xlim(0,1)
    plt.plot(x_list,analytical_sol_list,c='blue',label="Analytical solution")
    plt.plot(x_list,y1_list,c='red',label="Simple Euler method of 1st order")
    plt.plot(x_list,y2_list,c='orange',label="Simple Euler method with Taylor development up to 2nd order")
    plt.plot(x_list,y3_list,c='green',label="Modified Euler method")
    plt.plot(x_list,y4_list,c='purple',label="Improved Euler method")
    plt.xlabel("x")
    plt.ylabel("y(x)")
    plt.title("Comparison of the analytical solution and different Euler's methods with $h=$"+str(h_list[n]))
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()
    n+=1
print("Program execution time:",time.process_time()-start_time_program,"seconds.")