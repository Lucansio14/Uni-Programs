# -*- coding: utf-8 -*-
"""
Created on Wed Oct 6 12:07:15 2021
Solver of 1st Order Ordinary Differential Equations (ODEs): y'(x) = f(x,y) using the Runge-Kutta-Fehlberg method of 4th order with an error estimator
of 5th order (RKF45) with an adaptive (or not) step size "h" in the region [x0,xf] in the x-axis with initial condition y0 = y(x0).
More information on this website: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method.
@author: Lucas Romero Fernández
"""
import time
import matplotlib.pyplot as plt
import numpy as np
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants, functions and lists
h0=0.1#Initial (or fixed) value of the step size
h=h0
max_error=10**(-6)#Maximum error to be considered in the equation's solution
x0_init=1
x0=x0_init
y0=1
xf=10
h_list=[h]
x_list=[x0]
y_list=[y0]
y_prime_list=[y0]
error_appr_list=[0]#Error between RKF4 and RKF5 approximations
analytical_sol_list=[y0]
error_sol_list=[0]#Error between RKF5 and analytical solutions
def f(x,y):
    return ((-1)/(x**2))-4*(x-6)*np.exp(-2*(x-6)**2)
def analytical_sol(x):#Specific analytical solution for the case of f(x,y)=((-1)/(x^2))-4*(x-6)*e^(-2*(x-6)^2), for comparison purposes.
    return (y0+((1/x)+np.exp(-2*(x-6)**(2)))-((1/x0)+np.exp(-2*(x0-6)**2)))
flag=float(input("To utilize an adaptive step size, type '1'; for a fixed step size, type '0': "))
if flag!=0 and flag!=1:
    flag_error="Please, type '1' or '0' to decide between an adaptive or a fixed step size, respectively."
    raise Exception(flag_error)
#Solution evolution
while x0<=xf:
    while True:
        #RKF45 subroutine
        f0=f(x0,y0)
        f1=f(x0+h/4,y0+h/4*f0)
        f2=f(x0+3/8*h,y0+3*h/32*f0+9*h/32*f1)
        f3=f(x0+12*h/13,y0+1932*h/2197*f0-7200*h/2197*f1+7296*h/2197*f2)
        f4=f(x0+h,y0+439*h/216*f0-8*h*f1+3680*h/513*f2-845*h/4104*f3)
        f5=f(x0+h/2,y0-8*h/27*f0+2*h*f1-3544*h/2565*f2+1859*h/4104*f3-11*h/40*f4)
        y=y0+h*(25/216*f0+1408/2565*f2+2197/4104*f3-1/5*f4)
        y_prime=y0+h*(16/135*f0+6656/12825*f2+28561/56430*f3-9/50*f4+2/55*f5)
        anal_sol=analytical_sol(x0+h)
        error_appr=np.abs(y-y_prime)
        error_sol=np.abs(anal_sol-y_prime)
        if flag==1:#Adaptive step size case
            h_new=0.9*h*(h*max_error/error_appr)**(1/4)
            if h_new<=h:
                h=h_new
            if h_new>h:
                y_list.append(y)
                y_prime_list.append(y_prime)
                error_appr_list.append(error_appr)
                analytical_sol_list.append(anal_sol)
                error_sol_list.append(error_sol)
                h_list.append(h_new)
                x_list.append(x0)
                x0+=h
                y0=y
                h=h_new
                break
        if flag==0:#Fixed step size case
            y_list.append(y)
            y_prime_list.append(y_prime)
            error_appr_list.append(error_appr)
            analytical_sol_list.append(anal_sol)
            error_sol_list.append(error_sol)
            h_list.append(h)
            x_list.append(x0)
            x0+=h
            y0=y
            break
#Graphs
#Comparison between solutions
plt.figure(figsize=(9,5))
plt.xlim(x0_init-h0,xf+h0)
plt.plot(x_list,y_list,c="red",label="RKF4")
plt.plot(x_list,y_prime_list,c="green",label="RKF5")
plt.scatter(x_list,analytical_sol_list,c="blue",s=30,label="Analytical solution")
plt.xlabel("$x$")
plt.ylabel("$y(x)$")
plt.title("Comparison between the analytical solution and the numerical RKF45 solutions")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#Errors of the RKF45 method
plt.figure(figsize=(9,5))
plt.xlim(x0_init-h0,xf+h0)
plt.plot(x_list,error_appr_list,"-ro",label="Error between RKF4 and RKF5 solutions")
plt.plot(x_list,error_sol_list,"-go",label="Error between RKF5 and analytical solutions")
plt.xlabel("$x$")
plt.ylabel("$err[y(x)]$")
plt.title("Errors of the RKF45 method")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#Evolution of the step size 'h'
plt.figure(figsize=(9,5))
plt.xlim(x0_init-h0,xf+h0)
plt.plot(x_list,h_list,"-ro")
plt.xlabel("$x$")
plt.ylabel("$h$")
plt.title("Evolution of the step size")
plt.grid()
plt.tight_layout()
plt.show()
print("Program execution time:",time.process_time()-start_time_program,"seconds.")
