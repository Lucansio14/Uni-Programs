# -*- coding: utf-8 -*-
"""
Created on Thu Oct 7 11:23:54 2021
Solver of 2nd Order Ordinary Differential Equations (ODEs): x'(t,v) = v(t,x) = f_1(t,x,v) and v'(t,x,v) = f_2(t,x,v) using the Runge-Kutta-Fehlberg method of 4th order with
an error estimator of 5th order (RKF45) with an adaptive (or not) step size "h" in the region [t0,tf] in the t-axis with initial conditions x0 = x(t0) and v0 = v(t0,x0).
More information on these websites: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
                                    https://en.wikipedia.org/wiki/Van_der_Pol_oscillator
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
t0_init=0
t0=t0_init
x0=2
v0=0
tf=30
mu=4#Scalar parameter for the specific ODE of the Van der Pol oscillator (it indicates the nonlinearity and the strength of the damping)
h_list=[h]
t_list=[t0]
x_list=[x0]
x_prime_list=[x0]
v_list=[v0]
v_prime_list=[v0]
x_error_appr_list=[0]#Error of x between RKF4 and RKF5 approximations
v_error_appr_list=[0]#Error of v between RKF4 and RKF5 approximations
def f_1(t,x,v):
    return v
def f_2(t,x,v):#This specific case corresponds to the ODE of the Van der Pol oscillator
    return mu*(1-x**2)*v-x
flag=float(input("To utilize an adaptive step size, type '1'; for a fixed step size, type '0': "))
if flag!=0 and flag!=1:
    flag_error="Please, type '1' or '0' to decide between an adaptive or a fixed step size, respectively."
    raise Exception(flag_error)
#Solution evolution
while t0<=tf:
    while True:
        #RKF45 subroutine
        f0x=f_1(t0,x0,v0)
        f0v=f_2(t0,x0,v0)
        f1x=f_1(t0+h/4,x0+h/4*f0x,v0+h/4*f0v)
        f1v=f_2(t0+h/4,x0+h/4*f0x,v0+h/4*f0v)
        f2x=f_1(t0+3/8*h,x0+3*h/32*f0x+9*h/32*f1x,v0+3*h/32*f0v+9*h/32*f1v)
        f2v=f_2(t0+3/8*h,x0+3*h/32*f0x+9*h/32*f1x,v0+3*h/32*f0v+9*h/32*f1v)
        f3x=f_1(t0+12*h/13,x0+1932*h/2197*f0x-7200*h/2197*f1x+7296*h/2197*f2x,v0+1932*h/2197*f0v-7200*h/2197*f1v+7296*h/2197*f2v)
        f3v=f_2(t0+12*h/13,x0+1932*h/2197*f0x-7200*h/2197*f1x+7296*h/2197*f2x,v0+1932*h/2197*f0v-7200*h/2197*f1v+7296*h/2197*f2v)
        f4x=f_1(t0+h,x0+439*h/216*f0x-8*h*f1x+3680*h/513*f2x-845*h/4104*f3x,v0+439*h/216*f0v-8*h*f1v+3680*h/513*f2v-845*h/4104*f3v)
        f4v=f_2(t0+h,x0+439*h/216*f0x-8*h*f1x+3680*h/513*f2x-845*h/4104*f3x,v0+439*h/216*f0v-8*h*f1v+3680*h/513*f2v-845*h/4104*f3v)
        f5x=f_1(t0+h/2,x0-8*h/27*f0x+2*h*f1x-3544*h/2565*f2x+1859*h/4104*f3x-11*h/40*f4x,v0-8*h/27*f0v+2*h*f1v-3544*h/2565*f2v+1859*h/4104*f3v-11*h/40*f4v)
        f5v=f_2(t0+h/2,x0-8*h/27*f0x+2*h*f1x-3544*h/2565*f2x+1859*h/4104*f3x-11*h/40*f4x,v0-8*h/27*f0v+2*h*f1v-3544*h/2565*f2v+1859*h/4104*f3v-11*h/40*f4v)
        x=x0+h*(25/216*f0x+1408/2565*f2x+2197/4104*f3x-1/5*f4x)
        x_prime=x0+h*(16/135*f0x+6656/12825*f2x+28561/56430*f3x-9/50*f4x+2/55*f5x)
        x_error_appr=np.abs(x-x_prime)
        v=v0+h*(25/216*f0v+1408/2565*f2v+2197/4104*f3v-1/5*f4v)
        v_prime=v0+h*(16/135*f0v+6656/12825*f2v+28561/56430*f3v-9/50*f4v+2/55*f5v)
        v_error_appr=np.abs(v-v_prime)
        if flag==1:#Adaptive step size case
            h_x_new=0.9*h*(h*max_error/x_error_appr)**(1/4)
            h_v_new=0.9*h*(h*max_error/v_error_appr)**(1/4)
            h_new=min(h_x_new,h_v_new)
            if h_new<=h:
                h=h_new
            if h_new>h:
                x_list.append(x)
                x_prime_list.append(x_prime)
                x_error_appr_list.append(x_error_appr)
                v_list.append(v)
                v_prime_list.append(v_prime)
                v_error_appr_list.append(v_error_appr)
                h_list.append(h_new)
                t_list.append(t0)
                t0+=h
                x0=x
                v0=v
                h=h_new
                break
        if flag==0:#Fixed step size case
            x_list.append(x)
            x_prime_list.append(x_prime)
            x_error_appr_list.append(x_error_appr)
            v_list.append(v)
            v_prime_list.append(v_prime)
            v_error_appr_list.append(v_error_appr)
            h_list.append(h)
            t_list.append(t0)
            t0+=h
            x0=x
            v0=v
            break
#Graphs
#Comparison between solutions
#x vs. t
plt.figure(figsize=(9,5))
plt.xlim(t0_init-1,tf+1)
plt.plot(t_list,x_list,c="red",label="RKF4")
plt.plot(t_list,x_prime_list,c="green",label="RKF5")
plt.xlabel("$t$")
plt.ylabel("$x(t)$")
plt.title("Comparison between the numerical RKF45 solutions of $x$")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#v vs. t
plt.figure(figsize=(9,5))
plt.xlim(t0_init-1,tf+1)
plt.plot(t_list,v_list,c="red",label="RKF4")
plt.plot(t_list,v_prime_list,c="green",label="RKF5")
plt.xlabel("$t$")
plt.ylabel("$v(t)$")
plt.title("Comparison between the numerical RKF45 solutions of $v$")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#x vs. v (Phase Diagram)
plt.figure(figsize=(9,5))
plt.plot(x_list,v_list,c="red",label="RKF4")
plt.plot(x_prime_list,v_prime_list,c="green",label="RKF5")
plt.xlabel("$x$")
plt.ylabel("$v$")
plt.title("Comparison between the numerical RKF45 solutions (Phase Diagram)")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#Errors of the RKF45 method
plt.figure(figsize=(9,5))
plt.xlim(t0_init-1,tf+1)
plt.plot(t_list,x_error_appr_list,c="red",label="$x$")
plt.plot(t_list,v_error_appr_list,c="green",label="$v$")
plt.xlabel("$t$")
plt.ylabel("Error")
plt.title("Errors of the RKF45 method")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#Evolution of the step size 'h'
plt.figure(figsize=(9,5))
plt.xlim(t0_init-1,tf+1)
plt.plot(t_list,h_list,c="red")
plt.xlabel("$t$")
plt.ylabel("$h$")
plt.title("Evolution of the step size")
plt.grid()
plt.tight_layout()
plt.show()

print("Program execution time:",time.process_time()-start_time_program,"seconds.")

