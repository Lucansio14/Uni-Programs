# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 20:16:00 2021
Integration by trapezoidal rule (Integrate a function by adding the areas of the trapezoids below the function)
Simplified version (general and shorter algorithm without plots with trapezoids)
More information about the method in this website: https://en.wikipedia.org/wiki/Trapezoidal_rule
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import scipy as sp
#Definition of variables, constants, lists and functions
a=0#Lower bound of the integration range
b=np.log(3)#Upper bound of the integration range
n=int(10000)#(integer) Number of equal subdivisions of the integration range. Increase it to gain precision if necessary. Must be an integer and bigger than 1.
n_array=np.arange(n+1)
x_list=[]
Num_Sol=0#Initial value of the numerical solution with the general trapezoidal rule
def f(x):#Integrand of the definite integral
    return np.exp(x)
#main_program
start_time_program=time.process_time()#To calculate the program execution time
Sc_Sol=sp.integrate.quad(lambda x:f(x),a,b)[0]#Main process of the obtainment of the numerical solution with SciPy
#Main process of the obtainment of the numerical solution with the general trapezoidal rule
h=(b-a)/n#Step size
for i in range(2,len(n_array)):
    x_list.append(a+n_array[i]*h)
for i in range(2,len(x_list)):
    Num_Sol+=(1/2)*h*(f(x_list[i-1])+f(x_list[i]))
Error_Sol=abs(Num_Sol-Sc_Sol)
#Results of the integrations
print("Value of the numerical integration by trapezoidal rule =",Num_Sol)#Result of the numerical integration by trapezoidal rule
print("")
print("Value of the numerical integration with SciPy =",Sc_Sol)#Result of the numerical integration with SciPy
print("")
print("Value of the error =",Error_Sol)#Result of the error between methods
print("")
print("Value of the absolute relative error =",abs(((Error_Sol)/(Sc_Sol)))*100,"%")#Result of the absolute relative error between methods
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")