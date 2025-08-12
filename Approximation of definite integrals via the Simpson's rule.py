# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 17:50:23 2021
Calculate an approximation of an arbitrary definite integral by using the Simpson's rule.
More information about the method in this website: https://en.wikipedia.org/wiki/Simpson%27s_rule
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
#Definition of variables, constants, lists and functions
a=0#Lower bound of the integration range
b=np.pi#Upper bound of the integration range
n_min=int(1)#Minimum (integer) number of equal subdivisions of the integration range
n_max=int(20)#Maximum (integer) number of equal subdivisions of the integration range
n_step=int(1)#Spacing between the n values (must be an integer)
n_array=np.arange(n_min,n_max+n_step,n_step)
Num_Sol_list=[]
Error_Sol_list=[]
def f(x):#Integrand of the definite integral
    return (3/2)*(np.sin(x))**3
#main_program
start_time_program=time.process_time()#To calculate the program execution time
Sc_Sol=sp.integrate.quad(lambda x:f(x),a,b)[0]#Main process of the obtainment of the numerical solution with SciPy
#Main process of the obtainment of the numerical solution with the Simpson's rule
for n in n_array:
    h=(b-a)/(n)#Step size
    i=1
    j=1
    Sum1=0
    Sum2=0
    while i<=(n/2):
        Sum1+=f(a+(2*i-1)*h)
        i+=1
    while j<=((n/2)-1):
        Sum2+=f(a+2*j*h)
        j+=1
    Term3=4*Sum1
    Term4=2*Sum2
    Num_Sol=(h/3)*(f(a)+f(b)+Term3+Term4)
    Num_Sol_list.append(Num_Sol)
    Error_Sol=abs(Num_Sol-Sc_Sol)
    Error_Sol_list.append(Error_Sol)
#Graphs
#Approximations of the solution via Simpson's rule vs n
figure,axes=plt.subplots(figsize=(9,5))
axes.plot(n_array,Num_Sol_list,c='red')
axes.axhline(Sc_Sol,color="black",linestyle="dashed",alpha=0.7,label="ScyPy solution")
axes.set_title("Calculus of the definite integral with integrand f(x) using Simpson's rule")
axes.set_ylabel("Value of the definite integral")
axes.set_xlabel("n")
axes.set_xlim(n_min,n_max)
axes.xaxis.get_major_locator().set_params(integer=True)
axes.grid(True)
figure.set_tight_layout(True)
plt.legend()
plt.show()
#Error of the method vs n
figure,axes=plt.subplots(figsize=(9,5))
axes.plot(n_array,Error_Sol_list,c='red')
axes.set_title("Error of the definite integral's solution using Simpson's rule")
axes.set_ylabel("Error of the solution")
axes.set_xlabel("n")
axes.set_xlim(n_min,n_max)
axes.xaxis.get_major_locator().set_params(integer=True)
axes.grid(True)
figure.set_tight_layout(True)
plt.show()
#Results
print("Value of the definite integral of integrand f(x) using Simpson's rule (n = n_max = {})".format(n_max),"=",Num_Sol_list[-1])#Result of the definite integral using the Simpson's rule
print("")
print("Value of the definite integral of integrand f(x) using SciPy =",Sc_Sol)#Result of the definite integral using SciPy
print("")
print("Value of the final error =",Error_Sol_list[-1])#Result of the final error between methods
print("")
print("Value of the final absolute relative error =",abs(((Error_Sol_list[-1])/(Sc_Sol)))*100,"%")#Result of the final absolute relative error between methods
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")