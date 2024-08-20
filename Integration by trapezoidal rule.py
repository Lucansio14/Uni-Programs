# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 20:37:12 2020
Integration by trapezoidal rule (Integrate a function by adding the areas of the trapezoids below the function) with uniform grid.
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
import sympy as sy
#More information about the method in this website: https://en.wikipedia.org/wiki/Trapezoidal_rule
#main_program
#Definition of fixed variables
start_time_program=time.process_time()#To calculate the program execution time
def f(x):#Function to integrate
    return sy.sqrt(1+(4*x-1)**(2))*sy.exp(x)#Note: If a function like ln(x),e^(x),sqrt(x),etc. is put here, be sure to use the version of SymPy to avoid errors.
x0=0#Integration domain from x0 to xn
xn=5#Integration domain from x0 to xn
num=8#Number of points. Increase it to gain precision if necessary. Must be an integer and bigger than 1.
num_array=np.arange(0,num,1)#To identify each point
x=sy.Symbol("x")
print("Calculating the integral of",f(x),"from x =",x0,"to x =",xn,"using SymPy and the trapezoidal rule with",num,"points...")
#Main process of the numerical integration with SymPy
I_exact=sy.integrate(f(x),(x,x0,xn))
#Main process of the numerical integration by trapezoids of f(x) between x0 and xn
dx=(xn-x0)/num#Numerical value of the "step" (Length of the base of the trapezoid)
add=(f(x0)+f(xn))/2
x=x0
for i in range(1,num):
    x+=dx
    add+=f(x)
I=dx*add
E_t=abs(I_exact-I)#Error between methods
#Gathering values for the graph
x_array=np.linspace(x0,xn,num,dtype=float)#To represent f(x) obtained by the trapezoidal rule in a graph...
x_exact_array=np.linspace(x0,xn,1000,dtype=float)#To represent f(x) in a graph...
f_list=[]
f_exact_list=[]
for x in x_array:#To represent f(x) obtained by the trapezoidal rule in a graph...
    a=f(x)
    f_list.append(a)
for x in x_exact_array:#To represent f(x) in a graph...
    a=f(x)
    f_exact_list.append(a)
f_array=np.array(f_list,dtype=float)
f_exact_array=np.array(f_exact_list,dtype=float)
fx0_array=np.full(num,f_array[0],dtype=float)#Horizontal line of value f_array[0] for filling under the graph
#Graph
fig=plt.figure(figsize=(9,5))
plt.plot(x_array,f_array,marker="d",color="orange",label="Trapezoidal rule")
plt.plot(x_exact_array,f_exact_array,color="black",label="f(x)")
plt.hlines(f_array[0],x_array[0],x_array[-1],color="black",linestyle="dashed",alpha=0.6)
for i in num_array:
    plt.vlines(x_array[i],f_array[0],f_array[i],color="black",linestyle="dashed",alpha=0.6)
plt.fill_between(x_array,f_array,fx0_array,where=((x_array[0]-dx)<x_array)&(x_array<(x_array[-1]+dx)),color="orange",alpha=0.5)
plt.xlabel("$x$")
plt.ylabel("$f(x)$")
plt.grid()
plt.title("Integration by trapezoidal rule")
plt.legend()
plt.tight_layout()
plt.show()
#Results of the integrations
#Note:"().evalf()" is put when using SymPy in f(x) to always print a numerical (sometimes, very approximated) value, no matter the complexity of the result...
print("Value of the numerical integration with SymPy =",(I_exact).evalf())#Result of the numerical integration with SymPy
print("")
print("Value of the numerical integration by trapezoidal rule =",(I).evalf())#Result of the numerical integration by trapezoidal rule
print("")
print("Value of the error =",(E_t).evalf())#Result of the error between methods
print("")
print("Value of the absolute relative error =",(abs(((E_t)/(I_exact)))*100).evalf(),"%")#Result of the absolute relative error between methods
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")
