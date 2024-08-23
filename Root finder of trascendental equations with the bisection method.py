# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 13:31:33 2024
Root finder of trascendental equations with the bisection method
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
import sys
import sympy as sy
#More information about the method in this website: https://en.wikipedia.org/wiki/Bisection_method
#main_program
#Definition of fixed variables
start_time_program=time.process_time()#To calculate the program execution time
a=-1#Endpoint a of the interval
b=3#Endpoint b of the interval
x_array=np.linspace(a,b,1000,dtype=float)#To represent f(x) in a graph
f_list=[]#To represent f(x) in a graph
a_points_list=[]#To represent how the interval changes in the graph
b_points_list=[]#To represent how the interval changes in the graph
fa_points_list=[]#To represent how the interval changes in the graph
fb_points_list=[]#To represent how the interval changes in the graph
def f(x):#Equation to solve
    return 2*sy.exp(-x)-x+2**(-x)#Note: If a function like ln(x),e^(x),sqrt(x),etc. is put here, be sure to use the version of SymPy to avoid errors in the expression of the equation.
x=sy.Symbol("x")
#Main process of the bisection method
def bisecroot(f,a,b,verbose=False):
    """
    Returns a root of the function f in the interval [a,b] (assuming b>a)
    For this method to work, f(a) and f(b) must have different signs.
    It shows information about the process if verbose is active (for debugging reasons). 
    """
    if b<a:#In case a and b are reversed.
        a,b=b,a
    if verbose==True:
        print("Trying to find a root of",f(x),"between",a,"and",b)
    a_points_list.append(a)
    fa_points_list.append(f(a))
    b_points_list.append(b)
    fb_points_list.append(f(b))
    if f(a)==0:#For situations where a is already one of the roots...
        return a
    if f(b)==0:#For situations where b is already one of the roots...
        return b
    if f(a)*f(b)>0:#If f(a) and f(b) do not have different signs...
        print("f(a)=",f(a),"f(b)=",f(b))
        sys.exit("f(a) and f(b) must have different signs for the method to be succesful.")
    old_posib_sol=a
    while True:
        posib_sol=(a+b)/2.0
        if verbose==True: 
            print("Trying with middle point:",posib_sol)
        if f(posib_sol)==0 or posib_sol==old_posib_sol: 
            if verbose==True:
                if f(posib_sol)==0: 
                    print("Valued zero",f(posib_sol))
                else:
                    print("Not changing",posib_sol,old_posib_sol)
            break
        if f(a)*f(posib_sol)<0: 
            b=posib_sol#Solution is in the interval [a,m]
            b_points_list.append(b)
            fb_points_list.append(f(b))
        else:
            a=posib_sol#Solution is in the interval [m,b]
            a_points_list.append(a)
            fa_points_list.append(f(a))
        old_posib_sol=posib_sol
    return posib_sol
sol=bisecroot(f,a,b)
x_sol=sol
fx_sol=f(sol)
#Gathering values for the graph
for x1 in x_array:#To represent f(x) in a graph...
    r=f(x1)
    f_list.append(r)
f_array=np.array(f_list,dtype=float)
a_points_array=np.array(a_points_list,dtype=float)
fa_points_array=np.array(fa_points_list,dtype=float)
b_points_array=np.array(b_points_list,dtype=float)
fb_points_array=np.array(fb_points_list,dtype=float)
a_num_array=np.arange(0,len(a_points_array),1)#To identify each point
b_num_array=np.arange(0,len(b_points_array),1)#To identify each point
#Graph
plt.figure(figsize=(9,5))
plt.plot(x_array,f_array,color="black",alpha=0.7,label="f(x)")
plt.axhline(0,color="black",linestyle="dashed",alpha=0.6)
plt.scatter(a_points_array,fa_points_array,color="blue",marker="|",s=100,label="Endpoint a")
for i in a_num_array:
    plt.vlines(a_points_array[i],0,fa_points_array[i],color="black",linestyle="dashed",alpha=0.6)
plt.scatter(b_points_array,fb_points_array,color="red",marker="|",s=100,label="Endpoint b")
for i in b_num_array:
    plt.vlines(b_points_array[i],0,fb_points_array[i],color="black",linestyle="dashed",alpha=0.6)
plt.scatter(x_sol,fx_sol,color="green",marker="d",s=100,label="Solution")
plt.xlabel("$x$")
plt.ylabel("$f(x)$")
plt.grid()
plt.title("Root finder of trascendental equations with the bisection method")
plt.legend()
plt.tight_layout()
plt.show()
#Results
print("Solution of",f(x),"between a =",a,"and b =",b,": x =",sol)
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")