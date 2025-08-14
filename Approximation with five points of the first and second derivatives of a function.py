# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 20:42:33 2021
Approximation with five points of the first and second derivatives of a function at a point x0=0 (numerical differentiation with finite differences).
More information on this website: https://en.wikipedia.org/wiki/Numerical_differentiation
@author: Lucas Romero Fernández
"""
from sympy import *
#main_program
#Definition of the constants and function
x0=pi/4
h=1#Value of the real and positive step size. Decrease it to improve accuracy.
x=Symbol("x")
def f(x):#Note: If a function like ln(x),e^(x),sqrt(x),etc. is put here, be sure to use the version of SymPy to avoid errors in the expression of the equation.
    return sin(x)
#Main process of the ''exact'' numerical obtainment of the first and second derivatives of a function at a point x0
def diff_func_ex():
    df=diff(f(x),x)
    df_0=N(df.subs(x,x0))
    d2f=diff(df)
    d2f_0=N(d2f.subs(x,x0))
    return df_0,d2f_0
#Main process of the ''approximated'' numerical obtainment of the first and second derivatives of a function at a point x0
def diff_func_approx():
    df_0=(1/(12*h))*(f(x0-2*h)-8*f(x0-h)+8*f(x0+h)-f(x0+2*h))
    d2f_0=(1/(12*h**2))*(-f(x0-2*h)+16*f(x0-h)-30*f(x0)+16*f(x0+h)-f(x0+2*h))
    return df_0,d2f_0
#Results
print("The ''exact'' values of the first and second derivatives of the function",f(x),"calculated numerically at the point x0 =",x0,"with h =",h,"is, respectively,",float(diff_func_ex()[0]),"and",float(diff_func_ex()[1]))
print("")
print("The ''approximated'' values of the first and second derivatives of the function",f(x),"calculated numerically at the point x0 =",x0,"with h =",h,"is, respectively,",float(diff_func_approx()[0]),"and",float(diff_func_approx()[1]))