# -*- coding: utf-8 -*-
"""
Created on Sat Oct 2 20:10:35 2021
Approximation with three points of the first and second derivatives of a function f(x) at a point x0=0 (numerical differentiation with finite differences).
More information on this website: https://en.wikipedia.org/wiki/Numerical_differentiation
@author: Lucas Romero Fernández
"""
from sympy import *
#main_program
#Definition of the constants and function
x0=0
h=0.1#Value of the real and positive step size. Decrease it to improve accuracy.
x=Symbol("x")
def f(x):#Note: If a function like ln(x),e^(x),sqrt(x),etc. is put here, be sure to use the version of SymPy to avoid errors in the expression of the equation.
    return sin(x)
#Main process of the ''exact'' numerical obtainment of the first and second derivatives of a function at a point x0
def diff_func_ex():
    d1=diff(f(x),x)
    d1_0=N(d1.subs(x,x0))
    d2=diff(d1)
    d2_0=N(d2.subs(x,x0))
    return d1_0,d2_0
#Main process of the ''approximated'' numerical obtainment of the first and second derivatives of a function at a point x0
def diff_func_approx():
    d1_0=(f(x0+h)-f(x0-h))/(2*h)
    d2_0=(f(x0+h)-2*f(0)+f(x0-h))/(h**2)
    return d1_0,d2_0
#Results
print("The ''exact'' values of the first and second derivatives of the function",f(x),"calculated numerically at the point x0 =",x0,"is, respectively,",float(diff_func_ex()[0]),"and",float(diff_func_ex()[1]))
print("")
print("The ''approximated'' values of the first and second derivatives of the function",f(x),"calculated numerically at the point x0 =",x0,"is, respectively,",float(diff_func_approx()[0]),"and",float(diff_func_approx()[1]))