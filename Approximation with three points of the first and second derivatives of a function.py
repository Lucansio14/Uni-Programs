# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 20:10:35 2021
Approximation with three points of the first and second derivatives of a function at a point x0=0 (numerical differentiation with finite differences).
It is assumed that the function fits perfectly to the second order polynomial function ax^2+bx+c that goes through the three points.
More information on this website: https://en.wikipedia.org/wiki/Numerical_differentiation
@author: Lucas Romero Fernández
"""
from sympy import *
#main_program
#Definition of lists and constants
x0=0
x_values=[-1,0,1]#Values of x where the function is studied
h=1#Value of the real and positive step size, must match the differences between the x values
#Main process of the analytical obtainment of the first and second derivatives of a function at a point x0
def pol_2order_an():
    x=Symbol('x')
    a,b,c=symbols("a,b,c")
    d1=diff(a*x**(2)+b*x+c,x)
    d1_0=N(d1.subs(x,x0))
    d2=diff(d1,x)
    d2_0=N(d2.subs(x,x0))
    print("Analitically, the first derivative of the second order polynomial function at the point x0 =",x0,"is equal to",d1_0,"and the second derivative of the same function at that same point is equal to",d2_0)
    print("")
pol_2order_an()
#Main process of the numerical obtainment of the first and second derivatives of a function at a point x0
a=float(input("Value of a: "))
print("")
b=float(input("Value of b: "))
print("")
c=float(input("Value of c: "))#As can be seen, this value does not affect the derivatives. Left here so it can be verified personally
print("")
def pol_2order_ex():
    x=Symbol('x')
    d1=diff(a*x**(2)+b*x+c,x)
    d1_0=N(d1.subs(x,x0))
    d2=diff(d1)
    d2_0=N(d2.subs(x,x0))
    return d1_0,d2_0
def pol_2order_cal(x):#For the calculus of the approximations
    return a*x**(2)+b*x+c
#Results
print("The values of the first and second derivatives of the second order polynomial function calculated analitically at the point x0 =",x0,"is, respectively,",float(pol_2order_ex()[0]),"and",float(pol_2order_ex()[1]))
print("")
print("The approximated values of the first and second derivatives of the second order polynomial function calculated numerically at the point x0 =",x0,"is, respectively,",(pol_2order_cal(1)-pol_2order_cal(-1))/(2*h),"and",(pol_2order_cal(1)-2*pol_2order_cal(0)+pol_2order_cal(-1))/(h**2))
