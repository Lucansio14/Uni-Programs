# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 09:53:09 2020
Calculator for roots of polynomials up to seventh degree
@author: Lucas Romero Fernández
"""
import sys
from sympy import*
#main_program
print("Calculator for roots of polynomials (real or imaginary) up to seventh degree of type ax^7+bx^6+cx^5+dx^4+ex^3+fx^2+gx+h=0")
a=float(input("Value of a: "))
b=float(input("Value of b: "))
c=float(input("Value of c: "))
d=float(input("Value of d: "))
e=float(input("Value of e: "))
f=float(input("Value of f: "))
g=float(input("Value of g: "))
h=float(input("Value of h: "))
if (a==float(0)) and (b==float(0)) and (c==float(0)) and (d==float(0)) and (e==float(0)) and (f==float(0)) and (g==float(0)):
    sys.exit("There is no equation to solve...")
print("-----------------------------------------------------------")
print("Solving equation of the form:",a,"x^7 +",b,"x^6 +",c,"x^5 +",d,"x^4 +",e,"x^3 +",f,"x^2 +",g,"x +",h,"= 0:")
print("...")
print("...")
print("...")
x=Symbol('x')
def RootFinder(x):
    return solve(a*x**7+b*x**6+c*x**5+d*x**4+e*x**3+f*x**2+g*x+h,x)
print("Solutions:",RootFinder(x))