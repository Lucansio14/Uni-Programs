# -*- coding: utf-8 -*-
"""
Created on Wed Oct 2 13:22:20 2019
Square Root Calculator with the Babylonian Method
@author: Lucas Romero Fernández
"""
from timeit import timeit
#More information about the method in this website: https://blogs.sas.com/content/iml/2016/05/16/babylonian-square-roots.html
#main_program
def square_root(x):
    if x<0:
        x=-x
        v=1j
    else:
        v=1
    r=x
    t=0
    while t!=r:
        t=r
        r=1/2*((x/r)+r)
    return r*v
n=float(input("Number? "))
print("")
print("The square root of ",n," is ",square_root(n),".",sep="")
print("")
print("Program execution time:",timeit(number=1),"seconds.")
print("")
print("Goodbye.")