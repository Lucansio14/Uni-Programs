# -*- coding: utf-8 -*-
"""
Created on Wed Oct 2 13:22:20 2019
Square Root Calculator with the Babylonian Method
@author: Lucas Romero Fernández
"""
from timeit import timeit
import time
import numpy as np
#More information about the method in this website: https://blogs.sas.com/content/iml/2016/05/16/babylonian-square-roots.html
#main_program
start_time_program=time.process_time()#To calculate the program execution time
def square_root_BM(x):#Babylonian method
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
def square_root_SQRT(x):#sqrt method
    return np.emath.sqrt(n)#np.emath is used instead of something more normal like cmath because it is better not to show the imaginary part when it is not needed...
#For the timeit module
square_root_BM_string='''
def square_root_BM(x):
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
    '''
square_root_SQRT_string='''
def square_root_SQRT(x):
    return np.emath.sqrt(n)
    '''
#Results
n=float(input("Number? "))
print("")
print("With the babilonian method, the square root of ",n," is ",square_root_BM(n),".",sep="")
print("")
print("Execution time of the babilonian method:",timeit(square_root_BM_string,number=1),"seconds.")
print("")
print("With the built-in numpy function emath.sqrt, the square root of ",n," is ",square_root_SQRT(n),".",sep="")
print("")
print("Execution time of the numpy function emath.sqrt:",timeit(square_root_SQRT_string,"import numpy as np",number=1),"seconds.")
print("")
print("Value of the error between methods =",abs(square_root_SQRT(n)-square_root_BM(n)))
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")
print("")
print("Goodbye.")
