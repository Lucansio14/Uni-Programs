# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 18:15:01 2024
Calculus of the Gaussian curvature (k) in a geometry (2D only) with coordinates (x1 and x2) specified by the corresponding metric tensor (g_ab).
More information about the method in this website: https://randomphysics.com/wp-content/uploads/2023/03/ferrari_gualtieri_general_relativity.pdf (Chapter 1)
@author: Lucas Romero Fernández
"""
import time
from sympy import *
import numpy as np
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#First example (spherical surface). Change to your liking if necessary.
#Definition of symbols and constants
x0=Symbol("x0")#Constant. In this case, is r
x1=Symbol("x1")#In this case, the first coordinate is theta
x2=Symbol("x2")#In this case, the second coordinate is phi
g_11=simplify(x0**2)
g_22=simplify(x0**2*(sin(x1))**2)
g_12=simplify(0)
g_21=g_12#Applying symmetry of the metric tensor
g_ab=np.array([[g_11,g_12],[g_21,g_22]])
g=g_11*g_22-g_12**2
print("Metric tensor g_ab:")
print(g_ab)
print("")
print("Calculating k...")
print("")
#Main process of the obtainment of k in 2D
k=simplify((1/(2*g))*(2*diff(g_12,x1,x2)-diff(g_11,x2,x2)-diff(g_22,x1,x1))-(g_22/(4*g**2))*(diff(g_11,x1)*(2*diff(g_12,x2)-diff(g_22,x1))-(diff(g_11,x2))**2)+(g_12/(4*g**2))*(diff(g_11,x1)*diff(g_22,x2)-2*diff(g_11,x2)*diff(g_22,x1)+(2*diff(g_12,x1)-diff(g_11,x2))*(2*diff(g_12,x2)-diff(g_22,x1)))-(g_11/(4*g**2))*(diff(g_22,x2)*(2*diff(g_12,x1)-diff(g_11,x2))-(diff(g_22,x1))**2))
#Result
print("k =",k)
print("")
#Second example (Gauss-Bolyai-Lobachewski geometry). Change to your liking if necessary.
#Definition of symbols and constants
x0=Symbol("x0")#Constant. In this case, is r
x1=Symbol("x1")
x2=Symbol("x2")
g_11=simplify((x0**2*(1-(x2)**2))/((1-x1**2-x2**2)**2))
g_22=simplify((x0**2*(1-(x1)**2))/((1-x1**2-x2**2)**2))
g_12=simplify((x0**2*x1*x2)/((1-x1**2-x2**2)**2))
g_21=g_12#Applying symmetry of the metric tensor
g_ab=np.array([[g_11,g_12],[g_21,g_22]])
g=g_11*g_22-g_12**2
print("Metric tensor g_ab:")
print(g_ab)
print("")
print("Calculating k...")
print("")
#Main process of the obtainment of k in 2D
k=simplify((1/(2*g))*(2*diff(g_12,x1,x2)-diff(g_11,x2,x2)-diff(g_22,x1,x1))-(g_22/(4*g**2))*(diff(g_11,x1)*(2*diff(g_12,x2)-diff(g_22,x1))-(diff(g_11,x2))**2)+(g_12/(4*g**2))*(diff(g_11,x1)*diff(g_22,x2)-2*diff(g_11,x2)*diff(g_22,x1)+(2*diff(g_12,x1)-diff(g_11,x2))*(2*diff(g_12,x2)-diff(g_22,x1)))-(g_11/(4*g**2))*(diff(g_22,x2)*(2*diff(g_12,x1)-diff(g_11,x2))-(diff(g_22,x1))**2))
#Result
print("k =",k)
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")