# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 18:15:01 2024
Calculus of the metric tensor (g_ab_new) in new coordinate system (2D only), knowing, the metric tensor in 2D (g_ab) for the present/old coordinate system (optionally)
and the expressions of the old coordinates (xi1 and xi2) in respect of the new coordinates (x1 and x2), xi1(x1,x2) and xi2(x1,x2).
More information about the method in this website: https://randomphysics.com/wp-content/uploads/2023/03/ferrari_gualtieri_general_relativity.pdf (Chapter 1)
@author: Lucas Romero Fernández
"""
import time
from sympy import *
import numpy as np
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of symbols and constants
x1=Symbol("x1")#In this case, the new first coordinate is r (polar coordinates)
x2=Symbol("x2")#In this case, the new second coordinate is theta (polar coordinates)
xi1=x1*cos(x2)#In this case, cartesian coordinate x
xi2=x1*sin(x2)#In this case, cartesian coordinate y
g_11_old=simplify(1)#In this case, the old coordinates are cartesian coordinates
g_22_old=simplify(1)#In this case, the old coordinates are cartesian coordinates
g_12_old=simplify(0)#In this case, the old coordinates are cartesian coordinates
g_21_old=g_12_old#Applying symmetry of the metric tensor
g_old=np.array([[g_11_old,g_12_old],[g_21_old,g_22_old]])
print("Old metric tensor g_ab:")
print(g_old)
print("")
print("Calculating new metric tensor...")
print("")
#Main process of the obtainment of the metric tensor in 2D (main method) and result
g_11_new=simplify((diff(xi1,x1))**2+(diff(xi2,x1))**2)
g_22_new=simplify((diff(xi1,x2))**2+(diff(xi2,x2))**2)
g_12_new=simplify(diff(xi1,x1)*diff(xi1,x2)+diff(xi2,x1)*diff(xi2,x2))
g_21_new=g_12_new#Applying symmetry of the metric tensor
g_new=np.array([[g_11_new,g_12_new],[g_21_new,g_22_new]])
print("New metric tensor g_ab_new (main method):")
print(g_new)
print("")
#Main process of the obtainment of the metric tensor in 2D (alternative method) and result
g_11_new_alt=simplify(g_11_old*diff(xi1,x1)*diff(xi1,x1)+g_12_old*diff(xi1,x1)*diff(xi2,x1)+g_21_old*diff(xi2,x1)*diff(xi1,x1)+g_22_old*diff(xi2,x1)*diff(xi2,x1))
g_22_new_alt=simplify(g_11_old*diff(xi1,x2)*diff(xi1,x2)+g_12_old*diff(xi1,x2)*diff(xi2,x2)+g_21_old*diff(xi2,x2)*diff(xi1,x2)+g_22_old*diff(xi2,x2)*diff(xi2,x2))
g_12_new_alt=simplify(g_11_old*diff(xi1,x1)*diff(xi1,x2)+g_12_old*diff(xi1,x1)*diff(xi2,x2)+g_21_old*diff(xi2,x1)*diff(xi1,x2)+g_22_old*diff(xi2,x1)*diff(xi2,x2))
g_21_new_alt=g_12_new_alt#Applying symmetry of the metric tensor
g_new_alt=np.array([[g_11_new_alt,g_12_new_alt],[g_21_new_alt,g_22_new_alt]])
print("New metric tensor g_ab_new (alternative method):")
print(g_new_alt)
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")