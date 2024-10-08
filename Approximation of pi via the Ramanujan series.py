# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 12:08:35 2021
Calculate an aproximation of the number pi by using the Ramanujan series.
More information in this website: https://planetmath.org/ramanujansformulaforpi
@author: Lucas Romero Fernández
"""
import time
import math
import matplotlib.pyplot as plt
import numpy as np
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants and lists
k=0#Summation index
n=1#Summation upper bound, must be an integer
S_list=[]
Ram_pi_list=[]
error_Ram_pi_list=[]
#Main process of the obtainment of pi with the Ramanujan series
while k<=n:
    S=(math.factorial(4*k))*(1103+26390*k)/(((math.factorial(k))**(4))*(396)**(4*k))
    S_list.append(S)
    Ram_pi=(((2*np.sqrt(2))/9801)*sum(S_list))**(-1)
    print(Ram_pi)
    Ram_pi_list.append(Ram_pi)
    k+=1
#Main process of the obtainment of pi with NumPy
numpy_pi=np.pi
for i in range(0,len(Ram_pi_list)):#Error between methods
    error_Ram_pi_list.append(abs(numpy_pi-Ram_pi_list[i]))
iterations_list=np.arange(0,n+1,1)#For the graphs
#Graphs
#Approximations of the number pi vs number of iterations
plt.figure(figsize=(9,5))
plt.plot(iterations_list,Ram_pi_list,c='red')
plt.axhline(numpy_pi,color="black",linestyle="dashed",alpha=0.7,label="NumPy pi")
ax=plt.gca()
ax.yaxis.get_major_formatter().set_scientific(False)
ax.yaxis.get_major_formatter().set_useOffset(False)
plt.title("Calculus of the number pi with the Ramanujan series")
plt.ylabel("Approximations of the number pi")
plt.xlabel("Number of iterations")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#Error of the method vs number of iterations
plt.figure(figsize=(9,5))
plt.plot(iterations_list,error_Ram_pi_list,c='red')
plt.title("Error of the Ramanujan series")
plt.ylabel("Error of the expression")
plt.xlabel("Number of iterations")
plt.grid()
plt.tight_layout()
plt.show()
#Results
print("Value of pi with the Ramanujan series (n = {})".format(n),"=",Ram_pi_list[-1])#Result of the obtainment of pi with the Ramanujan series
print("")
print("Value of pi with NumPy =",numpy_pi)#Result of the obtainment of pi with NumPy
print("")
print("Value of the final error =",error_Ram_pi_list[-1])#Result of the final error between methods
print("")
print("Value of the final absolute relative error =",abs(((error_Ram_pi_list[-1])/(numpy_pi)))*100,"%")#Result of the final absolute relative error between methods
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")
