# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 21:19:31 2021
Calculate an approximation of the number pi by using the Leibniz formulae.
More information about the method in this website: https://en.wikipedia.org/wiki/Leibniz_formula_for_%CF%80
@author: Lucas Romero Fernández
"""
import time
import matplotlib.pyplot as plt
import numpy as np
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants and lists
k=0#Summation index
n=250#Summation upper bound/number of iterations, must be an integer
S_list=[]
Leib_pi_list=[]
error_Leib_pi_list=[]
#Main process of the obtainment of pi with the Leibniz formulae
while k<=n:
    S=((-1)**(k)/(2*k+1))
    S_list.append(S)
    Leib_pi=4*sum(S_list)
    Leib_pi_list.append(Leib_pi)
    k+=1
#Main process of the obtainment of pi with NumPy
numpy_pi=np.pi
for i in range(0,len(Leib_pi_list)):#Error between methods
    error_Leib_pi_list.append(abs(numpy_pi-Leib_pi_list[i]))
iterations_list=np.arange(0,n+1,1)#For the graphs
#Graphs
#Approximations of the number pi vs number of iterations
plt.figure(figsize=(9,5))
plt.plot(iterations_list,Leib_pi_list,c='red')
plt.axhline(numpy_pi,color="black",linestyle="dashed",alpha=0.7,label="NumPy pi")
plt.title("Calculus of the number pi with the Leibniz formulae")
plt.ylabel("Approximations of the number pi")
plt.xlabel("Number of iterations")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#Error of the method vs number of iterations
plt.figure(figsize=(9,5))
plt.plot(iterations_list,error_Leib_pi_list,c='red')
plt.title("Error of the Leibniz formulae")
plt.ylabel("Error of the expression")
plt.xlabel("Number of iterations")
plt.grid()
plt.tight_layout()
plt.show()
#Results
print("Value of pi with the Leibniz formulae (n = {})".format(n),"=",Leib_pi_list[-1])#Result of the obtainment of pi with the Leibniz formulae
print("")
print("Value of pi with NumPy =",numpy_pi)#Result of the obtainment of pi with NumPy
print("")
print("Value of the final error =",error_Leib_pi_list[-1])#Result of the final error between methods
print("")
print("Value of the final absolute relative error =",abs(((error_Leib_pi_list[-1])/(numpy_pi)))*100,"%")#Result of the final absolute relative error between methods
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")
