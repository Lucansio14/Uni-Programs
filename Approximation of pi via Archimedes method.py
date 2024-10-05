# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 11:37:23 2021
Calculate an aproximation of the number pi by using a variation of the Archimedes method.
More information about the method in this website: https://www.craig-wood.com/nick/articles/pi-archimedes/
@author: Lucas Romero Fernández
"""
import time
import matplotlib.pyplot as plt
import numpy as np
#main_program
start_time_program=time.process_time()#To calculate the program execution time
def f(n):#Formulae/expression to approximate pi with the Archimedes Method
    return (np.sin(2*np.pi/n)*n)/(2)
#Definition of lists and constants
n=96#Number of iterations. Increase it to gain precision. Must be a non zero positive integer bigger than 2
pi_list=[]
error_pi_list=[]
#Main process of the obtainment of pi with the Archimedes method
for k in range(3,n+1):#Note: This unusual start of the list is to eliminate the values of n corresponding to non poligons
    pi_list.append(f(k))
#Main process of the obtainment of pi with NumPy
numpy_pi=np.pi
for k in range(0,len(pi_list)):#Error between methods
    error_pi_list.append(numpy_pi-pi_list[k])
n_list=np.arange(3,n+1)#For the graphs
#Graphs
#Approximations of the number pi vs number of integrations (n)
plt.figure(figsize=(9,5))
plt.xlim(3,n)
plt.plot(n_list,pi_list,c='red')
plt.axhline(numpy_pi,color="black",linestyle="dashed",alpha=0.7,label="NumPy pi")
plt.title("Archimedes method to calculate pi")
plt.ylabel("Approximations of the number pi")
plt.xlabel("Number of integrations (n)")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#Error of the method vs number of integrations (n)
plt.figure(figsize=(9,5))
plt.xlim(3,n)
plt.plot(n_list,error_pi_list,c='red')
plt.title("Error of the Archimedes method")
plt.ylabel("Error of the method")
plt.xlabel("Number of integrations (n)")
plt.grid()
plt.tight_layout()
plt.show()
#Results
print("Value of pi with the Archimedes method (n = {})".format(n),"=",pi_list[-1])#Result of the obtainment of pi with the Archimedes method
print("")
print("Value of pi with NumPy =",numpy_pi)#Result of the obtainment of pi with NumPy
print("")
print("Value of the final error =",error_pi_list[-1])#Result of the final error between methods
print("")
print("Value of the final absolute relative error =",abs(((error_pi_list[-1])/(numpy_pi)))*100,"%")#Result of the final absolute relative error between methods
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")