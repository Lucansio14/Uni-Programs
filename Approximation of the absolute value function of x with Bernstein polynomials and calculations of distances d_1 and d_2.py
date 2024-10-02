# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 20:57:27 2021
Approximation of the absolute value function of x with Bernstein polynomials and calculations of distances d_1 and d_2.
More information about the Bernstein polynomials: https://en.wikipedia.org/wiki/Bernstein_polynomial
@author: Lucas Romero Fernández
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
#main_program
def f(x):#Function to approximate (analytical solution)
    return np.abs(x-1/2)#To represent graphically, the function has been moved to the interval [0,1]
def bernstein_polynomials(f,n,p):#Function used to obtain the Bernstein polynomials
    return np.sum([f(k/n)*st.binom.pmf(k,n,p) for k in np.arange(0,n+1)])
N=1000#Number of points in the x-axis
#Vectors of evenly spaced x-values in the interval [0,1] and corresponding y-values
x_vec=np.linspace(0,1,N)
y_vec=f(x_vec)
#For the approximation with Bernstein polynomials
bernstein3=lambda x:bernstein_polynomials(f,3,x)
bernstein3_vec=np.vectorize(bernstein3)
y3_vec=bernstein3_vec(x_vec)
bernstein7=lambda x:bernstein_polynomials(f,7,x)
bernstein7_vec=np.vectorize(bernstein7)
y7_vec=bernstein7_vec(x_vec)
bernstein15=lambda x:bernstein_polynomials(f,15,x)
bernstein15_vec=np.vectorize(bernstein15)
y15_vec=bernstein15_vec(x_vec)
bernstein30=lambda x:bernstein_polynomials(f,30,x)
bernstein30_vec=np.vectorize(bernstein30)
y30_vec=bernstein30_vec(x_vec)
bernstein60=lambda x:bernstein_polynomials(f,60,x)
bernstein60_vec=np.vectorize(bernstein60)
y60_vec=bernstein60_vec(x_vec)
#To calculate the distances
H_Order_poly=61#Highest order chosen for the Bernstein polynomials
n_list=np.arange(1,H_Order_poly)#Until the chosen Bernstein polynomial of highest order.
dist_1_list=[]
dist_2_list=[]
dist_inf_list=[]
for n in n_list:
    a=np.array(bernstein_polynomials(f,n,x_vec))
    b=np.array(f(x_vec))
    dist_1=np.sum(np.abs(a-b))
    dist_1_list.append(dist_1)
    dist_2=np.sqrt(np.sum(np.square(a-b)))
    dist_2_list.append(dist_2)
#Plots
#Approximation of the absolute value function of x with Bernstein polynomials
plt.figure(figsize=(8,8))
plt.plot(x_vec,y_vec,"blue",label="f(x)")
plt.plot(x_vec,y3_vec,"green",label="$B_{3}$")
plt.plot(x_vec,y7_vec,"red",label="$B_{7}$")
plt.plot(x_vec,y15_vec,"magenta",label="$B_{15}$")
plt.plot(x_vec,y30_vec,"brown",label="$B_{30}$")
plt.plot(x_vec,y60_vec,"yellow",label="$B_{60}$")
plt.axvline(0,color="black",linestyle="dashed",alpha=0.7)
plt.axvline(0.5,color="black",linestyle="dashed",alpha=0.7)
plt.axvline(1,color="black",linestyle="dashed",alpha=0.7)
plt.axhline(0,color="black",linestyle="dashed",alpha=0.7)
plt.xlabel("x")
plt.ylabel("y(x)")
plt.title(r'Approximations of |$x$|')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#Calculations of distances d_1 and d_2 vs order of the Bernstein polynomial
#d_1 vs order of the Bernstein polynomial
plt.figure(figsize=(8,8))
plt.plot(n_list,dist_1_list,"blue")
plt.axvline(1,color="black",linestyle="dashed",alpha=0.7)
plt.axvline(H_Order_poly-1,color="black",linestyle="dashed",alpha=0.7)
plt.xlabel("Order of the Bernstein polynomial")
plt.ylabel("$d_{1}$")
plt.title("Evolution of $d_{1}$")
plt.grid()
plt.tight_layout()
plt.show()
#d_2 vs order of the Bernstein polynomial
plt.figure(figsize=(8,8))
plt.plot(n_list,dist_2_list,"red")
plt.axvline(1,color="black",linestyle="dashed",alpha=0.7)
plt.axvline(H_Order_poly-1,color="black",linestyle="dashed",alpha=0.7)
plt.xlabel("Order of the Bernstein polynomial")
plt.ylabel("$d_{2}$")
plt.title("Evolution of $d_{2}$")
plt.grid()
plt.tight_layout()
plt.show()
#Comparison of d_1 and d_2 vs order of the Bernstein polynomial
plt.figure(figsize=(8,8))
plt.plot(n_list,dist_1_list,"blue",label="$d_{1}$")
plt.plot(n_list,dist_2_list,"red",label="$d_{2}$")
plt.axvline(1,color="black",linestyle="dashed",alpha=0.7)
plt.axvline(H_Order_poly-1,color="black",linestyle="dashed",alpha=0.7)
plt.xlabel("Order of the Bernstein polynomial")
plt.ylabel("Distances")
plt.title("Evolution of distances $d_{1}$ and $d_{2}$")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()