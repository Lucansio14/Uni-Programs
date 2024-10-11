# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 19:56:32 2021
Calculate an approximation of the piecewise function f(t) by using the corresponding expression in Fourier series S(t;n).
       1 if 0 < t < T/2,
f(t) = 0 if t=T/2,        where T = 2*pi, t = alpha*T and alpha = 0.01, 0.25, 0.49 (NOTE: The values of alpha used are this ones because the Fourier series is only correct for 0 < t < T/2).
      -1 if T/2 < t < T,
More information of Fourier series in this website: https://en.wikipedia.org/wiki/Fourier_series
@author: Lucas Romero Fernández
"""
import time
import matplotlib.pyplot as plt
import numpy as np
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants, functions and lists
T=2*np.pi
i=1#Summation index
n=100#Summation upper bound for the series, must be an integer
alpha_list=[0.01,0.25,0.49]
S_sum_list=[]
S_list001=[]
S_list025=[]
S_list049=[]
err_list=[]
n_list=np.arange(1,n+1)
def f(alpha):#Expressed in a simplified manner
    t=alpha*T
    if t<T/2:
        return 1
    if t==T/2:
        return 0
    if t>T/2:
        return -1
#Main process of the obtainment of the approximation of f(t) with the corresponding Fourier series
#alpha = 0.01
alpha=alpha_list[0]
while i<=n:
    S_sum=(1/(2*i-1))*np.sin(np.radians((2*(2*i-1)*np.pi*alpha)))
    S_sum_list.append(S_sum)
    S=(4/np.pi)*sum(S_sum_list)
    S_list001.append(S)
    i+=1
#alpha = 0.25
i=1
S_sum_list=[]
alpha=alpha_list[1]
while i<=n:
    S_sum=(1/(2*i-1))*np.sin(np.radians((2*(2*i-1)*np.pi*alpha)))
    S_sum_list.append(S_sum)
    S=(4/np.pi)*sum(S_sum_list)
    S_list025.append(S)
    i+=1
#alpha = 0.49
i=1
S_sum_list=[]
alpha=alpha_list[2]
while i<=n:
    S_sum=(1/(2*i-1))*np.sin(np.radians((2*(2*i-1)*np.pi*alpha)))
    S_sum_list.append(S_sum)
    S=(4/np.pi)*sum(S_sum_list)
    S_list049.append(S)
    i+=1
#Graphs
#S(t;n) vs n (alpha = 0.01)
plt.figure(figsize=(9,5))
plt.plot(n_list,S_list001,c='green')
plt.axhline(f(alpha_list[0]),color="black",linestyle="dashed",label="f(t)")
plt.title("Approximation of the piecewise function f(t) using the corresponding Fourier series S(t;n) with alpha = 0.01")
plt.ylabel("S(t;n)")
plt.xlabel("n")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#S(t;n) vs n (alpha = 0.25)
plt.figure(figsize=(9,5))
plt.plot(n_list,S_list025,c='green')
plt.axhline(f(alpha_list[1]),color="black",linestyle="dashed",label="f(t)")
plt.title("Approximation of the piecewise function f(t) using the corresponding Fourier series S(t;n) with alpha = 0.25")
plt.ylabel("S(t;n)")
plt.xlabel("n")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#S(t;n) vs n (alpha = 0.49)
plt.figure(figsize=(9,5))
plt.plot(n_list,S_list049,c='green')
plt.axhline(f(alpha_list[2]),color="black",linestyle="dashed",label="f(t)")
plt.title("Approximation of the piecewise function f(t) using the corresponding Fourier series S(t;n) with alpha = 0.49")
plt.ylabel("S(t;n)")
plt.xlabel("n")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
print("Program execution time:",time.process_time()-start_time_program,"seconds.")