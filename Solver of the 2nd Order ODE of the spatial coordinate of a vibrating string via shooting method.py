# -*- coding: utf-8 -*-
"""
Created on Frid Oct 10 17:58:23 2025
Solver of the 2nd Order Ordinary Differential Equation (ODE) of the spatial coordinate of a vibrating string: dy(x,y) = v(x,y), d2y(x,y) = dv(x,y) = yk^2; using the shooting method based on the improved Euler (or Heun's)
and bisection methods on several roots (e.g. eigenvalues k) of the vectors/eigenvalues equation in the region [x0,xf] in the x-axis with initial and final conditions y0 = y(x0) and yf = y(xf), respectively.                                                                                                  
More information on these websites: https://en.wikipedia.org/wiki/Shooting_method
                                    https://en.wikipedia.org/wiki/Heun%27s_method
                                    https://en.wikipedia.org/wiki/Bisection_method
                                    https://uomustansiriyah.edu.iq/media/lectures/5/5_2020_05_18!01_04_40_AM.pdf
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants, functions and lists
i=0#Counter for the mode number
h=0.0001#Step size of the Heun's method
prec=10**(-6)#Precision to be considered in the equation's solution
x0=0
xf=1
y0=0
yf=0
x_array=np.linspace(x0,xf,int((xf-x0)/h)+1)
k_res_list=[]#Resulting eigenvalues list
error_k_list=[]#Error list between numerical and analytical eigenvalues
dy=lambda x,y,dy:dy
d2y=lambda x,y,dy:-y*k**2
def Impr_Euler(x,y,dy,x_plus,h,v,dv):#Improved Euler (or Heun's) method for two variables
    #x_plus: Following value in the x-axis
    y1=y+h*v(x,y,dy)
    dy1=dy+h*dv(x,y,dy)
    y2=v(x_plus,y1,dy1)
    dy2=dv(x_plus,y1,dy1)
    y3=(v(x,y,dy)+y2)/2
    dy3=(dv(x,y,dy)+dy2)/2
    return y+h*y3,dy+h*dy3
#Procedure
for k in [3,7,9,13,16]:#List of initial test eigenvalues
    y_array=np.zeros(len(x_array))
    y_array[0]=y0
    k_anal=(i+1)*np.pi#Analytical eigenvalue
    #Preparation for the graph with all solution attempts for the specific initial test eigenvalue
    fig_k,ax_k=plt.subplots(figsize=(9,5))
    plt.xlim(x0-0.03,xf+0.03)
    plt.title("Solution attempts with initial test eigenvalue $k = %f$" % k)
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.grid()
    fig_k.set_tight_layout(True)
    while True:
        dy_array=np.zeros(len(x_array))
        dy_array[0]=k
        for n in range(0,len(x_array)-1):#Solution attempt
            y_and_dy_new=Impr_Euler(x_array[n],y_array[n],dy_array[n],x_array[n+1],h,dy,d2y)
            y_array[n+1]=y_and_dy_new[0]
            dy_array[n+1]=y_and_dy_new[1]
        plt.plot(x_array,y_array,c="red")#Plot the solution attempt
        k=(np.sqrt(abs(-d2y(x_array[-1],abs(y_array[-1]-yf),dy_array[-1])/(abs(y_array[-1]-yf))))+(i+1)*np.pi)*0.5#New eigenvalue obtained using the bisection method with yf (exclusive to the ODE of the spatial coordinate of a vibrating string)
        error_k=abs(k-k_anal)
        #Results for the specific initial test eigenvalue
        if error_k<prec:
            plt.show()#Show graph with all the solution attempts for the specific initial test eigenvalue
            y_anal=np.sin(np.pi*(i+1)*x_array)#Analytical solution
            k_res_list.append(k)
            error_k_list.append(error_k)
            #Graph for the comparison of the analytical and numerical solutions
            plt.figure(figsize=(9,5))
            plt.xlim(x0-0.03,xf+0.03)
            plt.plot(x_array,y_array,"r.-",label="Numerical")
            plt.plot(x_array,y_anal,c="green",label="Analytical")
            plt.xlabel("$x$")
            plt.ylabel("$y$")
            plt.title("Comparison between solutions of the equation for final eigenvalue $k = %f$" % k)
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.show()
            i+=1#Next mode number
            break
print("Obtained numerical eigenvalues:",k_res_list)
print("Error between numerical and analytical eigenvalues:",error_k_list)
print("Program execution time:",time.process_time()-start_time_program,"seconds.")