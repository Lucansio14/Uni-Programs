# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 09:46:14 2025
Polynomial interpolation using a polynomial of arbitrary order n and previously available or observed data.
More information on this website: https://en.wikipedia.org/wiki/Polynomial_interpolation
@author: Lucas Romero Fernández
"""
import time
import csv
import numpy as np
from sympy import symbols
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
#Extraction and transformation of data from file
X_list=[]#Positions x of the data (Units: [km])
Y_list=[]#Positions y of the data (Units: [km])
value_list=[]#Values of the data
file_name="data_t0_nobs_30.dat"
with open(file_name) as file:
    data=csv.reader(file,delimiter=" ")
    for row in data:
        X_list.append(row[0])
        Y_list.append(row[1])
        value_list.append(row[2])
X_list.pop(0)#Neglect headings
Y_list.pop(0)#Neglect headings
value_list.pop(0)#Neglect headings
X_array=np.array(X_list,dtype="float64")#The length of the data arrays must coincide between each other
Y_array=np.array(Y_list,dtype="float64")#The length of the data arrays must coincide between each other
value_T_array=np.array(value_list,dtype="float64")#The length of the data arrays must coincide between each other
value_array=np.transpose(value_T_array)
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Construction of the grid
num_points_x_grid=20#The number of points in the grid do not significantly affect the results
num_points_y_grid=25#The number of points in the grid do not significantly affect the results
X_g_array=np.linspace(min(X_array),max(X_array),num_points_x_grid)
Y_g_array=np.linspace(min(Y_array),max(Y_array),num_points_y_grid)
X_g_matrix,Y_g_matrix=np.meshgrid(X_g_array,Y_g_array)
X_g_matrix_resh=X_g_matrix.ravel(order="C")#Reshaping of the grid for the computations
Y_g_matrix_resh=Y_g_matrix.ravel(order="C")#Reshaping of the grid for the computations
#Polynomial interpolation method using a polynomial of this form: hat_phi_g = alpha_0*1 + alpha_1*x + alpha_2*x^2 + alpha_3*y + alpha_4*x*y + alpha_5*y^2 + ...
x=symbols("x")#To visualize the terms that are considered
y=symbols("y")#To visualize the terms that are considered
n=3#Integer polynomial order(>0)
m=0#Loop parameter
columns=((n+1)*(n+2))/2
P_obs_T_matrix=np.zeros((len(X_array),int(columns)))#Matrix of polynomial terms for the observed data values
P_g_T_matrix=np.zeros((len(X_g_array)*len(Y_g_array),int(columns)))#Matrix of polynomial terms for the grid values
Terms_list=[]#Considered terms of the polynomial
for j in range(0,n+1):#Computation of the transposed matrices P^T_obs and P^T_g
    for i in range(0,n-j+1):
        P_obs_T_matrix[:,int(m)]=X_array**i*(Y_array**j)
        P_g_T_matrix[:,int(m)]=X_g_matrix_resh**i*(Y_g_matrix_resh**j)
        Terms_list.append(x**i*y**j)
        m+=1
print("Terms considered:",Terms_list)
P_obs_matrix=np.transpose(P_obs_T_matrix)
alpha_array=np.dot(np.dot(np.linalg.inv((np.dot(P_obs_matrix,P_obs_T_matrix))),P_obs_matrix),value_array)#Computation of the polynomial coefficients (alpha)
inter_values_array=np.dot(P_g_T_matrix,alpha_array)#Computation of the matrix of interpolated values
point_values_array=np.dot(P_obs_T_matrix,alpha_array)#Computation of the matrix of point/data values
inter_values_matrix=np.reshape(inter_values_array,(len(Y_g_array),len(X_g_array)))#To obtain a matrix with the right shape for the plot
#Plot and error results
if abs(np.min(value_array))>abs(np.max(value_array)):#Color normalization
    m=-np.min(value_array)
else:
    m=np.max(value_array)
color_levels=np.arange(-int(m*10)/10,int(m*10)/10+0.1,0.1)#For the data point colors
norm=mcolors.BoundaryNorm(boundaries=color_levels,ncolors=256)#For the data point colors
fig,ax=plt.subplots(figsize=(11,7))
graph=ax.contourf(X_g_array,Y_g_array,inter_values_matrix,levels=color_levels,cmap="seismic")
ax.scatter(X_array,Y_array,c=value_array,s=30,norm=norm,cmap="seismic")
ax.set(title=f'Polynomical Interpolation of order $n$ = ${n}$',xlabel="$X$ [km]",ylabel="$Y$ [km]")
fig.colorbar(graph)
ax.grid()
fig.set_tight_layout(True)
plt.show()
Error_array=abs(value_array-point_values_array)
print("Error of each data point between measurement and interpolation:")
print(Error_array)
print("Program execution time:",time.process_time()-start_time_program,"seconds.")
