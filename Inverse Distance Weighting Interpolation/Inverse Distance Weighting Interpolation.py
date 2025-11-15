# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 12:39:11 2025
Inverse Distance Weighting (IDW) interpolation using weighted averages and previously available data.
More information on this website: https://en.wikipedia.org/wiki/Inverse_distance_weighting
@author: Lucas Romero Fernández
"""
import time
import csv
import numpy as np
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
#IDW interpolation method
L=25#Amplitude of the Gaussian Exponential (bigger value is equal to more points considered in the weight matrix, and viceversa, and it has to be between the minimum and maximum distance (or an average, ideally) between data points for good results).
def omega_g(x,y):#Weight function
    d=np.sqrt(x**2+y**2)#Distance between points
    return np.exp(-((d**2)/(2*L**2)))#Gaussian exponential type
Omega_T=np.zeros((len(X_g_matrix_resh),len(value_array)))#Transposed weight matrix
for i in range(0,len(X_g_matrix_resh)):
    for j in range(0,len(value_array)):
        Omega_T[i,j]=omega_g(X_g_matrix_resh[i]-X_array[j],Y_g_matrix_resh[i]-Y_array[j])#Computation of the transposed weight matrix Omega^T
for i in range(0,len(Omega_T)):
    Omega_T[i,:]=Omega_T[i,:]/np.sum(Omega_T[i,:])#Normalization of Omega^T (with the condition that the sum of the elements of each row = 1)
weighted_average_values_array=np.dot(Omega_T,value_array)#Computation of the matrix of weighted average values
weighted_average_values_matrix=np.reshape(weighted_average_values_array,(len(Y_g_array),len(X_g_array)))#To obtain a matrix with the right shape for the plot
#Plot
if abs(np.min(value_array))>abs(np.max(value_array)):#Color normalization
    m=-np.min(value_array)
else:
    m=np.max(value_array)
color_levels=np.arange(-int(m*10)/10,int(m*10)/10+0.1,0.1)#For the data/point colors
norm=mcolors.BoundaryNorm(boundaries=color_levels,ncolors=256)#For the data/point colors
fig,ax=plt.subplots(figsize=(11,7))
graph=ax.contourf(X_g_array,Y_g_array,weighted_average_values_matrix,levels=color_levels,cmap="seismic")
ax.scatter(X_array,Y_array,c=value_array,s=30,norm=norm,cmap="seismic")
ax.set(title=f'IDW interpolation with $L$ = ${L}$',xlabel="$X$ [km]",ylabel="$Y$ [km]")
fig.colorbar(graph)
ax.grid()
fig.set_tight_layout(True)
plt.show()
print("Program execution time:",time.process_time()-start_time_program,"seconds.")
