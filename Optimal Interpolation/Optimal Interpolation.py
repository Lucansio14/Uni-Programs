# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 09:46:14 2025
Optimal interpolation using observations and background data.
More information on these websites: https://www.atmosp.physics.utoronto.ca/PHY2509/ch3.pdf
                                    https://www.researchgate.net/publication/265626429_Introduction_to_Optimal_Interpolation_and_Variational_Analysis#pf5
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
#Computation of background or first guess values (corresponding to a polynomial interpolation of order n, in this case)
x=symbols("x")#To visualize the terms that are considered
y=symbols("y")#To visualize the terms that are considered
n=2#Integer polynomial order(>0)
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
print("Polynomial terms considered for the background data:",Terms_list)
P_obs_matrix=np.transpose(P_obs_T_matrix)#Computation of the matrix P_obs
alpha_array=np.dot(np.dot(np.linalg.inv((np.dot(P_obs_matrix,P_obs_T_matrix))),P_obs_matrix),value_array)#Computation of the polynomial coefficients (alpha)
inter_values_back_array=np.dot(P_g_T_matrix,alpha_array)#Computation of the matrix of interpolated values
point_values_back_array=np.dot(P_obs_T_matrix,alpha_array)#Computation of the matrix of point values
#Optimal Interpolation (Blend between Inverse Distance Weighting (IDW) interpolation and polynomial interpolation) with background values
L=30#Amplitude of the Gaussian Exponential (bigger value is equal to more points considered in the weight matrix, and viceversa, and it has to be between the minimum and maximum distance (or an average, ideally) between data points for good results).
eta=0.5#Mean error variance
sigma=2.5#Variance of grid values
Coeff_Noise_Signal=(eta**2)/(sigma**2)#Coefficient Noise-Signal \eta^2/\sigma^2 (how unreliable is the data given)
def rho_g(x,y):#Weight function
    d=np.sqrt(x**2+y**2)#Distance between points
    return np.exp(-((d**2)/(2*L**2)))#Gaussian exponential correlation function type
T_matrix=np.zeros((len(X_array),len(X_array)))#Correlation matrix
rho_T_matrix=np.zeros((len(X_g_matrix_resh),len(X_array)))#Observational transposed covariances matrix
rho_T_Point_matrix=np.zeros((len(X_array),len(X_array)))#Observational transposed covariances matrix for the point/data values
for i in range(0,len(X_g_matrix_resh)):
    for j in range(0,len(X_array)):
        rho_T_matrix[i,j]=rho_g(X_g_matrix_resh[i]-X_array[j],Y_g_matrix_resh[i]-Y_array[j])#Computation of the transposed correlation matrix rho^T_g
for i in range(0,len(X_array)):
    for j in range(0,len(X_array)):
        rho_T_Point_matrix[i,j]=rho_g(X_array[i]-X_array[j],Y_array[i]-Y_array[j])#Computation of the transposed correlation matrix rho^T_g for the point/data values
for i in range(0,len(X_array)):
    for j in range(0,len(X_array)):
        if i==j:
            T_matrix[i,j]=1+Coeff_Noise_Signal#Computation of the observational covariances matrix T
        else:
            T_matrix[i,j]=rho_g(X_array[i]-X_array[j],Y_array[i]-Y_array[j])#Computation of the observational covariances matrix T
rho_matrix=np.transpose(rho_T_matrix)#Computation of the correlation matrix rho
inter_values_array=inter_values_back_array+np.dot(np.dot(rho_T_matrix,np.linalg.inv(T_matrix)),(value_array-point_values_back_array))#Computation of the matrix of interpolated values, compared to the background values
point_values_array=value_array+np.dot(np.dot(rho_T_Point_matrix,np.linalg.inv(T_matrix)),(value_array-point_values_back_array))#Computation of the matrix of point/data values
error_analysis_matrix=(sigma**2)*(np.eye(len(np.dot(np.dot(rho_T_matrix,np.linalg.inv(T_matrix)),rho_matrix)))-np.dot(np.dot(rho_T_matrix,np.linalg.inv(T_matrix)),rho_matrix))#Computation of the error matrix
inter_values_matrix=np.reshape(inter_values_array,(len(Y_g_array),len(X_g_array)))#To obtain a matrix with the right shape for the plots
error_analysis_diag=np.diag(error_analysis_matrix)
error_analysis_matrix_resh=np.reshape(error_analysis_diag,(len(Y_g_array),len(X_g_array)))#To obtain a matrix with the right shape for the plots
error_point_array=abs(value_array-point_values_array)#Computation of the error matrix for the point/data values
#Plots
if abs(np.min(value_array))>abs(np.max(value_array)):#Color normalization
    m=-np.min(value_array)
else:
    m=np.max(value_array)
#Interpolation results
color_levels=np.arange(-int(m*10)/10,int(m*10)/10+0.1,0.1)#For the data/point colors
norm=mcolors.BoundaryNorm(boundaries=color_levels,ncolors=256)#For the data/point colors
fig,ax=plt.subplots(figsize=(11,7))
graph=ax.contourf(X_g_array,Y_g_array,inter_values_matrix,levels=color_levels,cmap="seismic")
ax.scatter(X_array,Y_array,c=point_values_array,s=30,norm=norm,cmap="seismic")
ax.set(title=f'Optimal interpolation with ($L$,$\eta^2/\sigma^2$) = (${L}$,${Coeff_Noise_Signal}$)',xlabel="$X$ [km]",ylabel="$Y$ [km]")
fig.colorbar(graph)
ax.grid()
fig.set_tight_layout(True)
plt.show()
#Error analysis
color_levels=np.arange(round(min(error_analysis_diag),0),round(max(error_analysis_diag),0),0.001)#For the data/point colors
norm=mcolors.BoundaryNorm(boundaries=color_levels,ncolors=6000)#For the data/point colors
fig,ax=plt.subplots(figsize=(11,7))
graph=ax.contourf(X_g_array,Y_g_array,error_analysis_matrix_resh,levels=color_levels,cmap="seismic")
ax.scatter(X_array,Y_array,c=error_point_array,s=30,norm=norm,cmap="seismic")
ax.set(title=f'Error analysis for the optimal interpolation with ($L$,$\eta^2/\sigma^2$) = (${L}$,${Coeff_Noise_Signal}$)',xlabel="$X$ [km]",ylabel="$Y$ [km]")
fig.colorbar(graph)
ax.grid()
fig.set_tight_layout(True)
plt.show()
print("Program execution time:",time.process_time()-start_time_program,"seconds.")