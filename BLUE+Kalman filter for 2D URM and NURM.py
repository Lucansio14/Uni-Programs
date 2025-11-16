# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 17:11:13 2025
BLUE("Best Linear Unbiased Estimation") + Kalman filter for the case of 2D URM ("Uniform Rectilinear Motion") (model) & NURM ("Non-Uniform Rectilinear Motion") with constant acceleration ("reality") using a variational method (linear version).
More information on these websites: https://fiveable.me/introduction-econometrics/unit-4/linear-unbiased-estimator-blue/study-guide/rXYVI3ihDJ5xd1DA
                                    https://en.wikipedia.org/wiki/Gauss%E2%80%93Markov_theorem
                                    https://en.wikipedia.org/wiki/Kalman_filter
                                    https://en.wikipedia.org/wiki/Linear_motion
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
#Definition of general constants, variables, matrices and lists
t0=0#Start time
T=4#End time
dt=0.01#Time-step size (not too big)
dt_obs=0.1#Time-step size (observations) (has to be multiple of dt and a factor of T) (not too big)
dim=2#Number of dimensions (2D = 2 dimensions)
Num_variables_URM=4#Number of variables considered for 2D URM (model)
Num_variables_NURM=6#Number of variables considered for 2D NURM ("reality")
t_array=np.around(np.arange(t0,T+dt,dt),2)#Values of time
t_obs_array=np.around(np.arange(t0+dt_obs,T+dt_obs,dt_obs),2)#Values of time for observed data (do not start at t0!)
x_0_URM_matrix=np.zeros([Num_variables_URM,1])#State vector x_n (current time value) for URM
x_n_URM_matrix=np.zeros([Num_variables_URM,1])#State vector x_n+1 (next time value) for URM
M_URM_matrix=np.array([[1,dt,0,0],[0,1,0,0],[0,0,1,dt],[0,0,0,1]])#Linear model operator M for 2D URM
M_URM_T_matrix=np.transpose(M_URM_matrix)#Transposed linear model operator M^T for URM
x_0_NURM_matrix=np.zeros([Num_variables_NURM,1])#State vector x_n (current time value) for NURM
x_n_NURM_matrix=np.zeros([Num_variables_NURM,1])#State vector x_n+1 (next time value) for NURM
M_NURM_matrix=np.array([[1,0,dt,0,0.5*(dt)**2,0],[0,1,0,dt,0,0.5*(dt)**2],[0,0,1,0,dt,0],[0,0,0,1,0,dt],[0,0,0,0,1,0],[0,0,0,0,0,1]])#Linear model operator M for NURM ("reality")
dxURM_list=[]#Position x for URM
dyURM_list=[]#Position y for URM
uURM_list=[]#Velocity in the x-direction for URM
vURM_list=[]#Velocity in the y-direction for URM
dxNURM_list=[]#Position x for NURM
dyNURM_list=[]#Position y for NURM
uNURM_list=[]#Velocity in the x-direction for NURM
vNURM_list=[]#Velocity in the y-direction for NURM
ax_list=[]#Acceleration in the x-direction for NURM
ay_list=[]#Acceleration in the y-direction for NURM
dxNURM_sel_list=[]#Selected values of position x for the observations
dyNURM_sel_list=[]#Selected values of position y for the observations
P_n_diag_x_list=[]#Position x values of the prediction matrix P_n
P_n_diag_y_list=[]#Position y values of the prediction matrix P_n
P_n_diag_u_list=[]#Velocity in the x-direction values of the prediction matrix P_n
P_n_diag_v_list=[]#Velocity in the y-direction values of the prediction matrix P_n
def P_evo(P_matrix):#Evolution function of the prediction matrix P_f
    P_f_matrix=np.dot(np.dot(M_URM_matrix,P_matrix),M_URM_T_matrix)+Q_matrix
    return P_f_matrix
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Initial conditions (or ICs) for every system variable
x_0_URM_matrix[0,0]=1#For dx (URM)
x_0_URM_matrix[1,0]=1#For u (URM)
x_0_URM_matrix[2,0]=2#For dy (URM)
x_0_URM_matrix[3,0]=0.5#For v (URM)
x_0_NURM_matrix[0,0]=1#For dx (NURM)
x_0_NURM_matrix[1,0]=2#For dy (NURM)
x_0_NURM_matrix[2,0]=1#For u (NURM)
x_0_NURM_matrix[3,0]=0.5#For v (NURM)
x_0_NURM_matrix[4,0]=2#For a_x (NURM)
x_0_NURM_matrix[5,0]=1#For a_y (NURM)
y_obs_array=np.array([[0],[0]])#Observation data (only for positions x and y) vector
P_n_matrix=0.15*np.eye(Num_variables_URM)#Prediction matrix P_n (with initial prediction values)
Q_matrix=np.zeros((Num_variables_URM,Num_variables_URM))#Covariance model errors matrix Q
H_matrix=np.array([[1,0,0,0],[0,0,1,0]])#Observation data matrix H
H_T_matrix=np.transpose(H_matrix)#Transposed observation data matrix H^T
#Computation of the "reality" values and selecting the "observed" values
for n in t_array:#For the "reality" values (NURM)
    x_n_NURM_matrix=np.dot(M_NURM_matrix,x_0_NURM_matrix)
    dxNURM_list.append(x_n_NURM_matrix[0,0])
    dyNURM_list.append(x_n_NURM_matrix[1,0])
    uNURM_list.append(x_n_NURM_matrix[2,0])
    vNURM_list.append(x_n_NURM_matrix[3,0])
    ax_list.append(x_n_NURM_matrix[4,0])
    ay_list.append(x_n_NURM_matrix[5,0])
    if n in t_obs_array:#For the "observed" values
        dxNURM_sel_list.append(x_n_NURM_matrix[0,0])
        dyNURM_sel_list.append(x_n_NURM_matrix[1,0])
    x_0_NURM_matrix=x_n_NURM_matrix#Update state vector
#Computation of the observation noise, the noise matrix, and the covariance matrix
sigma=1#Standard deviation of the noise
noise_x_array=np.random.uniform(-sigma,sigma,len(dxNURM_sel_list))#Random noise in position x
noise_y_array=np.random.uniform(-sigma,sigma,len(dyNURM_sel_list))#Random noise in position y
dxNURM_sel_list=dxNURM_sel_list+noise_x_array#Addition of the noise to the "observed" values
dyNURM_sel_list=dyNURM_sel_list+noise_y_array#Addition of the noise to the "observed" values
R_matrix=(sigma**2)*np.eye(dim)#Covariance matrix R
#Solution/Model evolution
i=0#Loop index
for n in t_array:#Model (URM)
    if n in t_obs_array:#For the "observed" values
        y_obs_array[0]=dxNURM_sel_list[i]
        y_obs_array[1]=dyNURM_sel_list[i]
        i+=1
        x_b_matrix=np.dot(M_URM_matrix,x_0_URM_matrix)#State vector without Kalman correction
        K_matrix=np.dot(P_n_matrix,np.dot(H_T_matrix,np.linalg.inv(R_matrix+np.dot(H_matrix,np.dot(P_n_matrix,H_T_matrix)))))#Kalman/Weight matrix
        P_n_matrix=np.dot((np.eye(len(K_matrix))-np.dot(K_matrix,H_matrix)),P_n_matrix)
        x_n_URM_matrix=x_b_matrix+np.dot(K_matrix,(y_obs_array-np.dot(H_matrix,x_b_matrix)))
    else:
        x_n_URM_matrix=np.dot(M_URM_matrix,x_0_URM_matrix)
    dxURM_list.append(x_n_URM_matrix[0,0])
    dyURM_list.append(x_n_URM_matrix[2,0])
    uURM_list.append(x_n_URM_matrix[1,0])
    vURM_list.append(x_n_URM_matrix[3,0])
    P_n_diag_x_list.append(P_n_matrix[0,0])
    P_n_diag_u_list.append(P_n_matrix[1,1])
    P_n_diag_y_list.append(P_n_matrix[2,2])
    P_n_diag_v_list.append(P_n_matrix[3,3])
    P_n_matrix=P_evo(P_n_matrix)
    x_0_URM_matrix=x_n_URM_matrix#Update state vector
#Plots
#x vs. t
plt.figure(figsize=(11,7))
plt.xlim(t0-dt_obs,T+dt_obs)
plt.plot(t_array,dxURM_list,label="Model")
plt.plot(t_array,dxNURM_list,label="'Reality'")
plt.scatter(t_obs_array,dxNURM_sel_list,c="black",label="Observations")
plt.xlabel("$t$")
plt.ylabel("$x$")
plt.grid()
plt.title("BLUE + Kalman filter method for linear motion: $x$ vs $t$")
plt.legend()
plt.tight_layout()
plt.show()
#y vs. t
plt.figure(figsize=(11,7))
plt.xlim(t0-dt_obs,T+dt_obs)
plt.plot(t_array,dyURM_list,label="Model")
plt.plot(t_array,dyNURM_list,label="'Reality'")
plt.scatter(t_obs_array,dyNURM_sel_list,c="black",label="Observations")
plt.xlabel("$t$")
plt.ylabel("$y$")
plt.grid()
plt.title("BLUE + Kalman filter method for linear motion: $y$ vs $t$")
plt.legend()
plt.tight_layout()
plt.show()
#u vs. t
plt.figure(figsize=(11,7))
plt.xlim(t0-dt_obs,T+dt_obs)
plt.plot(t_array,uURM_list,label="Model")
plt.plot(t_array,uNURM_list,label="'Reality'")
plt.xlabel("$t$")
plt.ylabel("$u$")
plt.grid()
plt.title("BLUE + Kalman filter method for linear motion: $u$ vs $t$")
plt.legend()
plt.tight_layout()
plt.show()
#v vs. t
plt.figure(figsize=(11,7))
plt.xlim(t0-dt_obs,T+dt_obs)
plt.plot(t_array,vURM_list,label="Model")
plt.plot(t_array,vNURM_list,label="'Reality'")
plt.xlabel("$t$")
plt.ylabel("$v$")
plt.grid()
plt.title("BLUE + Kalman filter method for linear motion: $v$ vs $t$")
plt.legend()
plt.tight_layout()
plt.show()
#Trajectory
plt.figure(figsize=(11,7))
plt.plot(dxURM_list,dyURM_list,label="Model")
plt.plot(dxNURM_list,dyNURM_list,label="'Reality'")
plt.scatter(dxNURM_sel_list,dyNURM_sel_list,c="black",label="Observations")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.grid()
plt.title("BLUE + Kalman filter method for linear motion: $x$ vs $y$")
plt.legend()
plt.tight_layout()
plt.show()
#P_n vs. t
plt.figure(figsize=(11,7))
plt.xlim(t0-dt_obs,T+dt_obs)
plt.plot(t_array,P_n_diag_x_list,label="x")
plt.plot(t_array,P_n_diag_y_list,"--",label="y")
plt.plot(t_array,P_n_diag_u_list,label="u")
plt.plot(t_array,P_n_diag_v_list,"--",label="v")
plt.xlabel("$t$")
plt.ylabel("$P_{n}$")
plt.grid()
plt.title("BLUE + Kalman filter method for linear motion: $P_{n}$ vs $t$")
plt.legend()
plt.tight_layout()
plt.show()
print("Program execution time:",time.process_time()-start_time_program,"seconds.")