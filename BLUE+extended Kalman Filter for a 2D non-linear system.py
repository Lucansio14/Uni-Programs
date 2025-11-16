# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 17:11:13 2025
BLUE("Best Linear Unbiased Estimation") + extended Kalman filter for the case of a 2D non-linear system consisting of the motion (with on-the-spot corrections on its path) of a bicycle (model) & the desired path chosen beforehand ("reality").
More information on these websites: https://fiveable.me/introduction-econometrics/unit-4/linear-unbiased-estimator-blue/study-guide/rXYVI3ihDJ5xd1DA
                                    https://en.wikipedia.org/wiki/Gauss%E2%80%93Markov_theorem
                                    https://en.wikipedia.org/wiki/Extended_Kalman_filter
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
#Definition of general constants, variables, matrices, lists, and functions
t0=0#Start time
T=4#End time
dt=0.01#Time-step size (not too big)
dt_obs=0.1#Time-step size (observations) (has to be multiple of dt and a factor of T) (not too big)
dim=2#Number of dimensions (2D = 2 dimensions)
Num_variables_bic=4#Number of variables of the bicycle model
t_array=np.around(np.arange(t0,T+dt,dt),2)#Values of time
t_obs_array=np.around(np.arange(t0+dt_obs,T+dt_obs,dt_obs),2)#Values of time for observed data (do not start at t0!)
x_0_reality_matrix=np.zeros([Num_variables_bic,1])#State vector x_n (current time value) for the "reality"
x_n_reality_matrix=np.zeros([Num_variables_bic,1])#State vector x_n+1 (next time value) for the "reality"
x_0_bic_matrix=np.zeros([Num_variables_bic,1])#State vector x_n (current time value) for the model
x_n_bic_matrix=np.zeros([Num_variables_bic,1])#State vector x_n+1 (next time value) for the model
L=1#Bicycle wheel diameter
delta=0.10472#Handlebar orientation (in radians: 1 degree is equal to 0.0174533 radians)
sigma=1#Standard deviation of the noise
omega=np.random.uniform(0,sigma)#Velocity noise
L_matrix=np.eye(Num_variables_bic)#Jacobian of the bicycle model matrix L
L_T_matrix=np.transpose(L_matrix)#Transposed jacobian of the bicycle model matrix L^T
H_matrix=np.zeros([dim,Num_variables_bic])#Jacobian of the observation data matrix H
H_T_matrix=np.transpose(H_matrix)#Transposed jacobian of the observation data matrix H^T
P_n_matrix=0.15*np.eye(Num_variables_bic)#Prediction matrix P_n (with initial prediction values)
Q_matrix=np.zeros([Num_variables_bic,Num_variables_bic])#Covariance model errors matrix (background error) Q
obs_points_x_list=[]#Observation or reference points for the position x
obs_points_y_list=[]#Observation or reference points for the position y
x_reality_list=[]#Position x for the "reality"
y_reality_list=[]#Position y for the "reality"
theta_reality_list=[]#Direction angle for the "reality"
v_reality_list=[]#Velocity (module) for the "reality"
x_bic_list=[]#Position x for the bicycle model
y_bic_list=[]#Position y for the bicycle model
theta_bic_list=[]#Direction angle for the bicycle model
v_bic_list=[]#Velocity (module) for the bicycle model
P_n_diag_x_list=[]#Position x values of the prediction matrix P_n
P_n_diag_y_list=[]#Position y values of the prediction matrix P_n
P_n_diag_theta_list=[]#Direction angle values of the prediction matrix P_n
P_n_diag_v_list=[]#Velocity (module) values of the prediction matrix P_n
def Model_Evol(x_0_matrix):#Evolution function of the state vector x_0
    x_n_matrix=np.zeros([Num_variables_bic,1])
    x_n_matrix[0,0]=x_0_matrix[0,0]+x_0_matrix[3,0]*np.cos(x_0_matrix[2,0])*dt
    x_n_matrix[1,0]=x_0_matrix[1,0]+x_0_matrix[3,0]*np.sin(x_0_matrix[2,0])*dt
    x_n_matrix[2,0]=x_0_matrix[2,0]+((x_0_matrix[3,0])/(L))*np.tan(delta)*dt
    x_n_matrix[3,0]=x_0_matrix[3,0]+omega
    return x_n_matrix
def dist_points(x_0_bic_matrix,Obs_points_list_x,Obs_points_list_y,i):#Distance between observation points function
    r=np.sqrt((x_0_bic_matrix[0,0]-Obs_points_list_x[i])**2+(x_0_bic_matrix[1,0]-Obs_points_list_y[i])**2)
    return r
def L_matrix_Evol(x_0_bic_matrix):#Evolution function of the jacobian of the bicycle model matrix L
    L_matrix[0,2]=-x_0_bic_matrix[3,0]*np.sin(x_0_bic_matrix[2,0])*dt
    L_matrix[0,3]=np.cos(x_0_bic_matrix[2,0])*dt
    L_matrix[1,2]=x_0_bic_matrix[3,0]*np.cos(x_0_bic_matrix[2,0])*dt
    L_matrix[1,3]=np.sin(x_0_bic_matrix[2,0])*dt
    L_matrix[2,3]=(1/L)*np.tan(delta)*dt
    return L_matrix
def H_matrix_Evol(x_0_bic_matrix,r_1,r_2,Obs_points_list_x,Obs_points_list_y,i,j):#Evolution function of the jacobian of the observation data matrix H
    H_matrix[0,0]=(-x_0_bic_matrix[0,0]+Obs_points_list_x[i])/(r_1)
    H_matrix[0,1]=(-x_0_bic_matrix[1,0]+Obs_points_list_y[i])/(r_1)
    H_matrix[1,0]=(-x_0_bic_matrix[0,0]+Obs_points_list_x[j])/(r_2)
    H_matrix[1,1]=(-x_0_bic_matrix[1,0]+Obs_points_list_y[j])/(r_2)
    return H_matrix
def P_evo(P_matrix):#Evolution function of the prediction matrix P_f
    P_f_matrix=np.dot(np.dot(L_matrix,P_matrix),L_T_matrix)+Q_matrix
    return P_f_matrix
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Initial conditions (or ICs) for every system variable (for the "reality")
x_0_reality_matrix[0,0]=0.25#For x ("reality")
x_0_reality_matrix[1,0]=0.25#For y ("reality")
x_0_reality_matrix[2,0]=0.15#For theta (in radians) ("reality")
x_0_reality_matrix[3,0]=0.75#For v ("reality")
#Computation of the "reality" values and selecting the "observed" values
for n in t_array:#For the "reality" values
    x_n_reality_matrix=Model_Evol(x_0_reality_matrix)
    x_reality_list.append(x_n_reality_matrix[0,0])
    y_reality_list.append(x_n_reality_matrix[1,0])
    theta_reality_list.append(x_n_reality_matrix[2,0])
    v_reality_list.append(x_n_reality_matrix[3,0])
    if n in t_obs_array:#For the "observed" values
        obs_points_x_list.append(x_n_reality_matrix[0,0])
        obs_points_y_list.append(x_n_reality_matrix[1,0])
    x_0_reality_matrix=x_n_reality_matrix#Update state vector
#Computation of the observation noise and the noise matrix, and the covariance matrix
noise_x_array=np.random.uniform(-sigma,sigma,len(obs_points_x_list))#Random noise in position x
noise_y_array=np.random.uniform(-sigma,sigma,len(obs_points_y_list))#Random noise in position y
obs_points_x_list=obs_points_x_list+noise_x_array#Addition of the noise to the "observed" values
obs_points_y_list=obs_points_y_list+noise_y_array#Addition of the noise to the "observed" values
R_matrix=(sigma**2)*np.eye(dim)#Covariance matrix R (Data reliability)
#Solution/Model evolution
#Initial conditions (or ICs) for every system variable (for the model)
x_0_bic_matrix[0,0]=0#For x
x_0_bic_matrix[1,0]=0#For y
x_0_bic_matrix[2,0]=0#For theta (in radians)
x_0_bic_matrix[3,0]=0.5#For v
y_obs_array=np.array([[0],[0]])#Observation data (only for positions x and y) vector
i=0#Loop index
i_final=len(obs_points_x_list)-1#Final iteration loop index
for n in t_array:#Model
    if n in t_obs_array:#For the "observed" values
        if i==i_final:#This last observation point always causes accuracy issues...
            r_1=dist_points(x_0_bic_matrix,obs_points_x_list,obs_points_y_list,i)
            r_2=dist_points(x_0_bic_matrix,obs_points_x_list,obs_points_y_list,i)
            H_matrix=H_matrix_Evol(x_0_bic_matrix,r_1,r_2,obs_points_x_list,obs_points_y_list,i,i)
        else:
            r_1=dist_points(x_0_bic_matrix,obs_points_x_list,obs_points_y_list,i)
            r_2=dist_points(x_0_bic_matrix,obs_points_x_list,obs_points_y_list,i+1)
            H_matrix=H_matrix_Evol(x_0_bic_matrix,r_1,r_2,obs_points_x_list,obs_points_y_list,i,i+1)
        y_obs_array[0]=r_1
        y_obs_array[1]=r_2
        H_T_matrix=np.transpose(H_matrix)
        i+=1
        x_b_matrix=Model_Evol(x_0_bic_matrix)#State vector without Kalman correction
        K_matrix=np.dot(P_n_matrix,np.dot(H_T_matrix,np.linalg.inv(R_matrix+np.dot(H_matrix,np.dot(P_n_matrix,H_T_matrix)))))#Kalman/Weight matrix
        P_n_matrix=np.dot((np.eye(len(K_matrix))-np.dot(K_matrix,H_matrix)),P_n_matrix)
        x_n_bic_matrix=x_b_matrix+np.dot(K_matrix,(y_obs_array-np.dot(H_matrix,x_b_matrix)))
    else:
        x_n_bic_matrix=Model_Evol(x_0_bic_matrix)
    x_bic_list.append(x_n_bic_matrix[0,0])
    y_bic_list.append(x_n_bic_matrix[1,0])
    theta_bic_list.append(x_n_bic_matrix[2,0])
    v_bic_list.append(x_n_bic_matrix[3,0])
    P_n_diag_x_list.append(P_n_matrix[0,0])
    P_n_diag_y_list.append(P_n_matrix[1,1])
    P_n_diag_theta_list.append(P_n_matrix[2,2])
    P_n_diag_v_list.append(P_n_matrix[3,3])
    L_matrix=L_matrix_Evol(x_0_bic_matrix)
    L_T_matrix=np.transpose(L_matrix)
    P_n_matrix=P_evo(P_n_matrix)
    x_0_bic_matrix=x_n_bic_matrix#Update state vector
#Plots
#x vs. t
plt.figure(figsize=(11,7))
plt.xlim(t0-dt_obs,T+dt_obs)
plt.plot(t_array,x_bic_list,label="Model")
plt.plot(t_array,x_reality_list,label="'Reality'")
plt.scatter(t_obs_array,obs_points_x_list,c="black",label="Observations")
plt.xlabel("$t$")
plt.ylabel("$x$")
plt.grid()
plt.title("BLUE + extended Kalman filter for bicycle system: $x$ vs $t$")
plt.legend()
plt.tight_layout()
plt.show()
#y vs. t
plt.figure(figsize=(11,7))
plt.xlim(t0-dt_obs,T+dt_obs)
plt.plot(t_array,y_bic_list,label="Model")
plt.plot(t_array,y_reality_list,label="'Reality'")
plt.scatter(t_obs_array,obs_points_y_list,c="black",label="Observations")
plt.xlabel("$t$")
plt.ylabel("$y$")
plt.grid()
plt.title("BLUE + extended Kalman filter for bicycle system: $y$ vs $t$")
plt.legend()
plt.tight_layout()
plt.show()
#theta vs. t
plt.figure(figsize=(11,7))
plt.xlim(t0-dt_obs,T+dt_obs)
plt.plot(t_array,theta_bic_list,label="Model")
plt.plot(t_array,theta_reality_list,label="'Reality'")
plt.xlabel("$t$")
plt.ylabel(r"$\theta$")
plt.grid()
plt.title(r"BLUE + extended Kalman filter for bicycle system: $\theta$ vs $t$")
plt.legend()
plt.tight_layout()
plt.show()
#|v| vs. t
plt.figure(figsize=(11,7))
plt.xlim(t0-dt_obs,T+dt_obs)
plt.plot(t_array,v_bic_list,label="Model")
plt.plot(t_array,v_reality_list,label="'Reality'")
plt.xlabel("$t$")
plt.ylabel("$|v|$")
plt.grid()
plt.title("BLUE + extended Kalman filter for bicycle system: $|v|$ vs $t$")
plt.legend()
plt.tight_layout()
plt.show()
#Trajectory
plt.figure(figsize=(11,7))
plt.plot(x_bic_list,y_bic_list,label="Model")
plt.plot(x_reality_list,y_reality_list,label="'Reality'")
plt.scatter(obs_points_x_list,obs_points_y_list,c="black",label="Observations")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.grid()
plt.title("BLUE + extended Kalman filter for bicycle system: $x$ vs $y$")
plt.legend()
plt.tight_layout()
plt.show()
#P_n vs. t
plt.figure(figsize=(11,7))
plt.xlim(t0-dt_obs,T+dt_obs)
plt.plot(t_array,P_n_diag_x_list,label="x")
plt.plot(t_array,P_n_diag_y_list,"--",label="y")
plt.plot(t_array,P_n_diag_theta_list,label=r"$\theta$")
plt.plot(t_array,P_n_diag_v_list,"--",label="v")
plt.xlabel("$t$")
plt.ylabel("$P_{n}$")
plt.grid()
plt.title("BLUE + extended Kalman filter for bicycle system: $P_{n}$ vs $t$")
plt.legend()
plt.tight_layout()
plt.show()
print("Program execution time:",time.process_time()-start_time_program,"seconds.")