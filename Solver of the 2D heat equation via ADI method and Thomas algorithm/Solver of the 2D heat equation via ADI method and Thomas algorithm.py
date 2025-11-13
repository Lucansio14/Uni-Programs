# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 12:09:49 2021
Solver of the 2D parabolic diffusion Parcial Derivative Equation (PDE) or the 2D heat equation (with diffusivity constant equal to 1): u_t = u_xx + u_yy using the Alternating Direction Implicit (ADI) method and (slightly modified) Thomas
algorithm in the numerical grid region [x0,xf] in the x-axis and [y0,yf] in the y-axis, with rectangular Dirichlet boundary conditions (BCs) u(t,x0,y) = u(t,x,y0) = u(t,xf,y) = u(t,x,yf) = 0 and initial condition u(t0,x,y) = u0(x,y).
More information on these websites: https://en.wikipedia.org/wiki/Heat_equation
                                    https://en.wikipedia.org/wiki/Alternating-direction_implicit_method
                                    https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
                                    https://www.akademisains.gov.my/asmsj/?mdocs-file=4319
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants, functions, arrays and lists
End_Time_Iter=200#Last time step (dt) iteration
Num_iter_shown=10#Number of time step (dt) iterations represented in the plots (NOT COUNTING t=0dt)
Time_Iter_array=np.linspace(0,End_Time_Iter,Num_iter_shown+1)#Each time step (dt) iteration represented in the plots
x0=0
y0=0
xf=320
yf=200
x_array=np.arange(x0,xf)#Numerical spatial x-array of grid points
y_array=np.arange(y0,yf)#Numerical spatial y-array of grid points
x_grid,y_grid=np.meshgrid(x_array,y_array)#Complete grid
nu=1.5#Grouped ratio between spatial (dx and dy) and time (dt) grid discretizations in both spatial axis (nu_x and nu_y)
theta=0.5#Parameter that performs a weighted average of two other numerical schemes (the FTCS (theta = 0) and BTCS (theta = 1) methods). An special case is obtained when theta = 0.5 (Crank Nicholson scheme).
R=xf-1#Number of 'active' (changing values) x-axis grid points (due to BCs)
S=yf-1#Number of 'active' (changing values) y-axis grid points (due to BCs)
def ReadIC():#Function to read and save the data from the file that contains the initial condition (THE NUMBER OF VALUES IN THE FILE MUST MATCH WITH THE NUMBER OF GRID POINTS)
	file_name="L(200x320).dat"#In this case, the initial condition has the shape of the letter 'L'
	with open(file_name) as file:
		content=file.read().splitlines()
	u_array=np.zeros((yf,xf))
	for i in range(y0,yf):
		j=0
		for x in content[i]:
			u_array[i,j]=x
			j+=1	
	return u_array[::-1,:]
u0_array=ReadIC()
u0_array_graph=u0_array#For the graphs
def Thomas(theta,nu,num_active_points,u_substep_array):#(slightly modified) Thomas or tridiagonal matrix algorithm
    x_array=np.zeros(num_active_points+1)
    e_array=np.zeros(num_active_points)
    f_array=np.zeros(num_active_points)
    a_array=np.full(num_active_points,-theta*nu)
    b_array=np.full(num_active_points,1+2*theta*nu)
    c_array=np.full(num_active_points,-theta*nu)
    for i in range(1,num_active_points):
        e_array[i]=-c_array[i]/(b_array[i]+a_array[i]*e_array[i-1])
        f_array[i]=(u_substep_array[i]-a_array[i]*f_array[i-1])/(b_array[i]+a_array[i]*e_array[i-1])
    x_array[num_active_points-1]=f_array[num_active_points-1]
    for j in range(num_active_points-2,0,-1):
        x_array[j]=f_array[j]+e_array[j]*x_array[j+1]
    return x_array
def ADI1(nu,R,S,u_substep_t_array):#Step 1 (substep 1 in time) of the ADI method (implicit schema in the x-direction (rows) and explicit schema in y-direction (columns))
    u_substep_x_array=np.zeros(R)
    for j in range(1,S):
        for i in range(1,R):
            u_substep_x_array[i]=1/2*nu*u_substep_t_array[j+1,i]+(1-nu)*u_substep_t_array[j,i]+1/2*nu*u_substep_t_array[j-1,i]
        u_substep_t_array[j,:]=Thomas(theta,nu,R,u_substep_x_array)
    return u_substep_t_array
def ADI2(nu,R,S,u_substep_t_array):#Step 2 (substep 2 in time) of the ADI method (implicit schema in the y-direction (columns) and explicit schema in x-direction (rows))
    u_substep_y_array=np.zeros(S)
    for i in range(1,R):
        for j in range(1,S):
            u_substep_y_array[j]=1/2*nu*u_substep_t_array[j,i+1]+(1-nu)*u_substep_t_array[j,i]+1/2*nu*u_substep_t_array[j,i-1]
        u_substep_t_array[:,i]=Thomas(theta,nu,S,u_substep_y_array)
    return u_substep_t_array
#Solution evolution and graphs for each time step iterations
Time_Iter=0
while Time_Iter<=End_Time_Iter:
    u_array=ADI2(nu,R,S,ADI1(nu,R,S,u0_array))
    if Time_Iter in Time_Iter_array:
        #2D (heat map)
        fig,axs=plt.subplots(figsize=(11,7))
        Cont=axs.contourf(np.arange(xf),np.arange(yf),u_array,cmap="hot")
        axs.set(title=f'Diffusion represented in 2D for $t={Time_Iter}$dt',xlabel="$x$",ylabel="$y$")
        fig.colorbar(Cont)
        axs.grid()
        fig.set_tight_layout(True)
        plt.show()
    u0_array=u_array
    Time_Iter+=1
print("Program execution time:",time.process_time()-start_time_program,"seconds.")