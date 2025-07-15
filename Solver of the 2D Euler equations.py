# -*- coding: utf-8 -*-
"""
Created on Tue Jun 9 20:19:17 2025
Solver of the 2D Euler equations: A basic "2D" (initial perturbations only in two directions (x and y), which means derivatives in other directions (z) equal to zero or no spatial variations in other directions) non linear equation system in the conservative and compressible form for pressure (p), density (\rho), velocity (v), and total energy (E), used primarly in hydrodynamics, with RK4 (time part), FD4 (spatial part) and Kreiss-Oliger dissipation filter.
The corresponding equations (considering mostly adimensional variables, and neglecting gravitational effects due to relatively small spatial domain):
    \partial_{t}\rho = -\div{\vec{v}\rho} (Continuity equation)
    \partial_{t}(\vec{v}\rho) = -(1/\rho)\grad{p} - \vec{v} · \grad{\vec{v}} (Momentum equation, because momentum (m) = density (\rho) * velocity (v) (this will be mayorly used in the computations))
    \partial_{t}p = -p\gamma\div{\vec{v}} - \vec{v} · \grad{p} ((Adiabatic) pressure equation)
    E = p/(\gamma - 1) + 0.5\vec{v}^2\rho (Ideal equation of state)
(where \gamma is the adiabatic index and, importantly for later, v_s = \sqrt{(p\gamma)/(\rho)} is the adiabatic sound velocity)
More information about these equations in this website: https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)
More information about CFL criteria for non linear equations in this website: https://scicomp.stackexchange.com/questions/31456/cfl-equation-for-non-linear-equation

Context about the computed solution (with initial perturbations in the velocity):
Without the use of more sophisticated methods (like shock capturing or hyperdifusion), the non linear case (not Amp_v,Amp_p<<1, such as Amp_v = 1 and Amp_p = 1) solutions are "corrupted" in some way by artifact spikes. The reason of this issue is completely numerical. For that, only the linear case (Amp_v,Amp_p<<1, such as Amp_v = 0.1 and Amp_p = 0.1) solution is correctly obtainable with this code.
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
import os
#Definition of general constants, variables, lists, arrays and functions
dx=0.02#Space-step size in the x coordinate
dy=0.02#Space-step size in the y coordinate
x0=-1#Start of the space domain in the x coordinate
xf=1#End of the space domain in the x coordinate
y0=-1#Start of the space domain in the y coordinate
yf=1#End of the space domain in the y coordinate
t0=0#Start time (Units = s)
tf=1.5#End time (Units = s) (approximately, due to variable CFL condition)
dPlot=5#For the plots (specifically, the number of iterations made between plots)
sigma_rho=1.5#Control parameter of the density Kreiss-Oliger dissipation filter term (Qd_rho) (sigma_rho >= 0 necessary!)
sigma_v_x=1.5#Control parameter of the x velocity component Kreiss-Oliger dissipation filter term (Qd_v_x) (sigma_v_x >= 0 necessary!)
sigma_v_y=1.5#Control parameter of the y velocity component Kreiss-Oliger dissipation filter term (Qd_v_y) (sigma_v_y >= 0 necessary!)
sigma_v_z=1.5#Control parameter of the z velocity component Kreiss-Oliger dissipation filter term (Qd_v_z) (sigma_v_z >= 0 necessary!)
sigma_p=1.5#Control parameter of the pressure Kreiss-Oliger dissipation filter term (Qd_p) (sigma_p >= 0 necessary!)
gamma=5/3#Constant adiabatic index
Amp_v=10**(-3)#Amplitude of the initial velocity gaussian perturbation (for linear case, Amp_v<<1 such as Amp_v = 0.1)
Amp_p=10**(-3)#Amplitude of the initial pressure gaussian perturbation (for linear case, Amp_p<<1 such as Amp_p = 0.1)
J=int((xf-x0)/dx)#Number of spatial points in the x coordinate
K=int((yf-y0)/dy)#Number of spatial points in the y coordinate
x_array=np.linspace(x0,xf,J)#Values of x
y_array=np.linspace(y0,yf,K)#Values of y
x_g_array,y_g_array=np.meshgrid(x_array,y_array)#Values of the x-y grid
Qd_rho_array=np.zeros((J,K))#Kreiss-Oliger dissipation filter term array for the continuity equation
Qd_v_x_array=np.zeros((J,K))#Kreiss-Oliger dissipation filter term array for the momentum equation in the x coordinate
Qd_v_y_array=np.zeros((J,K))#Kreiss-Oliger dissipation filter term array for the momentum equation in the y coordinate
Qd_v_z_array=np.zeros((J,K))#Kreiss-Oliger dissipation filter term array for the momentum equation in the z coordinate
Qd_p_array=np.zeros((J,K))#Kreiss-Oliger dissipation filter term array for the pressure equation
RHSContEq_array=np.zeros((J,K))#Right Hand Side array of the continuity equation
RHSMomEq_x_array=np.zeros((J,K))#Right Hand Side array of the x component of the momentum equation
RHSMomEq_y_array=np.zeros((J,K))#Right Hand Side array of the y component of the momentum equation
RHSMomEq_z_array=np.zeros((J,K))#Right Hand Side array of the z component of the momentum equation
RHSPresEq_array=np.zeros((J,K))#Right Hand Side array of the pressure equation
animation_frames_v_x_list=[]#Animation frames/artists for the x component of the velocity results
animation_frames_v_y_list=[]#Animation frames/artists for the y component of the velocity results
animation_frames_v_z_list=[]#Animation frames/artists for the z component of the velocity results
animation_frames_rho_list=[]#Animation frames/artists for the density results
animation_frames_p_list=[]#Animation frames/artists for the pressure results
animation_frames_E_list=[]#Animation frames/artists for the energy results
Num_frame=0#Number of solution frame for the animations
#Initial conditions for all the variables present in the equations (and general expression for v_s)
def ICv(x,y):#Initial condition for the velocity perturbation, in this case, a Gaussian function
    return Amp_v*np.exp(-(1/2)*((x-0)/(0.1))**2-(1/2)*((y-0)/(0.1))**2)
def ICv_zero(x,y):#Initial condition for the "equal to zero" velocity perturbation components
    return 0
def ICp(x,y):#Initial condition for the example of non perturbed pressure, in this case, unitary
    return 1
def ICp_pert(x,y):#Initial condition for the example with pressure perturbation, in this case, a Gaussian function
    return Amp_p*np.exp(-(1/2)*((x-0)/(0.1))**2-(1/2)*((y-0)/(0.1))**2)
def ICrho(x,y):#Initial condition for the density, in this case, unitary
    return 1
def ICE(p,rho,v_x,v_y,v_z):#Initial condition for the energy
    return ((p)/(gamma-1))+0.5*rho*(v_x**2+v_y**2+v_z**2)
def v_s(p,rho):#Expression of the sound velocity
    return np.sqrt((p*gamma)/(rho))
def InitialConditionVelocity_x(x_g_array,y_g_array):#Preparation for the implementation of the x component of the initial velocity condition
    v_x_array=np.zeros((J,K))
    for j in range(0,J):
        for k in range(0,K):
            v_x_array[j,k]=ICv(x_g_array[j,k],y_g_array[j,k])#Change here the initial condition to the desired case
    return v_x_array
def InitialConditionVelocity_y(x_g_array,y_g_array):#Preparation for the implementation of the y component of the initial velocity condition
    v_y_array=np.zeros((J,K))
    for j in range(0,J):
        for k in range(0,K):
            v_y_array[j,k]=ICv_zero(x_g_array[j,k],y_g_array[j,k])#Change here the initial condition to the desired case
    return v_y_array
def InitialConditionVelocity_z(x_g_array,y_g_array):#Preparation for the implementation of the z component of the initial velocity condition
    v_z_array=np.zeros((J,K))
    for j in range(0,J):
        for k in range(0,K):
            v_z_array[j,k]=ICv_zero(x_g_array[j,k],y_g_array[j,k])#Change here the initial condition to the desired case
    return v_z_array
def InitialConditionPressure(x_g_array,y_g_array):#Preparation for the implementation of the initial pressure condition
    p_array=np.zeros((J,K))
    for j in range(0,J):
        for k in range(0,K):
            p_array[j,k]=ICp(x_g_array[j,k],y_g_array[j,k])#Change here the initial condition to the desired case
    return p_array
def InitialConditionDensity(x_g_array,y_g_array):#Preparation for the implementation of the initial density condition
    rho_array=np.zeros((J,K))
    for j in range(0,J):
        for k in range(0,K):
            rho_array[j,k]=ICrho(x_g_array[j,k],y_g_array[j,k])
    return rho_array
def InitialConditionEnergy(p_array,rho_array,v_x_array,v_y_array,v_z_array):#Preparation for the implementation of the initial energy condition
    E_array=np.zeros((J,K))
    for j in range(0,J):
        for k in range(0,K):
            E_array[j,k]=ICE(p_array[j,k],rho_array[j,k],v_x_array[j,k],v_y_array[j,k],v_z_array[j,k])
    return E_array
def InitialConditionSoundVelocity(p_array,rho_array):#Preparation for the implementation of any component of the initial sound velocity condition (because of the isotropic property)
    v_s_array=np.zeros((J,K))
    for j in range(0,J):
        for k in range(0,K):
            v_s_array[j,k]=v_s(p_array[j,k],rho_array[j,k])
    return v_s_array
#Numerical methods utilized for each equation (general and specific ones)
def Finite_Difference_Derivative_Order4_x(u_array,j,k):#General spatial discretization for the x coordinate: Finite/discrete central differences (4rth order)
    dudx=(1/(12*dx))*(-u_array[j+2,k]+8*u_array[j+1,k]-8*u_array[j-1,k]+u_array[j-2,k])
    return dudx
def Finite_Difference_Derivative_Order4_y(u_array,j,k):#General spatial discretization for the y coordinate: Finite/discrete central differences (4rth order)
    dudy=(1/(12*dy))*(-u_array[j,k+2]+8*u_array[j,k+1]-8*u_array[j,k-1]+u_array[j,k-2])
    return dudy
def Periodic_BC(u_array):#General (circular) periodic boundary conditions
    u_array[0,:]=u_array[-4,:]
    u_array[1,:]=u_array[-3,:]
    u_array[-2,:]=u_array[2,:]
    u_array[-1,:]=u_array[3,:]
    u_array[:,0]=u_array[:,-4]
    u_array[:,1]=u_array[:,-3]
    u_array[:,-2]=u_array[:,2]
    u_array[:,-1]=u_array[:,3]
    return u_array
def RHSContEq(v_x_array,v_y_array,rho_array):#Right Hand Side of the continuity equation
    for j in range(2,J-2):
        for k in range(2,K-2):
            if j==2 or k==2 or j==J-3 or k==K-3:
                Qd_rho_array[j,k]=0
            else:
                Qd_rho_array[j,k]=(1/(64*dx))*(rho_array[j+3,k]-6*rho_array[j+2,k]+15*rho_array[j+1,k]-20*rho_array[j,k]+15*rho_array[j-1,k]-6*rho_array[j-2,k]+rho_array[j-3,k])+(1/(64*dy))*(rho_array[j,k+3]-6*rho_array[j,k+2]+15*rho_array[j,k+1]-20*rho_array[j,k]+15*rho_array[j,k-1]-6*rho_array[j,k-2]+rho_array[j,k-3])
            RHSContEq_array[j,k]=-rho_array[j,k]*(Finite_Difference_Derivative_Order4_x(v_x_array,j,k)+Finite_Difference_Derivative_Order4_y(v_y_array,j,k))-v_x_array[j,k]*Finite_Difference_Derivative_Order4_x(rho_array,j,k)-v_y_array[j,k]*Finite_Difference_Derivative_Order4_y(rho_array,j,k)+sigma_rho*Qd_rho_array[j,k]
    return RHSContEq_array
def RHSMomEq_x(v_x_array,v_y_array,rho_array,p_array):#Right Hand Side of the x component of the momentum equation
    for j in range(2,J-2):
        for k in range(2,K-2):
            if j==2 or k==2 or j==J-3 or k==K-3:
                Qd_v_x_array[j,k]=0
            else:
                Qd_v_x_array[j,k]=(1/(64*dx))*(v_x_array[j+3,k]-6*v_x_array[j+2,k]+15*v_x_array[j+1,k]-20*v_x_array[j,k]+15*v_x_array[j-1,k]-6*v_x_array[j-2,k]+v_x_array[j-3,k])+(1/(64*dy))*(v_x_array[j,k+3]-6*v_x_array[j,k+2]+15*v_x_array[j,k+1]-20*v_x_array[j,k]+15*v_x_array[j,k-1]-6*v_x_array[j,k-2]+v_x_array[j,k-3])
            RHSMomEq_x_array[j,k]=-((1)/(rho_array[j,k]))*Finite_Difference_Derivative_Order4_x(p_array,j,k)-v_x_array[j,k]*Finite_Difference_Derivative_Order4_x(v_x_array,j,k)-v_y_array[j,k]*Finite_Difference_Derivative_Order4_y(v_x_array,j,k)+sigma_v_x*Qd_v_x_array[j,k]
    return RHSMomEq_x_array
def RHSMomEq_y(v_x_array,v_y_array,rho_array,p_array):#Right Hand Side of the y component of the momentum equation
    for j in range(2,J-2):
        for k in range(2,K-2):
            if j==2 or k==2 or j==J-3 or k==K-3:
                Qd_v_y_array[j,k]=0
            else:
                Qd_v_y_array[j,k]=(1/(64*dx))*(v_y_array[j+3,k]-6*v_y_array[j+2,k]+15*v_y_array[j+1,k]-20*v_y_array[j,k]+15*v_y_array[j-1,k]-6*v_y_array[j-2,k]+v_y_array[j-3,k])+(1/(64*dy))*(v_y_array[j,k+3]-6*v_y_array[j,k+2]+15*v_y_array[j,k+1]-20*v_y_array[j,k]+15*v_y_array[j,k-1]-6*v_y_array[j,k-2]+v_y_array[j,k-3])
            RHSMomEq_y_array[j,k]=-((1)/(rho_array[j,k]))*Finite_Difference_Derivative_Order4_y(p_array,j,k)-v_x_array[j,k]*Finite_Difference_Derivative_Order4_x(v_y_array,j,k)-v_y_array[j,k]*Finite_Difference_Derivative_Order4_y(v_y_array,j,k)+sigma_v_y*Qd_v_y_array[j,k]
    return RHSMomEq_y_array
def RHSMomEq_z(v_x_array,v_y_array,v_z_array):#Right Hand Side of the z component of the momentum equation
    for j in range(2,J-2):
        for k in range(2,K-2):
            if j==2 or k==2 or j==J-3 or k==K-3:
                Qd_v_z_array[j,k]=0
            else:
                Qd_v_z_array[j,k]=(1/(64*dx))*(v_z_array[j+3,k]-6*v_z_array[j+2,k]+15*v_z_array[j+1,k]-20*v_z_array[j,k]+15*v_z_array[j-1,k]-6*v_z_array[j-2,k]+v_z_array[j-3,k])+(1/(64*dy))*(v_z_array[j,k+3]-6*v_z_array[j,k+2]+15*v_z_array[j,k+1]-20*v_z_array[j,k]+15*v_z_array[j,k-1]-6*v_z_array[j,k-2]+v_z_array[j,k-3])
            RHSMomEq_z_array[j,k]=-v_x_array[j,k]*Finite_Difference_Derivative_Order4_x(v_z_array,j,k)-v_y_array[j,k]*Finite_Difference_Derivative_Order4_y(v_z_array,j,k)+sigma_v_z*Qd_v_z_array[j,k]
    return RHSMomEq_z_array
def RHSPresEq(v_x_array,v_y_array,p_array):#Right Hand Side of the pressure equation
    for j in range(2,J-2):
        for k in range(2,K-2):
            if j==2 or k==2 or j==J-3 or k==K-3:
                Qd_p_array[j,k]=0
            else:
                Qd_p_array[j,k]=(1/(64*dx))*(p_array[j+3,k]-6*p_array[j+2,k]+15*p_array[j+1,k]-20*p_array[j,k]+15*p_array[j-1,k]-6*p_array[j-2,k]+p_array[j-3,k])+(1/(64*dy))*(p_array[j,k+3]-6*p_array[j,k+2]+15*p_array[j,k+1]-20*p_array[j,k]+15*p_array[j,k-1]-6*p_array[j,k-2]+p_array[j,k-3])
            RHSPresEq_array[j,k]=-gamma*p_array[j,k]*(Finite_Difference_Derivative_Order4_x(v_x_array,j,k)+Finite_Difference_Derivative_Order4_y(v_y_array,j,k))-v_x_array[j,k]*Finite_Difference_Derivative_Order4_x(p_array,j,k)-v_y_array[j,k]*Finite_Difference_Derivative_Order4_y(p_array,j,k)+sigma_p*Qd_p_array[j,k]
    return RHSPresEq_array
def RK4(v_x_array,v_y_array,v_z_array,rho_array,p_array,E_array):#Temporal discretization: Runge-Kutta of fourth order (RK4)
    k_1_ContEq=dt*RHSContEq(v_x_array,v_y_array,rho_array)
    k_1_MomEq_x=dt*RHSMomEq_x(v_x_array,v_y_array,rho_array,p_array)
    k_1_MomEq_y=dt*RHSMomEq_y(v_x_array,v_y_array,rho_array,p_array)
    k_1_MomEq_z=dt*RHSMomEq_z(v_x_array,v_y_array,v_z_array)
    k_1_PresEq=dt*RHSPresEq(v_x_array,v_y_array,p_array)
    rho_i_array=rho_array+0.5*k_1_ContEq
    v_x_i_array=v_x_array+0.5*k_1_MomEq_x
    v_y_i_array=v_y_array+0.5*k_1_MomEq_y
    v_z_i_array=v_z_array+0.5*k_1_MomEq_z
    p_i_array=p_array+0.5*k_1_PresEq
    rho_i_array=Periodic_BC(rho_i_array)#Imposing periodic BCs
    v_x_i_array=Periodic_BC(v_x_i_array)#Imposing periodic BCs
    v_y_i_array=Periodic_BC(v_y_i_array)#Imposing periodic BCs
    v_z_i_array=Periodic_BC(v_z_i_array)#Imposing periodic BCs
    p_i_array=Periodic_BC(p_i_array)#Imposing periodic BCs
    k_2_ContEq=dt*RHSContEq(v_x_i_array,v_y_i_array,rho_i_array)
    k_2_MomEq_x=dt*RHSMomEq_x(v_x_i_array,v_y_i_array,rho_i_array,p_i_array)
    k_2_MomEq_y=dt*RHSMomEq_y(v_x_i_array,v_y_i_array,rho_i_array,p_i_array)
    k_2_MomEq_z=dt*RHSMomEq_z(v_x_i_array,v_y_i_array,v_z_i_array)
    k_2_PresEq=dt*RHSPresEq(v_x_i_array,v_y_i_array,p_i_array)
    rho_i_array=rho_array+0.5*k_2_ContEq
    v_x_i_array=v_x_array+0.5*k_2_MomEq_x
    v_y_i_array=v_y_array+0.5*k_2_MomEq_y
    v_z_i_array=v_z_array+0.5*k_2_MomEq_z
    p_i_array=p_array+0.5*k_2_PresEq
    rho_i_array=Periodic_BC(rho_i_array)#Imposing periodic BCs
    v_x_i_array=Periodic_BC(v_x_i_array)#Imposing periodic BCs
    v_y_i_array=Periodic_BC(v_y_i_array)#Imposing periodic BCs
    v_z_i_array=Periodic_BC(v_z_i_array)#Imposing periodic BCs
    p_i_array=Periodic_BC(p_i_array)#Imposing periodic BCs
    k_3_ContEq=dt*RHSContEq(v_x_i_array,v_y_i_array,rho_i_array)
    k_3_MomEq_x=dt*RHSMomEq_x(v_x_i_array,v_y_i_array,rho_i_array,p_i_array)
    k_3_MomEq_y=dt*RHSMomEq_y(v_x_i_array,v_y_i_array,rho_i_array,p_i_array)
    k_3_MomEq_z=dt*RHSMomEq_z(v_x_i_array,v_y_i_array,v_z_i_array)
    k_3_PresEq=dt*RHSPresEq(v_x_i_array,v_y_i_array,p_i_array)
    rho_i_array=rho_array+k_3_ContEq
    v_x_i_array=v_x_array+k_3_MomEq_x
    v_y_i_array=v_y_array+k_3_MomEq_y
    v_z_i_array=v_z_array+k_3_MomEq_z
    p_i_array=p_array+k_3_PresEq
    rho_i_array=Periodic_BC(rho_i_array)#Imposing periodic BCs
    v_x_i_array=Periodic_BC(v_x_i_array)#Imposing periodic BCs
    v_y_i_array=Periodic_BC(v_y_i_array)#Imposing periodic BCs
    v_z_i_array=Periodic_BC(v_z_i_array)#Imposing periodic BCs
    p_i_array=Periodic_BC(p_i_array)#Imposing periodic BCs
    k_4_ContEq=dt*RHSContEq(v_x_i_array,v_y_i_array,rho_i_array)
    k_4_MomEq_x=dt*RHSMomEq_x(v_x_i_array,v_y_i_array,rho_i_array,p_i_array)
    k_4_MomEq_y=dt*RHSMomEq_y(v_x_i_array,v_y_i_array,rho_i_array,p_i_array)
    k_4_MomEq_z=dt*RHSMomEq_z(v_x_i_array,v_y_i_array,v_z_i_array)
    k_4_PresEq=dt*RHSPresEq(v_x_i_array,v_y_i_array,p_i_array)
    drhodt_array=rho_array+((1)/(6))*(k_1_ContEq+2*k_2_ContEq+2*k_3_ContEq+k_4_ContEq)
    dv_xdt_array=v_x_array+((1)/(6))*(k_1_MomEq_x+2*k_2_MomEq_x+2*k_3_MomEq_x+k_4_MomEq_x)
    dv_ydt_array=v_y_array+((1)/(6))*(k_1_MomEq_y+2*k_2_MomEq_y+2*k_3_MomEq_y+k_4_MomEq_y)
    dv_zdt_array=v_z_array+((1)/(6))*(k_1_MomEq_z+2*k_2_MomEq_z+2*k_3_MomEq_z+k_4_MomEq_z)
    dpdt_array=p_array+((1)/(6))*(k_1_PresEq+2*k_2_PresEq+2*k_3_PresEq+k_4_PresEq)
    return [drhodt_array,dv_xdt_array,dv_ydt_array,dv_zdt_array,dpdt_array]
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Solution evolution and plotting of the results
v_x_array=InitialConditionVelocity_x(x_g_array,y_g_array)
v_y_array=InitialConditionVelocity_y(x_g_array,y_g_array)
v_z_array=InitialConditionVelocity_z(x_g_array,y_g_array)
p_array=InitialConditionPressure(x_g_array,y_g_array)
rho_array=InitialConditionDensity(x_g_array,y_g_array)
E_array=InitialConditionEnergy(p_array,rho_array,v_x_array,v_y_array,v_z_array)
v_s_array=InitialConditionSoundVelocity(p_array,rho_array)
E_new_array=E_array
v_s_new_array=v_s_array
dt=0.5*(min(dx,dy)/v_s_array.max())#Variable time-step size (Units = s). For stability purposes, we use the CFL criteria (lambda=((max(v_s)*dt)/min(dx,dy))<=1) to obtain this value.
N=int((tf-t0)/dt)#Time discretitzation
#Preparation of the figures' environments
#v_x(x,y,t) vs x,y
fig_v_x,ax_v_x=plt.subplots(figsize=(9,5))
ax_v_x.set_xlim(x0,xf)
ax_v_x.set_ylim(y0,yf)
ax_v_x.set_title(r"Solution of $v_{x}$ in the 2D Euler equations")
ax_v_x.set_xlabel(r"$x$")
ax_v_x.set_ylabel(r"$y$")
ax_v_x.grid(True)
fig_v_x.set_tight_layout(True)
#v_y(x,y,t) vs x,y
fig_v_y,ax_v_y=plt.subplots(figsize=(9,5))
ax_v_y.set_xlim(x0,xf)
ax_v_y.set_ylim(y0,yf)
ax_v_y.set_title(r"Solution of $v_{y}$ in the 2D Euler equations")
ax_v_y.set_xlabel(r"$x$")
ax_v_y.set_ylabel(r"$y$")
ax_v_y.grid(True)
fig_v_y.set_tight_layout(True)
#v_z(x,y,t) vs x,y
fig_v_z,ax_v_z=plt.subplots(figsize=(9,5))
ax_v_z.set_xlim(x0,xf)
ax_v_z.set_ylim(y0,yf)
ax_v_z.set_title(r"Solution of $v_{z}$ in the 2D Euler equations")
ax_v_z.set_xlabel(r"$x$")
ax_v_z.set_ylabel(r"$y$")
ax_v_z.grid(True)
fig_v_z.set_tight_layout(True)
#\rho(x,y,t) vs x,y
fig_rho,ax_rho=plt.subplots(figsize=(9,5))
ax_rho.set_xlim(x0,xf)
ax_rho.set_ylim(y0,yf)
ax_rho.set_title(r"Solution of $\rho$ in the 2D Euler equations")
ax_rho.set_xlabel(r"$x$")
ax_rho.set_ylabel(r"$y$")
ax_rho.grid(True)
fig_rho.set_tight_layout(True)
#p(x,y,t) vs x,y
fig_p,ax_p=plt.subplots(figsize=(9,5))
ax_p.set_xlim(x0,xf)
ax_p.set_ylim(y0,yf)
ax_p.set_title(r"Solution of $p$ in the 2D Euler equations")
ax_p.set_xlabel(r"$x$")
ax_p.set_ylabel(r"$y$")
ax_p.grid(True)
fig_p.set_tight_layout(True)
#E(x,y,t) vs x,y
fig_E,ax_E=plt.subplots(figsize=(9,5))
ax_E.set_xlim(x0,xf)
ax_E.set_ylim(y0,yf)
ax_E.set_title(r"Solution of $E$ in the 2D Euler equations")
ax_E.set_xlabel(r"$x$")
ax_E.set_ylabel(r"$y$")
ax_E.grid()
fig_E.set_tight_layout(True)
for n in range(N+1):
    if n % dPlot==0:
        #v_x(x,y,t) vs x,y
        frame=ax_v_x.contourf(x_g_array,y_g_array,v_x_array,cmap=cm.seismic,antialiased=False)
        if Num_frame==0:
            cb_v_x=fig_v_x.colorbar(frame,shrink=0.5,aspect=5,ax=ax_v_x,pad=0.05,location="right")
        if Num_frame==1:
            cb_v_x.remove()
            cb_v_x=fig_v_x.colorbar(frame,shrink=0.5,aspect=5,ax=ax_v_x,pad=0.05,location="right")
        animation_frames_v_x_list.append([frame])
        #v_y(x,y,t) vs x,y
        frame=ax_v_y.contourf(x_g_array,y_g_array,v_y_array,cmap=cm.seismic,antialiased=False)
        if Num_frame==0:
            cb_v_y=fig_v_y.colorbar(frame,shrink=0.5,aspect=5,ax=ax_v_y,pad=0.05,location="right")
        if Num_frame==1:
            cb_v_y.remove()
            cb_v_y=fig_v_y.colorbar(frame,shrink=0.5,aspect=5,ax=ax_v_y,pad=0.05,location="right")
        animation_frames_v_y_list.append([frame])
        #v_z(x,y,t) vs x,y
        frame=ax_v_z.contourf(x_g_array,y_g_array,v_z_array,cmap=cm.seismic,antialiased=False)
        if Num_frame==0:
            cb_v_z=fig_v_z.colorbar(frame,shrink=0.5,aspect=5,ax=ax_v_z,pad=0.05,location="right")
        if Num_frame==1:
            cb_v_z.remove()
            cb_v_z=fig_v_z.colorbar(frame,shrink=0.5,aspect=5,ax=ax_v_z,pad=0.05,location="right")
        animation_frames_v_z_list.append([frame])
        #\rho(x,y,t) vs x,y
        frame=ax_rho.contourf(x_g_array,y_g_array,rho_array,cmap=cm.seismic,antialiased=False)
        if Num_frame==0:
            cb_rho=fig_rho.colorbar(frame,shrink=0.5,aspect=5,ax=ax_rho,pad=0.05,location="right")
        if Num_frame==1:
            cb_rho.remove()
            cb_rho=fig_rho.colorbar(frame,shrink=0.5,aspect=5,ax=ax_rho,pad=0.05,location="right")
        animation_frames_rho_list.append([frame])
        #p(x,y,t) vs x,y
        frame=ax_p.contourf(x_g_array,y_g_array,p_array,cmap=cm.seismic,antialiased=False)
        if Num_frame==0:
            cb_p=fig_p.colorbar(frame,shrink=0.5,aspect=5,ax=ax_p,pad=0.05,location="right")
        if Num_frame==1:
            cb_p.remove()
            cb_p=fig_p.colorbar(frame,shrink=0.5,aspect=5,ax=ax_p,pad=0.05,location="right")
        animation_frames_p_list.append([frame])
        #E(x,y,t) vs x,y
        frame=ax_E.contourf(x_g_array,y_g_array,E_array,cmap=cm.seismic,antialiased=False)
        if Num_frame==0:
            cb_E=fig_E.colorbar(frame,shrink=0.5,aspect=5,ax=ax_E,pad=0.05,location="right")
        if Num_frame==1:
            cb_E.remove()
            cb_E=fig_E.colorbar(frame,shrink=0.5,aspect=5,ax=ax_E,pad=0.05,location="right")
        animation_frames_E_list.append([frame])
        print("Frame ",Num_frame," computed at t = ",str(n*dt)," s.",sep="")
        Num_frame+=1
    RK4_list=RK4(v_x_array,v_y_array,v_z_array,rho_array,p_array,E_array)#Time evolution for the differential equations
    rho_new_array=RK4_list[0]
    v_x_new_array=RK4_list[1]
    v_y_new_array=RK4_list[2]
    v_z_new_array=RK4_list[3]
    p_new_array=RK4_list[4]
    rho_new_array=Periodic_BC(rho_new_array)#Imposing BCs
    v_x_new_array=Periodic_BC(v_x_new_array)#Imposing BCs
    v_y_new_array=Periodic_BC(v_y_new_array)#Imposing BCs
    v_z_new_array=Periodic_BC(v_z_new_array)#Imposing BCs
    p_new_array=Periodic_BC(p_new_array)#Imposing BCs
    E_new_array=((p_new_array)/(gamma-1))+0.5*rho_new_array*(v_x_new_array**2+v_y_new_array**2+v_z_new_array**2)#"Time" evolution of E
    v_s_new_array=v_s(p_new_array,rho_new_array)#"Time" evolution of v_s
    rho_array=rho_new_array#Update for the next time step
    v_x_array=v_x_new_array#Update for the next time step
    v_y_array=v_y_new_array#Update for the next time step
    v_z_array=v_z_new_array#Update for the next time step
    p_array=p_new_array#Update for the next time step
    E_array=E_new_array#Update for the next time step
    v_s_array=v_s_new_array#Update for the next time step
    dt=0.5*(min(dx,dy)/v_s_array.max())#Variable time-step size (Units = s). For stability purposes, we use the CFL criteria (lambda=((max(v_s)*dt)/min(dx,dy))<=1) to obtain this value.
    N=int((tf-t0)/dt)#Time discretitzation
print("All frames computed, creating animations...")
animation_v_x=animation.ArtistAnimation(fig=fig_v_x,artists=animation_frames_v_x_list,interval=200)#Animation of the x component of the velocity results
animation_v_x.save(filename="2D Euler Solution for v_x.mp4",writer="ffmpeg")#Animation save of the x component of the velocity results
animation_v_y=animation.ArtistAnimation(fig=fig_v_y,artists=animation_frames_v_y_list,interval=200)#Animation of the y component of the velocity results
animation_v_y.save(filename="2D Euler Solution for v_y.mp4",writer="ffmpeg")#Animation save of the y component of the velocity results
animation_v_z=animation.ArtistAnimation(fig=fig_v_z,artists=animation_frames_v_z_list,interval=200)#Animation of the z component of the velocity results
animation_v_z.save(filename="2D Euler Solution for v_z.mp4",writer="ffmpeg")#Animation save of the z component of the velocity results
animation_rho=animation.ArtistAnimation(fig=fig_rho,artists=animation_frames_rho_list,interval=200)#Animation of the density results
animation_rho.save(filename="2D Euler Solution for rho.mp4",writer="ffmpeg")#Animation save of the density results
animation_p=animation.ArtistAnimation(fig=fig_p,artists=animation_frames_p_list,interval=200)#Animation of the pressure results
animation_p.save(filename="2D Euler Solution for p.mp4",writer="ffmpeg")#Animation save of the pressure results
animation_E=animation.ArtistAnimation(fig=fig_E,artists=animation_frames_E_list,interval=200)#Animation of the energy results
animation_E.save(filename="2D Euler Solution for E.mp4",writer="ffmpeg")#Animation save of the energy results
plt.close()
print("Animations saved at ",os.getcwd(),".",sep="")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")