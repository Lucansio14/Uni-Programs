# -*- coding: utf-8 -*-
"""
Created on Thu May 15 11:29:45 2025
Solver of the 1D MHD equations: A "1D" (initial perturbations only in one direction (x), which means derivatives in other directions (y and z) equal to zero or no spatial variations in other directions) non linear equation system in the ideal and totally ionized plasma form for pressure (p), density (\rho), velocity (v), and magnetic field (B), used primarly in magnetohydrodynamics, with RK4 (time part), FD4 (spatial part) and Kreiss-Oliger dissipation filter.
The corresponding general equations in cartesian coordinates (considering adimensional variables and constants, thermodynamic equilibrium, Ohm's law for B, and neglecting gravitational effects due to relatively small spatial domain):
    \partial_{t}\rho = -\div{\vec{v}\rho} (Continuity equation)
    \partial_{t}(\vec{v}) = -(1/\rho)\grad{p} + (1/\mu\rho)\Rot{\vec{B}} x \vec{B} - (\vec{v} · \nabla)\vec{v} (Momentum equation, because momentum (m) = density (\rho) * velocity (v) (this will be mayorly used in the computations))
    \partial_{t}(\vec{B}) = \Rot{\vec{v} x \vec{B}} (Induction equation)
    \partial_{t}p = -p\gamma\div{\vec{v}} - (\vec{v} · \nabla)p ((Adiabatic) pressure equation)
    \div{B} = 0 (Solenoidal initial condition (because of the induction equation and constant initial magnetic field))
    (No state equation needed due to considering temperature as uniform, constant and not a variable in the system)
(where \gamma is the adiabatic index and, importantly for later, v_s = \sqrt{(p\gamma)/(\rho)} is the adiabatic sound velocity and v_A = \sqrt{(B^2)/(\mu\rho)} is the Alfvén velocity)
More information about these equations in these websites: https://warwick.ac.uk/fac/sci/physics/research/cfsa/people/valery/teaching/khu_mhd/KHU_mhd_handout.pdf
                                                          https://lweb.cfa.harvard.edu/~namurphy/Lectures/Ay253_01_IdealMHD.pdf
                                                          https://perswww.kuleuven.be/~u0016541/MHD_sheets_pdf/nsap430m.06.4.pdf
More information about CFL criteria for non linear equations in this website: https://scicomp.stackexchange.com/questions/31456/cfl-equation-for-non-linear-equation

Context about the computed solution (with initial perturbations in the velocity):
Without the use of more sophisticated methods (like shock capturing or hyperdifusion), the non linear case (not Amp_v<<1, such as Amp_v = 1) solutions are "corrupted" in some way by artifact spikes. The reason of this issue is completely numerical. For that, only the linear case (Amp_v<<1, such as Amp_v = 0.1) solution is correctly obtainable with this code.
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
#Definition of general constants, variables, lists, arrays and functions
dpos=0.02#Space-step size (unique, given the 1D dimensions and symmetric grid)
pos0=-1#Start of the space domain (unique, given the 1D dimensions and symmetric grid)
posf=1#End of the space domain (unique, given the 1D dimensions and symmetric grid)
t0=0#Start time (Units = s)
tf=1.5#End time (Units = s) (approximately, due to variable CFL condition)
dPlot=20#For the plots (specifically, the number of iterations made between plots)
sigma_rho=1.5#Control parameter of the density Kreiss-Oliger dissipation filter term (Qd_rho) (sigma_rho >= 0 necessary!)
sigma_v_x=1.5#Control parameter of the x velocity component Kreiss-Oliger dissipation filter term (Qd_v_x) (sigma_v_x >= 0 necessary!)
sigma_v_y=1.5#Control parameter of the y velocity component Kreiss-Oliger dissipation filter term (Qd_v_y) (sigma_v_y >= 0 necessary!)
sigma_v_z=1.5#Control parameter of the z velocity component Kreiss-Oliger dissipation filter term (Qd_v_z) (sigma_v_z >= 0 necessary!)
sigma_B_x=1.5#Control parameter of the x magnetic field component Kreiss-Oliger dissipation filter term (Qd_B_x) (sigma_B_x >= 0 necessary!)
sigma_B_y=1.5#Control parameter of the y magnetic field component Kreiss-Oliger dissipation filter term (Qd_B_y) (sigma_B_y >= 0 necessary!)
sigma_B_z=1.5#Control parameter of the z magnetic field component Kreiss-Oliger dissipation filter term (Qd_B_z) (sigma_B_z >= 0 necessary!)
sigma_p=1.5#Control parameter of the pressure Kreiss-Oliger dissipation filter term (Qd_p) (sigma_E >= 0 necessary!)
gamma=5/3#Constant adiabatic index
mu=1#Permeability constant (in this case, adimensional)
Amp_v=10**(-3)#Amplitude of the initial velocity gaussian perturbation (for linear case, Amp_v<<1 such as Amp_v = 0.1)
J=int((posf-pos0)/dpos)#Number of spatial points
pos_array=np.linspace(pos0,posf,J)#Values of x (or y or z, due to considering a symmetric grid in all directions in respect to the perturbation direction x)
Qd_rho_array=np.zeros(J)#Kreiss-Oliger dissipation filter term array for the continuity equation
Qd_v_x_array=np.zeros(J)#Kreiss-Oliger dissipation filter term array for the momentum equation in the x coordinate
Qd_v_y_array=np.zeros(J)#Kreiss-Oliger dissipation filter term array for the momentum equation in the y coordinate
Qd_v_z_array=np.zeros(J)#Kreiss-Oliger dissipation filter term array for the momentum equation in the z coordinate
Qd_B_x_array=np.zeros(J)#Kreiss-Oliger dissipation filter term array for the induction equation in the x coordinate
Qd_B_y_array=np.zeros(J)#Kreiss-Oliger dissipation filter term array for the induction equation in the y coordinate
Qd_B_z_array=np.zeros(J)#Kreiss-Oliger dissipation filter term array for the induction equation in the z coordinate
Qd_p_array=np.zeros(J)#Kreiss-Oliger dissipation filter term array for the pressure equation
RHSContEq_array=np.zeros(J)#Right Hand Side array of the continuity equation
RHSMomEq_x_array=np.zeros(J)#Right Hand Side array of the x component of the momentum equation
RHSMomEq_y_array=np.zeros(J)#Right Hand Side array of the y component of the momentum equation
RHSMomEq_z_array=np.zeros(J)#Right Hand Side array of the z component of the momentum equation
RHSIndEq_x_array=np.zeros(J)#Right Hand Side array of the x component of the induction equation
RHSIndEq_y_array=np.zeros(J)#Right Hand Side array of the y component of the induction equation
RHSIndEq_z_array=np.zeros(J)#Right Hand Side array of the z component of the induction equation
RHSPresEq_array=np.zeros(J)#Right Hand Side array of the pressure equation
animation_frames_v_x_list=[]#Animation frames/artists for the x component of the velocity results
animation_frames_v_y_list=[]#Animation frames/artists for the y component of the velocity results
animation_frames_v_z_list=[]#Animation frames/artists for the z component of the velocity results
animation_frames_B_x_list=[]#Animation frames/artists for the x component of the magnetic field results
animation_frames_B_y_list=[]#Animation frames/artists for the y component of the magnetic field results
animation_frames_B_z_list=[]#Animation frames/artists for the z component of the magnetic field results
animation_frames_rho_list=[]#Animation frames/artists for the density results
animation_frames_p_list=[]#Animation frames/artists for the pressure results
animation_frames_v_s_list=[]#Animation frames/artists for the sound velocity results (due to the isotropic property, all the components are equal)
animation_frames_v_A_x_list=[]#Animation frames/artists for the x component of the Alfvén velocity results
animation_frames_v_A_y_list=[]#Animation frames/artists for the y component of the Alfvén velocity results
animation_frames_v_A_z_list=[]#Animation frames/artists for the z component of the Alfvén velocity results
Num_frame=0#Number of solution frame for the animations
#Initial conditions for all the variables present in the equations (and general expression for v_s and v_A)
def ICv(x):#Initial condition for the velocity perturbation component, in this case, a Gaussian function
    return Amp_v*np.exp(-(1/2)*((x-0)/(0.1))**2)
def ICv_zero(x):#Initial condition for the "equal to zero" velocity components
    return 0
def ICp(x):#Initial condition for the pressure, in this case, unitary
    return 1
def ICrho(x):#Initial condition for the density, in this case, unitary
    return 1
def ICB(x):#Initial condition for the non-zero magnetic field components, in this case, unitary
    return 1
def ICB_zero(x):#Initial condition for the "equal to zero" magnetic field components
    return 0
def v_s(p,rho):#Expression of the sound velocity (1D because of isotropic property)
    return np.sqrt((p*gamma)/(rho))
def v_A(B,rho):#Expression of any component of the Alfvén velocity
    return np.sqrt((B**2)/(mu*rho))
def InitialConditionVelocity_x(pos_array):#Preparation for the implementation of the x component of the initial velocity condition
    j=0
    v_array=np.zeros(J)
    for pos in pos_array:
        v_array[j]=ICv(pos)#Change here the initial condition to the desired case
        j+=1
    return v_array
def InitialConditionVelocity_y(pos_array):#Preparation for the implementation of the y component of the initial velocity condition
    j=0
    v_array=np.zeros(J)
    for pos in pos_array:
        v_array[j]=ICv_zero(pos)#Change here the initial condition to the desired case
        j+=1
    return v_array
def InitialConditionVelocity_z(pos_array):#Preparation for the implementation of the z component of the initial velocity condition
    j=0
    v_array=np.zeros(J)
    for pos in pos_array:
        v_array[j]=ICv_zero(pos)#Change here the initial condition to the desired case
        j+=1
    return v_array
def InitialConditionPressure(pos_array):#Preparation for the implementation of the initial pressure condition
    j=0
    p_array=np.zeros(J)
    for pos in pos_array:
        p_array[j]=ICp(pos)
        j+=1
    return p_array
def InitialConditionDensity(pos_array):#Preparation for the implementation of the initial density condition
    j=0
    rho_array=np.zeros(J)
    for pos in pos_array:
        rho_array[j]=ICrho(pos)
        j+=1
    return rho_array
def InitialConditionMagneticField_x(pos_array):#Preparation for the implementation of the x component of the initial magnetic field condition
    j=0
    B_array=np.zeros(J)
    for pos in pos_array:
        B_array[j]=ICB_zero(pos)#Change here the initial condition to the desired case
        j+=1
    return B_array
def InitialConditionMagneticField_y(pos_array):#Preparation for the implementation of the y component of the initial magnetic field condition
    j=0
    B_array=np.zeros(J)
    for pos in pos_array:
        B_array[j]=ICB_zero(pos)#Change here the initial condition to the desired case
        j+=1
    return B_array
def InitialConditionMagneticField_z(pos_array):#Preparation for the implementation of the z component of the initial magnetic field condition
    j=0
    B_array=np.zeros(J)
    for pos in pos_array:
        B_array[j]=ICB(pos)#Change here the initial condition to the desired case
        j+=1
    return B_array
def InitialConditionSoundVelocity(p_array,rho_array):#Preparation for the implementation of any component of the initial sound velocity condition (because of the isotropic property)
    j=0
    v_s_array=np.zeros(J)
    while j<J:
        v_s_array[j]=v_s(p_array[j],rho_array[j])
        j+=1
    return v_s_array
def InitialConditionAlfvenVelocity(B_array,rho_array):#Preparation for the implementation of any component of the initial Alfvén velocity condition
    j=0
    v_A_array=np.zeros(J)
    while j<J:
        v_A_array[j]=v_A(B_array[j],rho_array[j])
        j+=1
    return v_A_array
#Numerical methods utilized for each equation (general and specific ones)
def Finite_Difference_Derivative_Order4(u_array,j):#General spatial discretization: Finite/discrete central differences (4rth order)
    dudr=(1/(12*dpos))*(-u_array[j+2]+8*u_array[j+1]-8*u_array[j-1]+u_array[j-2])
    return dudr
def Periodic_BC(u_array):#General (circular) periodic boundary conditions
    u_array[0]=u_array[-4]
    u_array[1]=u_array[-3]
    u_array[-2]=u_array[2]
    u_array[-1]=u_array[3]
    return u_array
def RHSContEq(v_x_array,rho_array):#Right Hand Side of the continuity equation
    for j in range(2,J-2):
        if j==2 or j==J-3:
            Qd_rho_array[j]=0
        else:
            Qd_rho_array[j]=(1/(64*dpos))*(rho_array[j+3]-6*rho_array[j+2]+15*rho_array[j+1]-20*rho_array[j]+15*rho_array[j-1]-6*rho_array[j-2]+rho_array[j-3])
        RHSContEq_array[j]=-rho_array[j]*Finite_Difference_Derivative_Order4(v_x_array,j)-v_x_array[j]*Finite_Difference_Derivative_Order4(rho_array,j)+sigma_rho*Qd_rho_array[j]
    return RHSContEq_array
def RHSMomEq_x(v_x_array,rho_array,p_array,B_y_array,B_z_array):#Right Hand Side of the x component of the momentum equation
    for j in range(2,J-2):
        if j==2 or j==J-3:
            Qd_v_x_array[j]=0
        else:
            Qd_v_x_array[j]=(1/(64*dpos))*(v_x_array[j+3]-6*v_x_array[j+2]+15*v_x_array[j+1]-20*v_x_array[j]+15*v_x_array[j-1]-6*v_x_array[j-2]+v_x_array[j-3])
        RHSMomEq_x_array[j]=-((1)/(rho_array[j]))*Finite_Difference_Derivative_Order4(p_array,j)+((1)/(mu*rho_array[j]))*(-B_z_array[j]*Finite_Difference_Derivative_Order4(B_z_array,j)-B_y_array[j]*Finite_Difference_Derivative_Order4(B_y_array,j))-v_x_array[j]*Finite_Difference_Derivative_Order4(v_x_array,j)+sigma_v_x*Qd_v_x_array[j]
    return RHSMomEq_x_array
def RHSMomEq_y(v_x_array,v_y_array,rho_array,B_x_array,B_y_array):#Right Hand Side of the y component of the momentum equation
    for j in range(2,J-2):
        if j==2 or j==J-3:
            Qd_v_y_array[j]=0
        else:
            Qd_v_y_array[j]=(1/(64*dpos))*(v_y_array[j+3]-6*v_y_array[j+2]+15*v_y_array[j+1]-20*v_y_array[j]+15*v_y_array[j-1]-6*v_y_array[j-2]+v_y_array[j-3])
        RHSMomEq_y_array[j]=((1)/(mu*rho_array[j]))*(B_x_array[j]*Finite_Difference_Derivative_Order4(B_y_array,j))-v_x_array[j]*Finite_Difference_Derivative_Order4(v_y_array,j)+sigma_v_y*Qd_v_y_array[j]
    return RHSMomEq_y_array
def RHSMomEq_z(v_x_array,v_z_array,rho_array,B_x_array,B_z_array):#Right Hand Side of the z component of the momentum equation
    for j in range(2,J-2):
        if j==2 or j==J-3:
            Qd_v_z_array[j]=0
        else:
            Qd_v_z_array[j]=(1/(64*dpos))*(v_z_array[j+3]-6*v_z_array[j+2]+15*v_z_array[j+1]-20*v_z_array[j]+15*v_z_array[j-1]-6*v_z_array[j-2]+v_z_array[j-3])
        RHSMomEq_z_array[j]=((1)/(mu*rho_array[j]))*(B_x_array[j]*Finite_Difference_Derivative_Order4(B_z_array,j))-v_x_array[j]*Finite_Difference_Derivative_Order4(v_z_array,j)+sigma_v_z*Qd_v_z_array[j]
    return RHSMomEq_z_array
def RHSIndEq_x(B_x_array):#Right Hand Side of the x component of the induction equation
    for j in range(2,J-2):
        if j==2 or j==J-3:
            Qd_B_x_array[j]=0
        else:
            Qd_B_x_array[j]=(1/(64*dpos))*(B_x_array[j+3]-6*B_x_array[j+2]+15*B_x_array[j+1]-20*B_x_array[j]+15*B_x_array[j-1]-6*B_x_array[j-2]+B_x_array[j-3])
        RHSIndEq_x_array[j]=0+sigma_B_x*Qd_B_x_array[j]
    return RHSIndEq_x_array
def RHSIndEq_y(v_x_array,v_y_array,B_x_array,B_y_array):#Right Hand Side of the y component of the induction equation
    for j in range(2,J-2):
        if j==2 or j==J-3:
            Qd_B_y_array[j]=0
        else:
            Qd_B_y_array[j]=(1/(64*dpos))*(B_y_array[j+3]-6*B_y_array[j+2]+15*B_y_array[j+1]-20*B_y_array[j]+15*B_y_array[j-1]-6*B_y_array[j-2]+B_y_array[j-3])
        RHSIndEq_y_array[j]=-Finite_Difference_Derivative_Order4(v_x_array*B_y_array-B_x_array*v_y_array,j)+sigma_B_y*Qd_B_y_array[j]
    return RHSIndEq_y_array
def RHSIndEq_z(v_x_array,v_z_array,B_x_array,B_z_array):#Right Hand Side of the z component of the induction equation
    for j in range(2,J-2):
        if j==2 or j==J-3:
            Qd_B_z_array[j]=0
        else:
            Qd_B_z_array[j]=(1/(64*dpos))*(B_z_array[j+3]-6*B_z_array[j+2]+15*B_z_array[j+1]-20*B_z_array[j]+15*B_z_array[j-1]-6*B_z_array[j-2]+B_z_array[j-3])
        RHSIndEq_z_array[j]=Finite_Difference_Derivative_Order4(v_z_array*B_x_array-B_z_array*v_x_array,j)+sigma_B_z*Qd_B_z_array[j]
    return RHSIndEq_z_array
def RHSPresEq(v_x_array,p_array):#Right Hand Side of the pressure equation
    for j in range(2,J-2):
        if j==2 or j==J-3:
            Qd_p_array[j]=0
        else:
            Qd_p_array[j]=(1/(64*dpos))*(p_array[j+3]-6*p_array[j+2]+15*p_array[j+1]-20*p_array[j]+15*p_array[j-1]-6*p_array[j-2]+p_array[j-3])
        RHSPresEq_array[j]=-gamma*p_array[j]*Finite_Difference_Derivative_Order4(v_x_array,j)-v_x_array[j]*Finite_Difference_Derivative_Order4(p_array,j)+sigma_p*Qd_p_array[j]
    return RHSPresEq_array
def RK4(v_x_array,v_y_array,v_z_array,rho_array,p_array,B_x_array,B_y_array,B_z_array):#Temporal discretization: Runge-Kutta of fourth order (RK4)
    k_1_ContEq=dt*RHSContEq(v_x_array,rho_array)
    k_1_MomEq_x=dt*RHSMomEq_x(v_x_array,rho_array,p_array,B_y_array,B_z_array)
    k_1_MomEq_y=dt*RHSMomEq_y(v_x_array,v_y_array,rho_array,B_x_array,B_y_array)
    k_1_MomEq_z=dt*RHSMomEq_z(v_x_array,v_z_array,rho_array,B_x_array,B_z_array)
    k_1_IndEq_x=dt*RHSIndEq_x(B_x_array)
    k_1_IndEq_y=dt*RHSIndEq_y(v_x_array,v_y_array,B_x_array,B_y_array)
    k_1_IndEq_z=dt*RHSIndEq_z(v_x_array,v_z_array,B_x_array,B_z_array)
    k_1_PresEq=dt*RHSPresEq(v_x_array,p_array)
    rho_i_array=rho_array+0.5*k_1_ContEq
    v_x_i_array=v_x_array+0.5*k_1_MomEq_x
    v_y_i_array=v_y_array+0.5*k_1_MomEq_y
    v_z_i_array=v_z_array+0.5*k_1_MomEq_z
    B_x_i_array=B_x_array+0.5*k_1_IndEq_x
    B_y_i_array=B_y_array+0.5*k_1_IndEq_y
    B_z_i_array=B_z_array+0.5*k_1_IndEq_z
    p_i_array=p_array+0.5*k_1_PresEq
    rho_i_array=Periodic_BC(rho_i_array)#Imposing periodic BCs
    v_x_i_array=Periodic_BC(v_x_i_array)#Imposing periodic BCs
    v_y_i_array=Periodic_BC(v_y_i_array)#Imposing periodic BCs
    v_z_i_array=Periodic_BC(v_z_i_array)#Imposing periodic BCs
    B_x_i_array=Periodic_BC(B_x_i_array)#Imposing periodic BCs
    B_y_i_array=Periodic_BC(B_y_i_array)#Imposing periodic BCs
    B_z_i_array=Periodic_BC(B_z_i_array)#Imposing periodic BCs
    p_i_array=Periodic_BC(p_i_array)#Imposing periodic BCs
    k_2_ContEq=dt*RHSContEq(v_x_i_array,rho_i_array)
    k_2_MomEq_x=dt*RHSMomEq_x(v_x_i_array,rho_i_array,p_i_array,B_y_i_array,B_z_i_array)
    k_2_MomEq_y=dt*RHSMomEq_y(v_x_i_array,v_y_i_array,rho_i_array,B_x_i_array,B_y_i_array)
    k_2_MomEq_z=dt*RHSMomEq_z(v_x_i_array,v_z_i_array,rho_i_array,B_x_i_array,B_z_i_array)
    k_2_IndEq_x=dt*RHSIndEq_x(B_x_i_array)
    k_2_IndEq_y=dt*RHSIndEq_y(v_x_i_array,v_y_i_array,B_x_i_array,B_y_i_array)
    k_2_IndEq_z=dt*RHSIndEq_z(v_x_i_array,v_z_i_array,B_x_i_array,B_z_i_array)
    k_2_PresEq=dt*RHSPresEq(v_x_i_array,p_i_array)
    rho_i_array=rho_array+0.5*k_2_ContEq
    v_x_i_array=v_x_array+0.5*k_2_MomEq_x
    v_y_i_array=v_y_array+0.5*k_2_MomEq_y
    v_z_i_array=v_z_array+0.5*k_2_MomEq_z
    B_x_i_array=B_x_array+0.5*k_2_IndEq_x
    B_y_i_array=B_y_array+0.5*k_2_IndEq_y
    B_z_i_array=B_z_array+0.5*k_2_IndEq_z
    p_i_array=p_array+0.5*k_2_PresEq
    rho_i_array=Periodic_BC(rho_i_array)#Imposing periodic BCs
    v_x_i_array=Periodic_BC(v_x_i_array)#Imposing periodic BCs
    v_y_i_array=Periodic_BC(v_y_i_array)#Imposing periodic BCs
    v_z_i_array=Periodic_BC(v_z_i_array)#Imposing periodic BCs
    B_x_i_array=Periodic_BC(B_x_i_array)#Imposing periodic BCs
    B_y_i_array=Periodic_BC(B_y_i_array)#Imposing periodic BCs
    B_z_i_array=Periodic_BC(B_z_i_array)#Imposing periodic BCs
    p_i_array=Periodic_BC(p_i_array)#Imposing periodic BCs
    k_3_ContEq=dt*RHSContEq(v_x_i_array,rho_i_array)
    k_3_MomEq_x=dt*RHSMomEq_x(v_x_i_array,rho_i_array,p_i_array,B_y_i_array,B_z_i_array)
    k_3_MomEq_y=dt*RHSMomEq_y(v_x_i_array,v_y_i_array,rho_i_array,B_x_i_array,B_y_i_array)
    k_3_MomEq_z=dt*RHSMomEq_z(v_x_i_array,v_z_i_array,rho_i_array,B_x_i_array,B_z_i_array)
    k_3_IndEq_x=dt*RHSIndEq_x(B_x_i_array)
    k_3_IndEq_y=dt*RHSIndEq_y(v_x_i_array,v_y_i_array,B_x_i_array,B_y_i_array)
    k_3_IndEq_z=dt*RHSIndEq_z(v_x_i_array,v_z_i_array,B_x_i_array,B_z_i_array)
    k_3_PresEq=dt*RHSPresEq(v_x_i_array,p_i_array)
    rho_i_array=rho_array+k_3_ContEq
    v_x_i_array=v_x_array+k_3_MomEq_x
    v_y_i_array=v_y_array+k_3_MomEq_y
    v_z_i_array=v_z_array+k_3_MomEq_z
    B_x_i_array=B_x_array+k_3_IndEq_x
    B_y_i_array=B_y_array+k_3_IndEq_y
    B_z_i_array=B_z_array+k_3_IndEq_z
    p_i_array=p_array+k_3_PresEq
    rho_i_array=Periodic_BC(rho_i_array)#Imposing periodic BCs
    v_x_i_array=Periodic_BC(v_x_i_array)#Imposing periodic BCs
    v_y_i_array=Periodic_BC(v_y_i_array)#Imposing periodic BCs
    v_z_i_array=Periodic_BC(v_z_i_array)#Imposing periodic BCs
    B_x_i_array=Periodic_BC(B_x_i_array)#Imposing periodic BCs
    B_y_i_array=Periodic_BC(B_y_i_array)#Imposing periodic BCs
    B_z_i_array=Periodic_BC(B_z_i_array)#Imposing periodic BCs
    p_i_array=Periodic_BC(p_i_array)#Imposing periodic BCs
    k_4_ContEq=dt*RHSContEq(v_x_i_array,rho_i_array)
    k_4_MomEq_x=dt*RHSMomEq_x(v_x_i_array,rho_i_array,p_i_array,B_y_i_array,B_z_i_array)
    k_4_MomEq_y=dt*RHSMomEq_y(v_x_i_array,v_y_i_array,rho_i_array,B_x_i_array,B_y_i_array)
    k_4_MomEq_z=dt*RHSMomEq_z(v_x_i_array,v_z_i_array,rho_i_array,B_x_i_array,B_z_i_array)
    k_4_IndEq_x=dt*RHSIndEq_x(B_x_i_array)
    k_4_IndEq_y=dt*RHSIndEq_y(v_x_i_array,v_y_i_array,B_x_i_array,B_y_i_array)
    k_4_IndEq_z=dt*RHSIndEq_z(v_x_i_array,v_z_i_array,B_x_i_array,B_z_i_array)
    k_4_PresEq=dt*RHSPresEq(v_x_i_array,p_i_array)
    drhodt_array=rho_array+((1)/(6))*(k_1_ContEq+2*k_2_ContEq+2*k_3_ContEq+k_4_ContEq)
    dv_xdt_array=v_x_array+((1)/(6))*(k_1_MomEq_x+2*k_2_MomEq_x+2*k_3_MomEq_x+k_4_MomEq_x)
    dv_ydt_array=v_y_array+((1)/(6))*(k_1_MomEq_y+2*k_2_MomEq_y+2*k_3_MomEq_y+k_4_MomEq_y)
    dv_zdt_array=v_z_array+((1)/(6))*(k_1_MomEq_z+2*k_2_MomEq_z+2*k_3_MomEq_z+k_4_MomEq_z)
    dB_xdt_array=B_x_array+((1)/(6))*(k_1_IndEq_x+2*k_2_IndEq_x+2*k_3_IndEq_x+k_4_IndEq_x)
    dB_ydt_array=B_y_array+((1)/(6))*(k_1_IndEq_y+2*k_2_IndEq_y+2*k_3_IndEq_y+k_4_IndEq_y)
    dB_zdt_array=B_z_array+((1)/(6))*(k_1_IndEq_z+2*k_2_IndEq_z+2*k_3_IndEq_z+k_4_IndEq_z)
    dpdt_array=p_array+((1)/(6))*(k_1_PresEq+2*k_2_PresEq+2*k_3_PresEq+k_4_PresEq)
    return [drhodt_array,dv_xdt_array,dv_ydt_array,dv_zdt_array,dB_xdt_array,dB_ydt_array,dB_zdt_array,dpdt_array]
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Solution evolution and plotting of the results
v_x_array=InitialConditionVelocity_x(pos_array)
v_y_array=InitialConditionVelocity_y(pos_array)
v_z_array=InitialConditionVelocity_z(pos_array)
p_array=InitialConditionPressure(pos_array)
rho_array=InitialConditionDensity(pos_array)
B_x_array=InitialConditionMagneticField_x(pos_array)
B_y_array=InitialConditionMagneticField_y(pos_array)
B_z_array=InitialConditionMagneticField_z(pos_array)
v_s_array=InitialConditionSoundVelocity(p_array,rho_array)#Due to the isotropic property, all the components are equal
v_A_x_array=InitialConditionAlfvenVelocity(B_x_array,rho_array)
v_A_y_array=InitialConditionAlfvenVelocity(B_y_array,rho_array)
v_A_z_array=InitialConditionAlfvenVelocity(B_z_array,rho_array)
v_s_new_array=v_s_array
v_A_x_new_array=v_A_x_array
v_A_y_new_array=v_A_y_array
v_A_z_new_array=v_A_z_array
dt=0.5*(dpos/max(3*gamma*((p_array)/(rho_array))+((1)/(mu*rho_array))*(B_x_array**2+B_y_array**2+B_z_array**2)))#Variable time-step size (Units = s). For stability purposes, we use the CFL criteria (lambda=((max(\sqrt{v_s^2+v_A^2})*dt)/dpos)<=1) to obtain this value.
N=int((tf-t0)/dt)#Time discretitzation
#Preparation of the figures' environments
#v_x(x,t) vs x
fig_v_x,ax_v_x=plt.subplots(figsize=(9,5))
plt.xlim(pos0,posf)
plt.title(r"Solution of $v_{x}$ in the 1D MHD equations")
plt.xlabel(r"$x$")
plt.ylabel(r"$v_{x}(x,t)$")
plt.grid()
fig_v_x.set_tight_layout(True)
#v_y(y,t) vs y
fig_v_y,ax_v_y=plt.subplots(figsize=(9,5))
plt.xlim(pos0,posf)
plt.title(r"Solution of $v_{y}$ in the 1D MHD equations")
plt.xlabel(r"$y$")
plt.ylabel(r"$v_{y}(y,t)$")
plt.grid()
fig_v_y.set_tight_layout(True)
#v_z(z,t) vs z
fig_v_z,ax_v_z=plt.subplots(figsize=(9,5))
plt.xlim(pos0,posf)
plt.title(r"Solution of $v_{z}$ in the 1D MHD equations")
plt.xlabel(r"$z$")
plt.ylabel(r"$v_{z}(z,t)$")
plt.grid()
fig_v_z.set_tight_layout(True)
#\rho(x,t) vs x
fig_rho,ax_rho=plt.subplots(figsize=(9,5))
plt.xlim(pos0,posf)
plt.ylim(max(rho_array)-0.5*Amp_v,max(rho_array)+0.5*Amp_v)
plt.title(r"Solution of $\rho$ in the 1D MHD equations")
plt.xlabel(r"$x$")
plt.ylabel(r"$\rho(x,t)$")
plt.grid()
fig_rho.set_tight_layout(True)
#p(x,t) vs x
fig_p,ax_p=plt.subplots(figsize=(9,5))
plt.xlim(pos0,posf)
plt.ylim(max(p_array)-Amp_v,max(p_array)+Amp_v)
plt.title(r"Solution of $p$ in the 1D MHD equations")
plt.xlabel(r"$x$")
plt.ylabel(r"$p(x,t)$")
plt.grid()
fig_p.set_tight_layout(True)
#B_x(x,t) vs x
fig_B_x,ax_B_x=plt.subplots(figsize=(9,5))
plt.xlim(pos0,posf)
plt.title(r"Solution of $B_{x}$ in the 1D MHD equations")
plt.xlabel(r"$x$")
plt.ylabel(r"$B_{x}(x,t)$")
plt.grid()
fig_B_x.set_tight_layout(True)
#B_y(y,t) vs y
fig_B_y,ax_B_y=plt.subplots(figsize=(9,5))
plt.xlim(pos0,posf)
plt.title(r"Solution of $B_{y}$ in the 1D MHD equations")
plt.xlabel(r"$y$")
plt.ylabel(r"$B_{y}(y,t)$")
plt.grid()
fig_B_y.set_tight_layout(True)
#B_z(z,t) vs z
fig_B_z,ax_B_z=plt.subplots(figsize=(9,5))
plt.xlim(pos0,posf)
plt.title(r"Solution of $B_{z}$ in the 1D MHD equations")
plt.xlabel(r"$z$")
plt.ylabel(r"$B_{z}(z,t)$")
plt.grid()
fig_B_z.set_tight_layout(True)
#v_s(x,t) vs x
fig_v_s,ax_v_s=plt.subplots(figsize=(9,5))
plt.xlim(pos0,posf)
plt.title(r"Solution of $v_{s}$ in the 1D MHD equations")
plt.xlabel(r"$x$")
plt.ylabel(r"$v_{s}(x,t)$")
plt.grid()
fig_v_s.set_tight_layout(True)
#v_A_x(x,t) vs x
fig_v_A_x,ax_v_A_x=plt.subplots(figsize=(9,5))
plt.xlim(pos0,posf)
plt.title(r"Solution of $v_{A,x}$ in the 1D MHD equations")
plt.xlabel(r"$x$")
plt.ylabel(r"$v_{A,x}(x,t)$")
plt.grid()
fig_v_A_x.set_tight_layout(True)
#v_A_y(y,t) vs y
fig_v_A_y,ax_v_A_y=plt.subplots(figsize=(9,5))
plt.xlim(pos0,posf)
plt.title(r"Solution of $v_{A,y}$ in the 1D MHD equations")
plt.xlabel(r"$y$")
plt.ylabel(r"$v_{A,y}(y,t)$")
plt.grid()
fig_v_A_y.set_tight_layout(True)
#v_A_z(z,t) vs z
fig_v_A_z,ax_v_A_z=plt.subplots(figsize=(9,5))
plt.xlim(pos0,posf)
plt.title(r"Solution of $v_{A,z}$ in the 1D MHD equations")
plt.xlabel(r"$z$")
plt.ylabel(r"$v_{A,z}(z,t)$")
plt.grid()
fig_v_A_z.set_tight_layout(True)
for n in range(N+1):
    if n % dPlot==0:
        #v_x(x,t) vs x
        frame=ax_v_x.plot(pos_array,v_x_array,color="blue")
        animation_frames_v_x_list.append(frame)
        #v_y(y,t) vs y
        frame=ax_v_y.plot(pos_array,v_y_array,color="blue")
        animation_frames_v_y_list.append(frame)
        #v_z(z,t) vs z
        frame=ax_v_z.plot(pos_array,v_z_array,color="blue")
        animation_frames_v_z_list.append(frame)
        #\rho(x,t) vs x
        frame=ax_rho.plot(pos_array,rho_array,color="blue")
        animation_frames_rho_list.append(frame)
        #p(x,t) vs x
        frame=ax_p.plot(pos_array,p_array,color="blue")
        animation_frames_p_list.append(frame)
        #B_x(x,t) vs x
        frame=ax_B_x.plot(pos_array,B_x_array,color="blue")
        animation_frames_B_x_list.append(frame)
        #B_y(y,t) vs y
        frame=ax_B_y.plot(pos_array,B_y_array,color="blue")
        animation_frames_B_y_list.append(frame)
        #B_z(z,t) vs z
        frame=ax_B_z.plot(pos_array,B_z_array,color="blue")
        animation_frames_B_z_list.append(frame)
        #v_s(x,t) vs x
        frame=ax_v_s.plot(pos_array,v_s_array,color="blue")
        animation_frames_v_s_list.append(frame)
        #v_A_x(x,t) vs x
        frame=ax_v_A_x.plot(pos_array,v_A_x_array,color="blue")
        animation_frames_v_A_x_list.append(frame)
        #v_A_y(y,t) vs y
        frame=ax_v_A_y.plot(pos_array,v_A_y_array,color="blue")
        animation_frames_v_A_y_list.append(frame)
        #v_A_z(z,t) vs z
        frame=ax_v_A_z.plot(pos_array,v_A_z_array,color="blue")
        animation_frames_v_A_z_list.append(frame)
        print("Frame ",Num_frame," computed at t = ",str(n*dt)," s.",sep="")
        Num_frame+=1
    RK4_list=RK4(v_x_array,v_y_array,v_z_array,rho_array,p_array,B_x_array,B_y_array,B_z_array)#Time evolution for the differential equations
    rho_new_array=RK4_list[0]
    v_x_new_array=RK4_list[1]
    v_y_new_array=RK4_list[2]
    v_z_new_array=RK4_list[3]
    B_x_new_array=RK4_list[4]
    B_y_new_array=RK4_list[5]
    B_z_new_array=RK4_list[6]
    p_new_array=RK4_list[7]
    rho_new_array=Periodic_BC(rho_new_array)#Imposing BCs
    v_x_new_array=Periodic_BC(v_x_new_array)#Imposing BCs
    v_y_new_array=Periodic_BC(v_y_new_array)#Imposing BCs
    v_z_new_array=Periodic_BC(v_z_new_array)#Imposing BCs
    B_x_new_array=Periodic_BC(B_x_new_array)#Imposing BCs
    B_y_new_array=Periodic_BC(B_y_new_array)#Imposing BCs
    B_z_new_array=Periodic_BC(B_z_new_array)#Imposing BCs
    p_new_array=Periodic_BC(p_new_array)#Imposing BCs
    v_s_new_array=v_s(p_new_array,rho_new_array)#"Time" evolution of v_s (due to the isotropic property, all the components are equal)
    v_A_x_new_array=v_A(B_x_new_array,rho_new_array)#"Time" evolution of the x component of v_A
    v_A_y_new_array=v_A(B_y_new_array,rho_new_array)#"Time" evolution of the y component of v_A
    v_A_z_new_array=v_A(B_z_new_array,rho_new_array)#"Time" evolution of the z component of v_A
    rho_array=rho_new_array#Update for the next time step
    v_x_array=v_x_new_array#Update for the next time step
    v_y_array=v_y_new_array#Update for the next time step
    v_z_array=v_z_new_array#Update for the next time step
    B_x_array=B_x_new_array#Update for the next time step
    B_y_array=B_y_new_array#Update for the next time step
    B_z_array=B_z_new_array#Update for the next time step
    p_array=p_new_array#Update for the next time step
    v_s_array=v_s_new_array#Update for the next time step
    v_A_x_array=v_A_x_new_array#Update for the next time step
    v_A_y_array=v_A_y_new_array#Update for the next time step
    v_A_z_array=v_A_z_new_array#Update for the next time step
    dt=0.5*(dpos/max(3*gamma*((p_array)/(rho_array))+((1)/(mu*rho_array))*(B_x_array**2+B_y_array**2+B_z_array**2)))#Variable time-step size (Units = s). For stability purposes, we use the CFL criteria (lambda=((max(\sqrt{v_s^2+v_A^2})*dt)/dpos)<=1) to obtain this value.
    N=int((tf-t0)/dt)#Time discretitzation
print("All frames computed, creating animations...")
animation_v_x=animation.ArtistAnimation(fig=fig_v_x,artists=animation_frames_v_x_list,interval=200)#Animation of the x component of the velocity results
animation_v_x.save(filename="1D MHD Solution for v_x.mp4",writer="ffmpeg")#Animation save of the x component of the velocity results
animation_v_y=animation.ArtistAnimation(fig=fig_v_y,artists=animation_frames_v_y_list,interval=200)#Animation of the y component of the velocity results
animation_v_y.save(filename="1D MHD Solution for v_y.mp4",writer="ffmpeg")#Animation save of the y component of the velocity results
animation_v_z=animation.ArtistAnimation(fig=fig_v_z,artists=animation_frames_v_z_list,interval=200)#Animation of the z component of the velocity results
animation_v_z.save(filename="1D MHD Solution for v_z.mp4",writer="ffmpeg")#Animation save of the z component of the velocity results
animation_rho=animation.ArtistAnimation(fig=fig_rho,artists=animation_frames_rho_list,interval=200)#Animation of the density results
animation_rho.save(filename="1D MHD Solution for rho.mp4",writer="ffmpeg")#Animation save of the density results
animation_B_x=animation.ArtistAnimation(fig=fig_B_x,artists=animation_frames_B_x_list,interval=200)#Animation of the x component of the magnetic field results
animation_B_x.save(filename="1D MHD Solution for B_x.mp4",writer="ffmpeg")#Animation save of the x component of the magnetic field results
animation_B_y=animation.ArtistAnimation(fig=fig_B_y,artists=animation_frames_B_y_list,interval=200)#Animation of the y component of the magnetic field results
animation_B_y.save(filename="1D MHD Solution for B_y.mp4",writer="ffmpeg")#Animation save of the y component of the magnetic field results
animation_B_z=animation.ArtistAnimation(fig=fig_B_z,artists=animation_frames_B_z_list,interval=200)#Animation of the z component of the magnetic field results
animation_B_z.save(filename="1D MHD Solution for B_z.mp4",writer="ffmpeg")#Animation save of the z component of the magnetic field results
animation_p=animation.ArtistAnimation(fig=fig_p,artists=animation_frames_p_list,interval=200)#Animation of the pressure results
animation_p.save(filename="1D MHD Solution for p.mp4",writer="ffmpeg")#Animation save of the pressure results
animation_v_s=animation.ArtistAnimation(fig=fig_v_s,artists=animation_frames_v_s_list,interval=200)#Animation of the sound velocity results (due to the isotropic property, all the components are equal)
animation_v_s.save(filename="1D MHD Solution for v_s.mp4",writer="ffmpeg")#Animation save of the sound velocity results (due to the isotropic property, all the components are equal)
animation_v_A_x=animation.ArtistAnimation(fig=fig_v_A_x,artists=animation_frames_v_A_x_list,interval=200)#Animation of the x component of the Alfvén velocity results
animation_v_A_x.save(filename="1D MHD Solution for v_A_x.mp4",writer="ffmpeg")#Animation save of the x component of the Alfvén velocity results
animation_v_A_y=animation.ArtistAnimation(fig=fig_v_A_y,artists=animation_frames_v_A_y_list,interval=200)#Animation of the y component of the Alfvén velocity results
animation_v_A_y.save(filename="1D MHD Solution for v_A_y.mp4",writer="ffmpeg")#Animation save of the y component of the Alfvén velocity results
animation_v_A_z=animation.ArtistAnimation(fig=fig_v_A_z,artists=animation_frames_v_A_z_list,interval=200)#Animation of the z component of the Alfvén velocity results
animation_v_A_z.save(filename="1D MHD Solution for v_A_z.mp4",writer="ffmpeg")#Animation save of the z component of the Alfvén velocity results
plt.close()
print("Animations saved at ",os.getcwd(),".",sep="")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")