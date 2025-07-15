# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 10:21:24 2025
Solver of the 1D Euler equations: A basic 1D non linear equation system in the conservative and compressible form for pressure (p), density (\rho), velocity (v), and total energy (E), used primarly in hydrodynamics, with RK4 (time part), FD4 (spatial part) and Kreiss-Oliger dissipation filter.
The corresponding equations (considering mostly adimensional variables, and neglecting gravitational effects due to relatively small spatial domain):
    \partial_{t}\rho = -\partial_{x}(v\rho) (Continuity equation)
    \partial_{t}(v\rho) = -\partial_{x}(v^2\rho + p) ---> \partial_{t}v = -((1)/(\rho))\partial_{p} - v\partial_{x}v (Momentum equation, because momentum (m) = density (\rho) * velocity (v) (this will be mayorly used in the computations))
    \partial_{t}E = -\partial_{x}(v(E + p)) ((Adiabatic) energy equation)
    p = (\gamma - 1)(E - 0.5v^2\rho) (Ideal equation of state)
(where \gamma is the adiabatic index and, importantly for later, v_s = \sqrt{(p\gamma)/(\rho)} is the adiabatic sound velocity)
More information about these equations in these websites: https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)
                                                          https://www.clawpack.org/riemann_book/html/Euler.html
                                                          https://eprints.maths.manchester.ac.uk/150/1/JuH-EE.pdf
More information about CFL criteria for non linear equations in this website: https://scicomp.stackexchange.com/questions/31456/cfl-equation-for-non-linear-equation

Context about the computed solution (with initial perturbations in the velocity):
Without the use of more sophisticated methods (like shock capturing), the non linear case (not Amp_v<<1, such as Amp_v = 1) solutions are "corrupted" in some way by artifact spikes. The reason of this issue is completely numerical. For that, only the linear case (Amp_v<<1, such as Amp_v = 0.1) solution is correctly obtainable with this code.
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
#Definition of general constants, variables, lists, arrays and functions
dx=0.02#Space-step size
x0=-1#Start of the space domain
xf=1#End of the space domain
t0=0#Start time (Units = s)
tf=1.5#End time (Units = s) (approximately, due to variable CFL condition)
dPlot=5#For the plots (specifically, the number of iterations made between plots)
sigma_rho=1.5#Control parameter of the density Kreiss-Oliger dissipation filter term (Qd_rho) (sigma_rho >= 0 necessary!)
sigma_v=1.5#Control parameter of the velocity Kreiss-Oliger dissipation filter term (Qd_v) (sigma_v >= 0 necessary!)
sigma_E=1.5#Control parameter of the energy Kreiss-Oliger dissipation filter term (Qd_E) (sigma_E >= 0 necessary!)
gamma=5/3#Constant adiabatic index
Amp_v=10**(-3)#Amplitude of the initial velocity gaussian perturbation (for linear case, Amp_v<<1 such as Amp_v = 0.1)
J=int((xf-x0)/dx)#Number of spatial points
x_array=np.linspace(x0,xf,J)#Values of x
Qd_rho_array=np.zeros(J)#Kreiss-Oliger dissipation filter term array for the continuity equation
Qd_v_array=np.zeros(J)#Kreiss-Oliger dissipation filter term array for the momentum equation
Qd_E_array=np.zeros(J)#Kreiss-Oliger dissipation filter term array for the energy equation
RHSContEq_array=np.zeros(J)#Right Hand Side array of the continuity equation
RHSMomEq_array=np.zeros(J)#Right Hand Side array of the momentum equation
RHSEnergyEq_array=np.zeros(J)#Right Hand Side array of the energy equation
animation_frames_v_list=[]#Animation frames/artists for the velocity results
animation_frames_rho_list=[]#Animation frames/artists for the density results
animation_frames_p_list=[]#Animation frames/artists for the pressure results
animation_frames_E_list=[]#Animation frames/artists for the energy results
Num_frame=0#Number of solution frame for the animations
#Initial conditions for all the variables present in the equations (and general expression for v_s)
def ICv(x):#Initial condition for the velocity perturbation, in this case, a Gaussian function
    return Amp_v*np.exp(-(1/2)*((x-0)/(0.1))**2)
def ICp(x):#Initial condition for the pressure, in this case, unitary
    return 1
def ICrho(x):#Initial condition for the density, in this case, unitary
    return 1
def ICE(p,rho,v):#Initial condition for the energy
    return ((p)/(gamma-1))+0.5*rho*v**2
def v_s(p,rho):#Expression of the sound velocity
    return np.sqrt((p*gamma)/(rho))
def InitialConditionVelocity(x_array):#Preparation for the implementation of the initial velocity condition
    j=0
    v_array=np.zeros(J)
    for x in x_array:
        v_array[j]=ICv(x)
        j+=1
    return v_array
def InitialConditionPressure(x_array):#Preparation for the implementation of the initial pressure condition
    j=0
    p_array=np.zeros(J)
    for x in x_array:
        p_array[j]=ICp(x)
        j+=1
    return p_array
def InitialConditionDensity(x_array):#Preparation for the implementation of the initial density condition
    j=0
    rho_array=np.zeros(J)
    for x in x_array:
        rho_array[j]=ICrho(x)
        j+=1
    return rho_array
def InitialConditionEnergy(p_array,rho_array,v_array):#Preparation for the implementation of the initial energy condition
    j=0
    E_array=np.zeros(J)
    while j<J:
        E_array[j]=ICE(p_array[j],rho_array[j],v_array[j])
        j+=1
    return E_array
def InitialConditionSoundVelocity(p_array,rho_array):#Preparation for the implementation of the initial sound velocity condition
    j=0
    v_s_array=np.zeros(J)
    while j<J:
        v_s_array[j]=v_s(p_array[j],rho_array[j])
        j+=1
    return v_s_array
#Numerical methods utilized for each equation (general and specific ones)
def Finite_Difference_Derivative_Order4(u_array,j):#General spatial discretization: Finite/discrete central differences (4rth order)
    dudr=(1/(12*dx))*(-u_array[j+2]+8*u_array[j+1]-8*u_array[j-1]+u_array[j-2])
    return dudr
def Periodic_BC(u_array):#General (circular) periodic boundary conditions
    u_array[0]=u_array[-4]
    u_array[1]=u_array[-3]
    u_array[-2]=u_array[2]
    u_array[-1]=u_array[3]
    return u_array
def RHSContEq(v_array,rho_array):#Right Hand Side of the continuity equation
    for j in range(2,J-2):
        if j==2 or j==J-3:
            Qd_rho_array[j]=0
        else:
            Qd_rho_array[j]=(1/(64*dx))*(rho_array[j+3]-6*rho_array[j+2]+15*rho_array[j+1]-20*rho_array[j]+15*rho_array[j-1]-6*rho_array[j-2]+rho_array[j-3])
        RHSContEq_array[j]=-Finite_Difference_Derivative_Order4(rho_array*v_array,j)+sigma_rho*Qd_rho_array[j]
    return RHSContEq_array
def RHSMomEq(v_array,rho_array,p_array):#Right Hand Side of the momentum equation
    for j in range(2,J-2):
        if j==2 or j==J-3:
            Qd_v_array[j]=0
        else:
            Qd_v_array[j]=(1/(64*dx))*(v_array[j+3]-6*v_array[j+2]+15*v_array[j+1]-20*v_array[j]+15*v_array[j-1]-6*v_array[j-2]+v_array[j-3])
        RHSMomEq_array[j]=-((1)/(rho_array[j]))*Finite_Difference_Derivative_Order4(p_array,j)-v_array[j]*Finite_Difference_Derivative_Order4(v_array,j)+sigma_v*Qd_v_array[j]
    return RHSMomEq_array
def RHSEnergyEq(v_array,E_array,p_array):#Right Hand Side of the energy equation
    for j in range(2,J-2):
        if j==2 or j==J-3:
            Qd_E_array[j]=0
        else:
            Qd_E_array[j]=(1/(64*dx))*(E_array[j+3]-6*E_array[j+2]+15*E_array[j+1]-20*E_array[j]+15*E_array[j-1]-6*E_array[j-2]+E_array[j-3])
        RHSEnergyEq_array[j]=-Finite_Difference_Derivative_Order4(v_array*(E_array+p_array),j)+sigma_E*Qd_E_array[j]
    return RHSEnergyEq_array 
def RK4(v_array,rho_array,E_array,p_array):#Temporal discretization: Runge-Kutta of fourth order (RK4)
    k_1_ContEq=dt*RHSContEq(v_array,rho_array)
    k_1_MomEq=dt*RHSMomEq(v_array,rho_array,p_array)
    k_1_EnergyEq=dt*RHSEnergyEq(v_array,E_array,p_array)
    rho_i_array=rho_array+0.5*k_1_ContEq
    v_i_array=v_array+0.5*k_1_MomEq
    E_i_array=E_array+0.5*k_1_EnergyEq
    rho_i_array=Periodic_BC(rho_i_array)#Imposing periodic BCs
    v_i_array=Periodic_BC(v_i_array)#Imposing periodic BCs
    E_i_array=Periodic_BC(E_i_array)#Imposing periodic BCs
    p_i_array=(gamma-1)*(E_i_array-0.5*v_i_array**2*rho_i_array)
    k_2_ContEq=dt*RHSContEq(v_i_array,rho_i_array)
    k_2_MomEq=dt*RHSMomEq(v_i_array,rho_i_array,p_i_array)
    k_2_EnergyEq=dt*RHSEnergyEq(v_i_array,E_i_array,p_i_array)
    rho_i_array=rho_array+0.5*k_2_ContEq
    v_i_array=v_array+0.5*k_2_MomEq
    E_i_array=E_array+0.5*k_2_EnergyEq
    rho_i_array=Periodic_BC(rho_i_array)#Imposing periodic BCs
    v_i_array=Periodic_BC(v_i_array)#Imposing periodic BCs
    E_i_array=Periodic_BC(E_i_array)#Imposing periodic BCs
    p_i_array=(gamma-1)*(E_i_array-0.5*v_i_array**2*rho_i_array)
    k_3_ContEq=dt*RHSContEq(v_i_array,rho_i_array)
    k_3_MomEq=dt*RHSMomEq(v_i_array,rho_i_array,p_i_array)
    k_3_EnergyEq=dt*RHSEnergyEq(v_i_array,E_i_array,p_i_array)
    rho_i_array=rho_array+k_3_ContEq
    v_i_array=v_array+k_3_MomEq
    E_i_array=E_array+k_3_EnergyEq
    rho_i_array=Periodic_BC(rho_i_array)#Imposing periodic BCs
    v_i_array=Periodic_BC(v_i_array)#Imposing periodic BCs
    E_i_array=Periodic_BC(E_i_array)#Imposing periodic BCs
    p_i_array=(gamma-1)*(E_i_array-0.5*v_i_array**2*rho_i_array)
    k_4_ContEq=dt*RHSContEq(v_i_array,rho_i_array)
    k_4_MomEq=dt*RHSMomEq(v_i_array,rho_i_array,p_i_array)
    k_4_EnergyEq=dt*RHSEnergyEq(v_i_array,E_i_array,p_i_array)
    drhodt_array=rho_array+((1)/(6))*(k_1_ContEq+2*k_2_ContEq+2*k_3_ContEq+k_4_ContEq)
    dvdt_array=v_array+((1)/(6))*(k_1_MomEq+2*k_2_MomEq+2*k_3_MomEq+k_4_MomEq)
    dEdt_array=E_array+((1)/(6))*(k_1_EnergyEq+2*k_2_EnergyEq+2*k_3_EnergyEq+k_4_EnergyEq)
    return [drhodt_array,dvdt_array,dEdt_array]
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Solution evolution and plotting of the results
v_array=InitialConditionVelocity(x_array)
p_array=InitialConditionPressure(x_array)
rho_array=InitialConditionDensity(x_array)
E_array=InitialConditionEnergy(p_array,rho_array,v_array)
v_s_array=InitialConditionSoundVelocity(p_array,rho_array)
p_new_array=p_array
v_s_new_array=v_s_array
dt=0.5*(dx/max(v_s_array))#Variable time-step size (Units = s). For stability purposes, we use the CFL criteria (lambda=((max(v_s)*dt)/dx)<=1) to obtain this value.
N=int((tf-t0)/dt)#Time discretitzation
#Preparation of the figures' environments
#v(x,t) vs x
fig_v,ax_v=plt.subplots(figsize=(9,5))
plt.xlim(x0,xf)
plt.ylim(-Amp_v-0.1*Amp_v,Amp_v+0.1*Amp_v)
plt.title(r"Solution of $v$ in the 1D Euler equations")
plt.xlabel(r"$x$")
plt.ylabel(r"$v(x,t)$")
plt.grid()
fig_v.set_tight_layout(True)
#\rho(x,t) vs x
fig_rho,ax_rho=plt.subplots(figsize=(9,5))
plt.xlim(x0,xf)
plt.ylim(max(rho_array)-Amp_v,max(rho_array)+Amp_v)
plt.title(r"Solution of $\rho$ in the 1D Euler equations")
plt.xlabel(r"$x$")
plt.ylabel(r"$\rho(x,t)$")
plt.grid()
fig_rho.set_tight_layout(True)
#p(x,t) vs x
fig_p,ax_p=plt.subplots(figsize=(9,5))
plt.xlim(x0,xf)
plt.ylim(max(p_array)-2*Amp_v,max(p_array)+2*Amp_v)
plt.title(r"Solution of $p$ in the 1D Euler equations")
plt.xlabel(r"$x$")
plt.ylabel(r"$p(x,t)$")
plt.grid()
fig_p.set_tight_layout(True)
#E(x,t) vs x
fig_E,ax_E=plt.subplots(figsize=(9,5))
plt.xlim(x0,xf)
plt.ylim(max(E_array)-2*Amp_v,max(E_array)+2*Amp_v)
plt.title(r"Solution of $E$ in the 1D Euler equations")
plt.xlabel(r"$x$")
plt.ylabel(r"$E(x,t)$")
plt.grid()
fig_E.set_tight_layout(True)
for n in range(N+1):
    if n % dPlot==0:
        #v(x,t) vs x
        frame=ax_v.plot(x_array,v_array,color="blue")
        animation_frames_v_list.append(frame)
        #\rho(x,t) vs x
        frame=ax_rho.plot(x_array,rho_array,color="blue")
        animation_frames_rho_list.append(frame)
        #p(x,t) vs x
        frame=ax_p.plot(x_array,p_array,color="blue")
        animation_frames_p_list.append(frame)
        #E(x,t) vs x
        frame=ax_E.plot(x_array,E_array,color="blue")
        animation_frames_E_list.append(frame)
        print("Frame ",Num_frame," computed at t = ",str(n*dt)," s.",sep="")
        Num_frame+=1
    RK4_list=RK4(v_array,rho_array,E_array,p_array)#Time evolution for the differential equations
    rho_new_array=RK4_list[0]
    v_new_array=RK4_list[1]
    E_new_array=RK4_list[2]
    rho_new_array=Periodic_BC(rho_new_array)#Imposing BCs
    v_new_array=Periodic_BC(v_new_array)#Imposing BCs
    E_new_array=Periodic_BC(E_new_array)#Imposing BCs
    p_new_array=(gamma-1)*(E_new_array-0.5*v_new_array**2*rho_new_array)#"Time" evolution of p
    v_s_new_array=v_s(p_new_array,rho_new_array)#"Time" evolution of v_s
    rho_array=rho_new_array#Update for the next time step
    v_array=v_new_array#Update for the next time step
    E_array=E_new_array#Update for the next time step
    p_array=p_new_array#Update for the next time step
    v_s_array=v_s_new_array#Update for the next time step
    dt=0.5*(dx/max(v_s_array))#Variable time-step size (Units = s). For stability purposes, we use the CFL criteria (lambda=((max(v_s)*dt)/dx)<=1) to obtain this value.
    N=int((tf-t0)/dt)#Time discretitzation
print("All frames computed, creating animations...")
animation_v=animation.ArtistAnimation(fig=fig_v,artists=animation_frames_v_list,interval=200)#Animation of the velocity results
animation_v.save(filename="1D Euler Solution for v.mp4",writer="ffmpeg")#Animation save of the velocity results
animation_rho=animation.ArtistAnimation(fig=fig_rho,artists=animation_frames_rho_list,interval=200)#Animation of the density results
animation_rho.save(filename="1D Euler Solution for rho.mp4",writer="ffmpeg")#Animation save of the density results
animation_p=animation.ArtistAnimation(fig=fig_p,artists=animation_frames_p_list,interval=200)#Animation of the pressure results
animation_p.save(filename="1D Euler Solution for p.mp4",writer="ffmpeg")#Animation save of the pressure results
animation_E=animation.ArtistAnimation(fig=fig_E,artists=animation_frames_E_list,interval=200)#Animation of the energy results
animation_E.save(filename="1D Euler Solution for E.mp4",writer="ffmpeg")#Animation save of the energy results
plt.close()
print("Animations saved at ",os.getcwd(),".",sep="")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")