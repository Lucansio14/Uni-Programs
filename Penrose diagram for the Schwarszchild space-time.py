# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 13:23:47 2025
Penrose diagram for the, in essence, Schwarszchild space-time 1+1 (t,r) using the flat adimensional Kruskal coordinate system (T,X).
The relation between (t,r) and (T,X) in a simple way has been obtained with the help of: https://crul.github.io/CursoRelatividadGeneralJavierGarcia/docs/Miguel%20Cañizares%20-%2048%20Resumen-Diagrama%20Penrose%20Kruskal.pdf (beware, it is in spanish)
The notation of this program and the procedure following that can be found in: https://cosmo.nyu.edu/yacine/teaching/GR_2018/lectures/lecture24.pdf
Thus, refer to this website for further information.
@author: Lucas Romero Fernández
"""
import numpy as np
import matplotlib.pyplot as plt
#Definition of initial parameters
plot_scale=np.pi#Plot limits
dr=0.1#Space-step size
dt=0.3#Time-step size
rf=15#End of the space domain
tf=9#End of the time domain
#main_program
#Plot creation and LaTeX implementation
plt.rcParams.update({"text.usetex":True,"font.family":"serif","font.serif":"DejaVu Serif","mathtext.fontset":"cm","font.size":12}) #Use Latex in plots
fig,ax=plt.subplots(figsize=(9,5),dpi=600)#Use dpi=600 for superb image quality
#Zone I: Right "diamond"
#r=const orbits plotting (represented in blue)
for r in np.arange(1.0000001,rf+dr,dr):
    tilde_T_list=[]#Auxiliar time variable for the transformation to a limited space-time ("Plottable" infinities)
    tilde_X_list=[]#Auxiliar space variable for the transformation to a limited space-time ("Plottable" infinities)
    for t in np.arange(-tf,tf+dt,dt):
        T=np.sqrt(r-1)*np.exp(0.5*r)*np.sinh(0.5*t)#Change to Kruskal coordinates (depends on the region)
        X=np.sqrt(r-1)*np.exp(0.5*r)*np.cosh(0.5*t)#Change to Kruskal coordinates (depends on the region)
        tilde_T=np.arctan(T+X)+np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_X=np.arctan(T+X)-np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_T_list.append(tilde_T)
        tilde_X_list.append(tilde_X)
    ax.plot(tilde_X_list,tilde_T_list,color="blue",alpha=0.7)
#t=const orbits plotting (represented in red)
for t in np.arange(-tf,tf+dt,dt):
    tilde_T_list=[]#Auxiliar time variable for the transformation to a limited space-time ("Plottable" infinities)
    tilde_X_list=[]#Auxiliar space variable for the transformation to a limited space-time ("Plottable" infinities)
    for r in np.arange(1.00001,rf+dr,dr):
        T=np.sqrt(r-1)*np.exp(0.5*r)*np.sinh(0.5*t)#Change to Kruskal coordinates (depends on the region)
        X=np.sqrt(r-1)*np.exp(0.5*r)*np.cosh(0.5*t)#Change to Kruskal coordinates (depends on the region)
        tilde_T=np.arctan(T+X)+np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_X=np.arctan(T+X)-np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_T_list.append(tilde_T)
        tilde_X_list.append(tilde_X)
    ax.plot(tilde_X_list,tilde_T_list,color="red",alpha=0.7)
#Zone II: Upper half "diamond"
#r=const orbits plotting (represented in blue)
for r in np.arange(0.999999,0.00001,-dr):
    tilde_T_list=[]#Auxiliar time variable for the transformation to a limited space-time ("Plottable" infinities)
    tilde_X_list=[]#Auxiliar space variable for the transformation to a limited space-time ("Plottable" infinities)
    for t in np.arange(-tf,tf+dt,dt):
        T=np.sqrt(1-r)*np.exp(0.5*r)*np.cosh(0.5*t)#Change to Kruskal coordinates (depends on the region)
        X=np.sqrt(1-r)*np.exp(0.5*r)*np.sinh(0.5*t)#Change to Kruskal coordinates (depends on the region)
        tilde_T=np.arctan(T+X)+np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_X=np.arctan(T+X)-np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_T_list.append(tilde_T)
        tilde_X_list.append(tilde_X)
    ax.plot(tilde_X_list,tilde_T_list,color="blue",alpha=0.7)
#t=const orbits plotting (represented in red)
for t in np.arange(-tf,tf+dt,dt):
    tilde_T_list=[]#Auxiliar time variable for the transformation to a limited space-time ("Plottable" infinities)
    tilde_X_list=[]#Auxiliar space variable for the transformation to a limited space-time ("Plottable" infinities)
    for r in np.arange(0.999999,0.00001,-dr):
        T=np.sqrt(1-r)*np.exp(0.5*r)*np.cosh(0.5*t)#Change to Kruskal coordinates (depends on the region)
        X=np.sqrt(1-r)*np.exp(0.5*r)*np.sinh(0.5*t)#Change to Kruskal coordinates (depends on the region)
        tilde_T=np.arctan(T+X)+np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_X=np.arctan(T+X)-np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_T_list.append(tilde_T)
        tilde_X_list.append(tilde_X)
    ax.plot(tilde_X_list,tilde_T_list,color="red",alpha=0.7)
#Zone III: Left "diamond"
#r=const orbits plotting (represented in blue)
for r in np.arange(1.00001,rf+dr,dr):
    tilde_T_list=[]#Auxiliar time variable for the transformation to a limited space-time ("Plottable" infinities)
    tilde_X_list=[]#Auxiliar space variable for the transformation to a limited space-time ("Plottable" infinities)
    for t in np.arange(-tf,tf+dt,dt):
        T=-np.sqrt(r-1)*np.exp(0.5*r)*np.sinh(-0.5*t)#Change to Kruskal coordinates (depends on the region)
        X=-np.sqrt(r-1)*np.exp(0.5*r)*np.cosh(-0.5*t)#Change to Kruskal coordinates (depends on the region)
        tilde_T=np.arctan(T+X)+np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_X=np.arctan(T+X)-np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_T_list.append(tilde_T)
        tilde_X_list.append(tilde_X)
    ax.plot(tilde_X_list,tilde_T_list,color="blue",alpha=0.7)
#t=const orbits plotting (represented in red)
for t in np.arange(-tf,tf+dt,dt):
    tilde_T_list=[]#Auxiliar time variable for the transformation to a limited space-time ("Plottable" infinities)
    tilde_X_list=[]#Auxiliar space variable for the transformation to a limited space-time ("Plottable" infinities)
    for r in np.arange(1.00001,rf+dr,dr):
        T=-np.sqrt(r-1)*np.exp(0.5*r)*np.sinh(-0.5*t)#Change to Kruskal coordinates (depends on the region)
        X=-np.sqrt(r-1)*np.exp(0.5*r)*np.cosh(-0.5*t)#Change to Kruskal coordinates (depends on the region)
        tilde_T=np.arctan(T+X)+np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_X=np.arctan(T+X)-np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_T_list.append(tilde_T)
        tilde_X_list.append(tilde_X)
    ax.plot(tilde_X_list,tilde_T_list,color="red",alpha=0.7)
#Zone IV: Lower half "diamond"
#r=const orbits plotting (represented in blue)
for r in np.arange(0.999999,0.00001,-dr):
    tilde_T_list=[]#Auxiliar time variable for the transformation to a limited space-time ("Plottable" infinities)
    tilde_X_list=[]#Auxiliar space variable for the transformation to a limited space-time ("Plottable" infinities)
    for t in np.arange(-tf,tf+dt,dt):
        T=-np.sqrt(1-r)*np.exp(0.5*r)*np.cosh(0.5*t)#Change to Kruskal coordinates (depends on the region)
        X=-np.sqrt(1-r)*np.exp(0.5*r)*np.sinh(0.5*t)#Change to Kruskal coordinates (depends on the region)
        tilde_T=np.arctan(T+X)+np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_X=np.arctan(T+X)-np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_T_list.append(tilde_T)
        tilde_X_list.append(tilde_X)
    ax.plot(tilde_X_list,tilde_T_list,color="blue",alpha=0.7)
#t=const orbits plotting (represented in red)
for t in np.arange(-tf,tf+dt,dt):
    tilde_T_list=[]#Auxiliar time variable for the transformation to a limited space-time ("Plottable" infinities)
    tilde_X_list=[]#Auxiliar space variable for the transformation to a limited space-time ("Plottable" infinities)
    for r in np.arange(0.999999,0.00001,-dr):
        T=-np.sqrt(1-r)*np.exp(0.5*r)*np.cosh(0.5*t)#Change to Kruskal coordinates (depends on the region)
        X=-np.sqrt(1-r)*np.exp(0.5*r)*np.sinh(0.5*t)#Change to Kruskal coordinates (depends on the region)
        tilde_T=np.arctan(T+X)+np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_X=np.arctan(T+X)-np.arctan(T-X)#Direct full transformation to a limited space-time ("Plottable" infinities)
        tilde_T_list.append(tilde_T)
        tilde_X_list.append(tilde_X)
    ax.plot(tilde_X_list,tilde_T_list,color="red",alpha=0.7)
#Singularities plotting (represented in green)
tilde_X_sing_array=np.arange(-np.pi/2,np.pi/2+0.1,0.1)#\tilde{X} values of the singularities
tilde_T_sing_upper_list=[np.pi/2]*len(tilde_X_sing_array)#\tilde{T} values of the upper singularity
tilde_T_sing_low_list=[-np.pi/2]*len(tilde_X_sing_array)#\tilde{T} values of the lower singularity
ax.plot(tilde_X_sing_array,tilde_T_sing_upper_list,color="green",marker="*")
ax.plot(tilde_X_sing_array,tilde_T_sing_low_list,color="green",marker="*")
#Infinities indicators
ax.text(3.2,-0.05,r"$i^{0}$",fontsize=14)#Spacelike
ax.text(1.55,1.7,r"$i^{+}$",fontsize=14)#Future timelike
ax.text(1.55,-1.8,r"$i^{-}$",fontsize=14)#Past timelike
ax.text(2.4,0.9,r"$\mathcal{I}^{+}$",fontsize=14)#Future null
ax.text(2.4,-1.05,r"$\mathcal{I}^{-}$",fontsize=14)#Past null
ax.text(-3.35,-0.05,r"$i^{0}$",fontsize=14)#Spacelike (extended)
ax.text(-1.6,1.7,r"$i^{+}$",fontsize=14)#Future timelike (extended)
ax.text(-1.6,-1.8,r"$i^{-}$",fontsize=14)#Past timelike (extended)
ax.text(-2.6,0.9,r"$\mathcal{I}^{+}$",fontsize=14)#Future null (extended)
ax.text(-2.6,-1.0,r"$\mathcal{I}^{-}$",fontsize=14)#Past null (extended)
#Plot configuration
plt.gca().set_aspect("equal",adjustable="box")
ax.spines["left"].set_position("zero")
ax.spines["right"].set_color("none")
ax.spines["bottom"].set_position("zero")
ax.spines["top"].set_color("none")
ax.set_xlim(-plot_scale,plot_scale)
ax.set_ylim(-plot_scale/2,plot_scale/2)
ax.set_xlabel(r"$\mathrm{\tilde{X}}$",loc="right")
ax.set_ylabel(r"$\mathrm{\tilde{T}}$",loc="top",rotation="horizontal")
ax.grid(False)
fig.set_tight_layout(True)
plt.savefig("Penrose diagram for the Schwarszchild space-time")
plt.show()