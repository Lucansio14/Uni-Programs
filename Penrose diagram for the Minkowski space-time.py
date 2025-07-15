# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 13:18:27 2025
Penrose diagram for the Minkowski space-time 1+1 (t,r) using spherical polar coordinates (T,R).
The notation and procedure followed for its obtainment can be found in: https://cosmo.nyu.edu/yacine/teaching/GR_2018/lectures/lecture24.pdf
Thus, refer to this website for further information.
@author: Lucas Romero Fernández
"""
import numpy as np
import matplotlib.pyplot as plt
#main_program
#Plot creation and LaTeX implementation
plt.rcParams.update({"text.usetex":True,"font.family":"serif","font.serif":"DejaVu Serif","mathtext.fontset":"cm","font.size":12}) #Use Latex in plots
fig,ax=plt.subplots(figsize=(9,5),dpi=600)#Use dpi=600 for superb image quality
#r=const orbits plotting (represented in blue)
for r in np.arange(-50,50.2,0.2):
    T_list=[]#Auxiliar time variable for the transformation to a limited space-time ("Plottable" infinities)
    R_list=[]#Auxiliar space variable for the transformation to a limited space-time ("Plottable" infinities)
    for t in np.arange(-100,100.2,0.2):
        T=np.arctan(t+r)+np.arctan(t-r)#Direct full transformation to a limited space-time ("Plottable" infinities)
        R=np.arctan(t+r)-np.arctan(t-r)#Direct full transformation to a limited space-time ("Plottable" infinities)
        T_list.append(T)
        R_list.append(R)
    ax.plot(R_list,T_list,color="blue",alpha=0.7)
#t=const orbits plotting (represented in red)
for t in np.arange(-50,50.2,0.2):
    T_list=[]#Auxiliar time variable for the transformation to a limited space-time ("Plottable" infinities)
    R_list=[]#Auxiliar space variable for the transformation to a limited space-time ("Plottable" infinities)
    for r in np.arange(-50,50.2,0.2):
        T=np.arctan(t+r)+np.arctan(t-r)#Direct full transformation to a limited space-time ("Plottable" infinities)
        R=np.arctan(t+r)-np.arctan(t-r)#Direct full transformation to a limited space-time ("Plottable" infinities)
        T_list.append(T)
        R_list.append(R)
    ax.plot(R_list,T_list,color="red",alpha=0.7)
#Infinities indicators
ax.text(np.pi+0.10,-0.10,r"$i^{0}$",fontsize=14)#Spacelike
ax.text(-0.05,np.pi+0.10,r"$i^{+}$",fontsize=14)#Future timelike
ax.text(-0.10,-np.pi-0.30,r"$i^{-}$",fontsize=14)#Past timelike
ax.text(np.pi/2,np.pi/2+0.20,r"$\mathcal{I}^{+}$",fontsize=14)#Future null
ax.text(np.pi/2+0.10,-np.pi/2-0.20,r"$\mathcal{I}^{-}$",fontsize=14)#Past null
ax.text(-np.pi/2-0.20,np.pi/2+0.20,r"$\mathcal{I}^{+}$",fontsize=14)#Future timelike (extended)
ax.text(-np.pi/2-0.40,-np.pi/2-0.20,r"$\mathcal{I}^{-}$",fontsize=14)#Past timelike (extended)
ax.text(-np.pi-0.30,-0.10,r"$i^{0}$",fontsize=14)#Spacelike (extended)
#Plot configuration
plt.gca().set_aspect("equal",adjustable="box")
ax.spines["left"].set_position("zero")
ax.spines["right"].set_color("none")
ax.spines["bottom"].set_position("zero")
ax.spines["top"].set_color("none")
ax.set_xlim(-np.pi,np.pi)
ax.set_ylim(-np.pi,np.pi)
ax.set_xlabel(r"$\mathrm{R}$",loc="right")
ax.set_ylabel(r"$\mathrm{T}$",loc="top",rotation="horizontal")
ax.grid(False)
fig.set_tight_layout(True)
plt.savefig("Penrose diagram for the Minkowski space-time")
plt.show()