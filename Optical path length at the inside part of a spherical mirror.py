# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 13:15:23 2024
Optical path length of a ray that undergoes reflection at the inside part of a spherical mirror.
Context: At the inside part of a spherical mirror, considering n=1 (refraction index equals 1 (void)) 
and an homogenous medium, a light ray coming from point A, with coordinates (xA,yA), bounces off the 
mirror's surface and passes through point B, with coordinates (−xA,yA).
Example of spherical mirror: https://avantierinc.com/wp-content/uploads/2024/02/Figure-2.-Pole-of-Spherical-Mirrors-768x343.png
@author: Lucas Romero Fernández
"""
import matplotlib.pyplot as plt
import numpy as np
#Definitions of functions and the values of parameters
def y(R,x):#y as a function of x for points on a circumference of radius R
	return np.sqrt(R**2-x**2)
def L(x,R,xA,yA):#Optical path length
	return np.sqrt((x-xA)**2+(y(R,x)-yA)**2)+np.sqrt((x+xA)**2+(y(R,x)-yA)**2)
#main_program
#Parameter values
R=50#Mirror radius (units=mm)
xA=-30#x-coordinate of left point (units=mm)
yA=10#y-coordinate of both points (units=mm)
N=301#Number of points in the x-axis
#Vectors of evenly spaced x-values between -R and +R and corresponding optical path length
x_vec=np.linspace(-R,R,N)
L_vec=L(x_vec,R,xA,yA)
#To calculate the incidence and reflection angles
aA=np.sqrt((x_vec-xA)**2+(y(R,x_vec)-yA)**2)
aB=np.sqrt((x_vec+xA)**2+(y(R,x_vec)-yA)**2)
b=R
c=np.sqrt(xA**2+yA**2)
thetai_vec=np.arccos((aA**2+b**2-c**2)/(2*aA*b))#Incidence angle
thetat_vec=np.arccos((aB**2+b**2-c**2)/(2*aB*b))#Reflection angle
#Plots
#L (with respect to alpha) vs x
plt.figure(figsize=(8,10))
plt.subplot(2,1,1)
plt.plot(x_vec,L_vec,'ok',markersize=3)
plt.axvline(-R,color="blue",linestyle="dashed",alpha=0.7,label=r'$-R,+R$')
plt.axvline(R,color="blue",linestyle="dashed",alpha=0.7)
plt.xlabel(r'$x$ [mm]')
plt.ylabel(r'$L$ [mm]')
plt.title(r'Optical path length at the inside part of a spherical mirror ($x_A$ = '+str(xA)+' mm'+r', $y_A$ = '+str(yA)+' mm)')
plt.grid()
plt.legend()
plt.tight_layout()
#Incidence and reflection angles vs x
plt.subplot(2,1,2)
plt.plot(x_vec,thetai_vec,"-k",label=r'$\theta_i$')
plt.plot(x_vec,thetat_vec,"-r",label=r'$\theta_t$')
plt.axvline(-R,color="blue",linestyle="dashed",alpha=0.7,label=r'$-R,+R$')
plt.axvline(R,color="blue",linestyle="dashed",alpha=0.7)
plt.xlabel(r'$x$ [mm]')
plt.ylabel(r'$\theta_i, \theta_t$ [rad]')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()
