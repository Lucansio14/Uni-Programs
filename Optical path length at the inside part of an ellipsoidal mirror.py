# -*- coding: utf-8 -*-
"""
Created on Wed Oct 2 11:55:48 2024
Optical path length of a ray that undergoes reflection at the inside part of an ellipsoidal mirror.
Context: At the inside part of an ellipsoidal mirror, considering n=1 (refraction index equals 1 (void)) 
and an homogenous medium, a light ray coming from point A, with coordinates (xA,yA), bounces off the 
mirror's surface and passes through point B, with coordinates (−xA,yA).
Example of ellipsoidal mirror: https://www.meetoptics.com/assets/ellipsoidal-mirrors-C2dpGKMK.svg
@author: Lucas Romero Fernández
"""
import matplotlib.pyplot as plt
import numpy as np
#Definitions of functions and the values of parameters
def y(a,b,x):#y as a function of x for points on an ellipse with semi-major and semi-minor axes a and b, respectively
	return b*np.sqrt(1-x**2/a**2)
def L(x,a,b,xA,yA):#Optical path length
	return np.sqrt((x-xA)**2+(y(a,b,x)-yA)**2)+np.sqrt((x+xA)**2+(y(a,b,x)-yA)**2)
#main_program
#Parameter values
a=50#Length of the ellipse semi-major axis (units=mm)
b=30#Length of the ellipse semi-minor axis (units=mm)
xA=-30#x-coordinate of left point (units=mm)
yA=0#y-coordinate of both points (units=mm)
N=301#Number of points in the x-axis
#Vectors of evenly spaced x-values between -a and +a and corresponding optical path length
x_vec=np.linspace(-a,a,N)
L_vec=L(x_vec,a,b,xA,yA)
#Plot
#L (with respect to alpha) vs x
plt.figure(figsize=(8,10))
plt.plot(x_vec,L_vec,'ok',markersize=3)
plt.axvline(-a,color="blue",linestyle="dashed",alpha=0.7,label=r'$-a,+a$')
plt.axvline(a,color="blue",linestyle="dashed",alpha=0.7)
plt.xlabel(r'$x$ [mm]')
plt.ylabel(r'$L$ [mm]')
plt.title(r'Optical path length at the inside part of an ellipsoidal mirror ($x_A$ ='+str(xA)+' mm'+r', $y_A$ = '+str(yA)+' mm)')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()
