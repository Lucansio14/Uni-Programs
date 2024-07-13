# coding: utf-8
"""
Created on Wed Sep 16 13:21:32 2020

@author: Lucas Romero Fernández
"""
# Ballistics (with friction)
# 
# In the following program, the equations of ballistics (taking into account the forces of friction in a simple manner), 
# also known as the equations of parabolic motion, are integrated in a simple model with different variable parameters to 
# adjust that make possible to calculate almost every essential scenario in ballistics.
# 
# Since this projectile does not make contact with the ground in the majority of the movement, it will be assumed that the 
# only force of friction that exists is that of the medium that is travelling through (which is proportional to the speed of
# the projectile, contrary in direction to the movement of it and of value $F_{Fr}=-knv$ (or $F_{Fr}=-mbv$), which is 
# similar to the study of fluid dynamics. These expressions can be used because of the self-imposed relatively small speed
# of the projectile.
# 
# Equations of motion:
# 
# $x(t)=x_0+\frac{v_0x}{b}*(1-e^{-b*t})$ (x-axis)
# $y(t)=y_0+\frac{1}{b}*(\frac{g}{b}+v_0y)*(1-e^{-b*t})-\frac{g}{b}*t$ (y-axis)
# 
# Example of initial conditions and parameters of the motion:
# (in International System units (S.I.), the angle of launch in degrees and it will be assumed that the projectile has the 
# shape of a sphere)
# 
# Shape of the projectile: $k=6*\pi*r$ (In this case, a sphere)
# 
# Friction coefficient (which mostly depends on the projectile): $b=\frac{k*n}{m}$
# 
# Radius of the projectile: $r=0.11145$
# 
# Mass of the projectile: $m=0.45$
# 
# Initial position in the x-axis: $x_0=0$
# 
# Initial position in the y-axis: $y_0=0$
# 
# Viscosity (cinematic) coefficient (which depends on the fluid): $n=1.51*10^{-5}$ (air at $293K$ and $1 atm$ of atmosferic pressure)
# 
# Force of the launching of the projectile: $F=1000$
# 
# Time to accelerate to the initial velocity: $t_c=0.01$
# 
# Angle of launch: $\theta=30$
# 
# Initial velocity: $v_0=\frac{F}{m}*t_c$
# 
# Gravity at Earth's surface: $g=9.8$
# 
# Gravity at Moon's surface: $g=1.6$
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
#main_program
def projectile(x0=0,y0=0,m=0.45,F=1000,tc=0.01,theta=30,g=9.8,r=0.11145,n=1.51*10**(-5)):
    v0=(F/m)*tc
    v0x=v0*np.cos(np.radians(theta))#x component of initial velocity
    v0y=v0*np.sin(np.radians(theta))#y component of initial velocity
    k=6*np.pi*r
    b=(k*n)/m
    tFlight=float(b**2*y0+b*v0y+g*LambertW(-(b*v0y+g)*exp(-(b**2*y0+b*v0y+g)/g)/g)+g)/(b*g)#Time of Flight (obtained by solving with Sympy y(t) for t)
    n=100#Number of points in t-axis
    t=np.linspace(0,tFlight,n)
    x=x0+(v0x/b)*(1-np.exp(-b*t))
    y=y0+(float(1/b)*(float(g/b)+v0y)*(1-np.exp(-b*t)))-(float(g/b)*t)
    y_max=y.max()#Maximum height achieved
    x_max=x.max()#Maximum distance traveled/Range
    i=np.where(y==y_max)[0][0]#Position of y_max in y array
    x_y_max=x[i]#x value where y=y_max
    return x,y,x_max,y_max,x_y_max,tFlight
#Examples
#g
plt.plot(projectile()[0],projectile()[1],label="Example on Earth(g=9.8)")
plt.plot(projectile(g=1.6)[0],projectile(g=1.6)[1],label="Example on the Moon (g=1.6)")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.grid()
plt.title("Trajectory of the projectile")
plt.tight_layout()
plt.legend()
plt.show()
print("Earth Example (g=9.8): x_max =",projectile()[2],", (x_y_max,y_max) = ",(projectile()[4],projectile()[3]),", tFlight =",projectile()[5])
print("Moon Example (g=1.6): x_max =",projectile(g=1.6)[2],", (x_y_max,y_max) = ",(projectile(g=1.6)[4],projectile(g=1.6)[3]),", tFlight =",projectile(g=1.6)[5])
#\theta
for angle in np.arange(40,52,2):
    plt.plot(projectile(theta=angle)[0],projectile(theta=angle)[1],label="Earth Example: $\\theta = $"+str(angle))
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.grid()
plt.title("Trajectory of the projectile")
plt.tight_layout()
plt.legend()
plt.show()
for angle in np.arange(40,52,2):
    print("\\theta =",angle)
    print("x_max =",projectile(theta=angle)[2],", (x_y_max,y_max) = ",(projectile(theta=angle)[4],projectile(theta=angle)[3]),", tFlight =",projectile(theta=angle)[5])
#y0
for height in [0,1,2.5,5,10]:
    plt.plot(projectile(y0=height)[0],projectile(y0=height)[1],label="$y_{0} = $"+str(height))
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.grid()
plt.title("Trajectory of the projectile")
plt.tight_layout()
plt.legend()
plt.show()
for height in [0,1,2.5,5,10]:
    print("y_0 = ",height)
    print("x_max =",projectile(y0=height)[2],", (x_y_max,y_max) = ",(projectile(y0=height)[4],projectile(y0=height)[3]),", tFlight =",projectile(y0=height)[5])
#F
for force in [300,500,600,1000,2000]:
    plt.plot(projectile(F=force)[0],projectile(F=force)[1],label="$F = $"+str(force))
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.grid()
plt.title("Trajectory of the projectile")
plt.tight_layout()
plt.legend()
plt.show()
for force in [300,500,600,1000,2000]:
    print("$F = $",force)
    print("x_max =",projectile(F=force)[2],", (x_y_max,y_max) = ",(projectile(F=force)[4],projectile(F=force)[3]),", tFlight =",projectile(F=force)[5])
#tc
for tc in [0.005,0.01,0.015,0.02,0.025]:
    plt.plot(projectile(tc=tc)[0],projectile(tc=tc)[1],label="$t_{c} = $"+str(tc))
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.grid()
plt.title("Trajectory of the projectile")
plt.tight_layout()
plt.legend()
plt.show()
for tc in [0.005,0.01,0.015,0.02,0.025]:
    print("t_c = ",tc)
    print("x_max =",projectile(tc=tc)[2],", (x_y_max,y_max) = ",(projectile(tc=tc)[4],projectile(tc=tc)[3]),", tFlight =",projectile(tc=tc)[5])
#m
for mass in [0.2,0.3,0.5,0.7,1]:
    plt.plot(projectile(m=mass)[0],projectile(m=mass)[1],label="$m = $"+str(mass))
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.grid()
plt.title("Trajectory of the projectile")
plt.tight_layout()
plt.legend()
plt.show()
for mass in [0.2,0.3,0.5,0.7,1]:
    print("m = ",mass)
    print("x_max =",projectile(m=mass)[2],", (x_y_max,y_max) = ",(projectile(m=mass)[4],projectile(m=mass)[3]),", tFlight =",projectile(m=mass)[5])
#x0
for distance in [0,1,2.5,5,10]:
    plt.plot(projectile(x0=distance)[0],projectile(x0=distance)[1],label="$x_{0} = $"+str(distance))
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.grid()
plt.title("Trajectory of the projectile")
plt.tight_layout()
plt.legend()
plt.show()
for distance in [0,1,2.5,5,10]:
    print("x0 = ",distance)
    print("x_max =",projectile(x0=distance)[2],", (x_y_max,y_max) = ",(projectile(x0=distance)[4],projectile(x0=distance)[3]),", tFlight =",projectile(x0=distance)[5])
#r
for radius in [0.11145,5,10,15,20]:
    plt.plot(projectile(r=radius)[0],projectile(r=radius)[1],label="$r = $"+str(radius))
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.grid()
plt.title("Trajectory of the projectile")
plt.tight_layout()
plt.legend()
plt.show()
for radius in [0.11145,5,10,15,20]:
    print("r = ",radius)
    print("x_max =",projectile(r=radius)[2],", (x_y_max,y_max) = ",(projectile(r=radius)[4],projectile(r=radius)[3]),", tFlight =",projectile(r=radius)[5])
#n
for viscosity in [1.33*10**(-5),1.51*10**(-5),1.60*10**(-5),1.69*10**(-5),1.011*10**(-6)]:
    plt.plot(projectile(n=viscosity)[0],projectile(n=viscosity)[1],label="$n = $"+str(viscosity))
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.grid()
plt.title("Trajectory of the projectile")
plt.tight_layout()
plt.legend()
plt.show()
for viscosity in [1.33*10**(-5),1.51*10**(-5),1.60*10**(-5),1.69*10**(-5),1.011*10**(-6)]:#Air:$293K$,$303K$,$313K$;Water:$293K$ (respectively).
    print("n = ",viscosity)
    print("x_max =",projectile(n=viscosity)[2],", (x_y_max,y_max) = ",(projectile(n=viscosity)[4],projectile(n=viscosity)[3]),", tFlight =",projectile(n=viscosity)[5])
