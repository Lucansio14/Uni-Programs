# coding: utf-8
"""
Created on Mon Sep 18 13:23:46 2020

@author: Lucas Romero Fernández
"""
# Frictionless Ballistics
# 
# In the following program, the equations of frictionless ballistics, also known as the equations of parabolic motion, are
# integrated with different variable parameters to adjust that make possible to calculate almost every essential scenario in
# frictionless ballistics.
#
# Equations of motion:
#
# $x(t)=x_0+v_0*cos(\theta)*t$ (x-axis)
# $y(t)=y_0+v_0*sin(\theta)*t-\frac{1}{2}*g*t^2$ (y-axis)
# 
# Example of initial conditions and parameters of the motion: 
# (in International System units (S.I.) and the angle of launch in degrees)
# 
# Initial position in the x-axis: $x_0=0$
# 
# Initial position in the y-axis: $y_0=0$ 
# 
# Mass of the projectile: $m=0.3$
# 
# Force of the launching of the projectile: $F=500$
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
import numpy as np
import matplotlib.pyplot as plt
#main_program
def projectile(x0=0,y0=0,m=0.3,F=500,tc=0.01,theta=30,g=9.8):
    v0=(F/m)*tc
    v0x=v0*np.cos(np.radians(theta))#x component of initial velocity
    v0y=v0*np.sin(np.radians(theta))#y component of initial velocity
    tFlight=(v0y+np.sqrt(2*g*y0+v0y**2))/(g)#Time of Flight
    n=100#Number of points in t-axis
    t=np.linspace(0,tFlight,n)
    x=x0+v0x*t
    y=y0+v0y*t-(1/2)*g*t**2
    y_max=y.max()#Maximum height achieved
    x_max=x.max()#Maximum distance traveled/Range
    t_y_max=(v0y+np.sqrt((v0y**(2))-2*g*(y_max-y0)))/(g)#t value where y=y_max
    x_y_max=x0+v0x*t_y_max#x value where y=y_max
    return x,y,x_max,y_max,x_y_max,tFlight
#Examples
#g
plt.plot(projectile()[0],projectile()[1],label="Example on Earth (g=9.8)")
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
    print("\\theta = ",angle)
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
