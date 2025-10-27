# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 21:39:59 2022
Numerical resolution of the Thomas-Fermi(TF) model/equation: d2φ(x)/dx2 = φ''(x) = (φ(x)^(3/2))/(x^(1/2)), with adimensional variables φ(x) (universal function) and x (coordinate), for free neutral atoms with the imposed Initial and Boundary Conditions (ICs and BCs): φ(x = x0 = 0) = 1, φ'(x = x0 = 0) = –1.5880710226, φ(x = xf -> ∞) = 0, from x = 0 outwards with number of electrons N = [Ni,Nf]. With the assistance of the pdf file attached with this code. All units of variables are in atomic units [au].
More information on these websites: https://en.wikipedia.org/wiki/Thomas%E2%80%93Fermi_model
                                    https://en.wikipedia.org/wiki/Atomic_units
@author: Lucas Romero Fernández
"""
import time
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants, arrays, lists, and functions
Ni=10
N_change=50#Where, approximately, the behavior of the energies changes
Nf=200
x0=np.exp(-200)#e^(-200) ≈ 0
x_change=5#Where, approximately, the behavior of the universal function and its derivative change
xf=35#φ(x = 35) ≈ φ(x -> ∞), φ'(x = 35) ≈ φ'(x -> ∞)
x_array1=np.linspace(x0,x_change,75)
x_array2=np.linspace(x_change+0.5,xf,100)
x_array=np.concatenate((x_array1,x_array2))
N_array=np.arange(Ni,Nf+1)#N = Z (for neutral atoms)
ICs_list=[1,-1.5880710226]#Initial conditions for φ(x) and φ'(x), respectively
x_tabulated_list=[0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.20,1.40,1.60,1.80,2.00,2.20,2.40,2.60,2.80,3.00,3.50,4.00,4.50,5.00,5.50,6.00,7.00,8.00,9.00,10.00,12.00,14.00,16.00,18.00,20.00,25.00,30.00,35.00]#From the pdf file attached to this code
phi_tabulated_list=[1.000000000,0.971976639,0.946962774,0.923826768,0.902153033,0.881697077,0.862291781,0.843813275,0.826164915,0.809268576,0.793059432,0.720639476,0.659541161,0.606986383,0.561162024,0.520791457,0.484930988,0.452858715,0.424008052,0.374241230,0.332901370,0.298097707,0.268469510,0.243008507,0.220949979,0.201702701,0.184802149,0.169878264,0.156632673,0.129369597,0.108404257,0.091948134,0.078807779,0.068160362,0.059422949,0.046097819,0.036587255,0.029590936,0.024314293,0.017063922,0.012478407,0.009424081,0.007304845,0.005784941,0.003473753,0.002255835,0.001550930]#From the pdf file attached to this code
phi_prime_tabulated_list=[1.588071023,1.309304963,1.199103451,1.117740319,1.051608480,0.995354646,0.946194872,0.902454221,0.863028591,0.827142829,0.794227009,0.661799780,0.564642444,0.489411613,0.429171872,0.379794745,0.338607156,0.303775756,0.273989052,0.225908594,0.189041426,0.160115008,0.136998438,0.118243192,0.102830976,0.090026276,0.079285763,0.070200388,0.062457131,0.047501047,0.036943758,0.029271448,0.023560075,0.019221348,0.015867549,0.011142532,0.008088603,0.006033075,0.004602882,0.002830536,0.001844501,0.001257435,0.000888831,0.000647254,0.000324043,0.000180670,0.00010891]#From the pdf file attached to this code
Total_E_list=[]
T_list=[]
Uen_list=[]
Uee_list=[]
U_list=[]
def Total_E(N):#Total energy of the electrons
    return -0.7687*N**(7/3)
def T(N):#Kinetic energy of the electrons
    return -Total_E(N)
def Uen(N):#Potential energy of the electrons due to the electric attraction of the positively charged nucleus
    return -1.79*N**(7/3)
def Uee(N):#Potential energy of the electrons due to their mutual electric repulsion
    return -(1/7)*Uen(N)
def U(N):#Total potential energy of the electrons
    return Uen(N)+Uee(N)
def Equation(phi,x):#Equation of the model reduced to a system of ordinary differential equations
    return (phi[1],(phi[0]**(3/2))/(x**(1/2)))
#Main algorithm of the resolution of the equation
phi_int_array=odeint(Equation,ICs_list,x_array)
phi_array=phi_int_array[:,0]
phi_prime_array=-phi_int_array[:,1]
for i in N_array:
    Total_E_list.append(Total_E(i))
    T_list.append(T(i))
    Uen_list.append(Uen(i))
    Uee_list.append(Uee(i))
    U_list.append(U(i))
#Graphs
plt.rcParams["figure.figsize"]=(9,5)
#Comparison between universal functions of Thomas-Fermi
#For small x values
plt.xlim(x0,x_change)
plt.plot(x_array,phi_array,color='red',label="Numerical")
plt.plot(x_tabulated_list,phi_tabulated_list,c='green',linestyle="dashed",label="Tabulated")
plt.xlabel("$x$[ua]")
plt.ylabel("$\phi(x)$[ua]")
plt.title("Comparison between results of $\phi(x)$ for small $x$ values")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#General
plt.xlim(x0-1,xf+1)
plt.plot(x_array,phi_array,color='red',label="Numerical")
plt.plot(x_tabulated_list,phi_tabulated_list,c='green',linestyle="dashed",label="Tabulated")
plt.xlabel("$x$[ua]")
plt.ylabel("$\phi(x)$[ua]")
plt.title("Comparison between results of $\phi(x)$")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#Comparison between the (minus) derivatives of the universal functions of Thomas-Fermi
#For small x values
plt.xlim(x0-0.1,x_change+0.1)
plt.plot(x_array,phi_prime_array,color='red',label="Numerical")
plt.plot(x_tabulated_list,phi_prime_tabulated_list,c='green',linestyle="dashed",label="Tabulated")
plt.xlabel("$x$[ua]")
plt.ylabel("$-\phi′(x)$[ua]")
plt.title("Comparison between results of $-\phi'(x)$ for small $x$ values")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#General
plt.xlim(x0-1,xf+1)
plt.plot(x_array,phi_prime_array,color='red',label="Numerical")
plt.plot(x_tabulated_list,phi_prime_tabulated_list,c='green',linestyle="dashed",label="Tabulated")
plt.xlabel("$x$[ua]")
plt.ylabel("$-\phi′(x)$[ua]")
plt.title("Comparison between results of $-\phi'(x)$")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#Evolution of the different energies
#With respect to small N = Z values
plt.xlim(Ni-0.5,N_change+0.5)
plt.ylim(-20000,10000)
plt.plot(N_array,Total_E_list,color='red',label="$E_{Total}$")
plt.plot(N_array,T_list,color='green',label="$T$")
plt.plot(N_array,Uen_list,color='blue',label="$U_{en}$")
plt.plot(N_array,Uee_list,color='orange',label="$U_{ee}$")
plt.plot(N_array,U_list,color='purple',label="$U$")
plt.xlabel("$N$")
plt.ylabel("Energies [ua]")
plt.title("Energy evolution with respect to small $N$ values")
plt.legend(loc=(0.05,0.05))
plt.grid()
plt.tight_layout()
plt.show()
#With respect to all N = Z values considered
plt.xlim(Ni-0.5,Nf+0.5)
plt.plot(N_array,Total_E_list,color='red',label="$E_{Total}$")
plt.plot(N_array,T_list,color='green',label="$T$")
plt.plot(N_array,Uen_list,color='blue',label="$U_{en}$")
plt.plot(N_array,Uee_list,color='orange',label="$U_{ee}$")
plt.plot(N_array,U_list,color='purple',label="$U$")
plt.xlabel("$N$")
plt.ylabel("Energies [ua]")
plt.title("Energy evolution with respect to $N$")
plt.legend(loc=(0.05,0.05))
plt.grid()
plt.tight_layout()
plt.show()
print("Program execution time:",time.process_time()-start_time_program,"seconds.")