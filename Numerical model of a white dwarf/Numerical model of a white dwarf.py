# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 12:07:05 2023
Numerical model of a white dwarf (based on the pdf file attached with this code) to approximate the already known behavior and structure by searching for various parameters important for its study. Specially, the relationship between
the star's radius R and the star's mass M with respect to its central density rho_c and the consequent relationship between M and R will be analyzed. Moreover, the pressure profiles P and density rho with respect to R will be
investigated as well, and these numerical results will be compared with those of a polytropic (solutions of the Lane-Emden) model (with n = 3) and contrasted with relativistic and non-relativistic theoretical results. All units of 
variables are in SI except some marked of length and mass, which are in cm and g, respectively (for convenience purposes).
More information on these websites: https://en.wikipedia.org/wiki/White_dwarf
                                    https://www.astro.princeton.edu/~burrows/classes/403/white.dwarfs.pdf
                                    https://en.wikipedia.org/wiki/Polytrope
                                    https://en.wikipedia.org/wiki/Lane%E2%80%93Emden_equation
                                    https://en.wikipedia.org/wiki/Newton%27s_method
                                    https://en.wikipedia.org/wiki/Midpoint_method
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants, initial conditions (ICs), arrays, lists, and functions
M_sun=2*10**(33)#[g]
R_sun=6.96*10**(10)#[cm]
G=6.6743015*10**(-8)#Gravitational constant in [(cm^3)/(g*s^2)] units
mu_e=2#For the case of a highly electronically degenerate white dwarf
M_CH=1.44*M_sun#Chandrasekhar limit mass in [g] units
r_0=0#Initial radius in [cm] units
M_r_0=0#Mass in r_0 in [g] units
R_list=[]
M_list=[]
rho_c_values_list=[10**5,10**6,10**7,10**8,10**9,10**10,10**11]#[(g)/(cm^3)]
log_rho_c_values_list=np.log10(rho_c_values_list)
dr_values_list=[300*10**5,200*10**5,150*10**5,80*10**5,50*10**5,25*10**5,12*10**5]#The radial width of each concentric spherical layer of the star used for the computations in [cm] units
R_polytrope=10**(7)#[m] (corresponding to a white dwarf of R = 10000 [km])
M_polytrope=2*10**(30)#[kg] (corresponding to a white dwarf of M = 1 M_sun)
xi_polytrope=6.89685#Dimensionless radius related to polytropes
lambda_polytrope=R_polytrope/xi_polytrope#Parameter related to the Lane-Emden equation in [m] units
u_0=[1,0]#Initial conditions for the solution of the Lane-Emden equation
x_polytrope_list=np.linspace(10**(-6),7,100)#Dimensionless parameter important for solving the equations that build the numerical model
R_non_relativistic_list=[]
R_relativistic_list=[]
def rho(x):
    return 9.82*10**(5)*mu_e*x**3
def x_from_rho(rho):
    return ((rho)/(9.82*10**(5)*mu_e))**(1/3)
def f(x):#Function important for the solution of the Lane-Emden equation
    return x*(2*(x**(2))-3)*((1+x**(2))**(1/2))+3*np.log(x+((1+x**(2))**(1/2)))
def P(x):
    return 6.01*10**(22)*f(x)
def K(P):#Constant of proportionality (or f(x) in the document)
    return ((P)/(6.01*10**(22)))
def df(x):
    return ((8*x**(4))/((1+x**(2))**(1/2)))
def P_1(P_c,rho_c,dr):#Pressure value at the first step of the method
    return P_c-((2/3)*G*np.pi*rho_c**(2)*dr**(2))
def M_r_1(M_r_0,rho_c,dr):#M value (at radial coordinate r inside the star) at the first step (or first value of r)
    return M_r_0+((4/3)*np.pi*rho_c*dr**(3))
def x_Newton(x_j,P,prec):#Newton's method used to find the x in the whole white dwarf
    def x_j2(x_j,P):
        return x_j-((f(x_j)-K(P))/(df(x_j)))
    x_j2_list=[x_j,x_j2(x_j,P)]
    while True:
        if np.abs(x_j2_list[-1]-x_j2_list[-2])<prec:
            break
        else:
            x_j2_list.append(x_j2(x_j2_list[-1],P))
    return x_j2_list[-1]
#Midpoint method functions
#First medium step
def r_midpoint(r_i,dr):
    return r_i+0.5*dr
def P_midpoint(P_i,dr,rho_i,M_r_i,r_i):
    return (P_i-((0.5*dr*rho_i*G*M_r_i)/(r_i**(2))))
def M_r_midpoint(M_r_i,dr,rho_i,r_i):
    return M_r_i+0.5*dr*4*np.pi*rho_i*r_i**(2)
def x_midpoint(x_j,P_midpoint,prec):
    return x_Newton(x_j,P_midpoint,prec)
def rho_midpoint(x_midpoint):
    return 1.964*10**(6)*x_midpoint**(3)
#Second medium step
def r_midpoint2(r_i,dr):
    return r_i+dr
def P_midpoint2(P_i,dr,rho_midpoint,M_r_midpoint,r_midpoint):
    return (P_i-((dr*rho_midpoint*G*M_r_midpoint)/(r_midpoint**(2))))
def M_r_midpoint2(M_r_i,dr,rho_midpoint,r_midpoint):
    return M_r_i+dr*4*np.pi*rho_midpoint*r_midpoint**(2)
def x_midpoint2(x_midpoint,P_midpoint2,prec):
    return x_Newton(x_midpoint,P_midpoint2,prec)
def rho_midpoint2(x_midpoint2):
    return 1.964*10**(6)*x_midpoint2**(3)
#Main algorithm of the numerical method
def main_program(rho_c,dr,save,prec):#"prec" denotes the imposed accuracy of the method
    P_c=P(x_from_rho(rho_c))#Central pressure value
    r_1=r_midpoint2(r_0,dr)#r value at the first step of the method
    x_0=x_from_rho(rho_c)#Initial x
    x_1=x_Newton(x_0,P_c,prec)#x value at the first step of the method
    rho_1=rho(x_1)#rho value at the first step of the method
    i=0#Iteration number of the numerical method
    P_list=[P_c,P_1(P_c,rho_c,dr)]
    M_r_list=[M_r_0,M_r_1(M_r_0,rho_c,dr)]
    r_list=[r_0,r_1]
    x_list=[x_0,x_1]
    rho_list=[rho_c,rho_1]
    i_list=[i]
    while True:
        i+=1
        i_list.append(i)
        if P_midpoint(P_list[-1],dr,rho_list[-1],M_r_list[-1],r_list[-1])<0:
            break
        r_int=r_midpoint(r_list[-1],dr)
        P_int=P_midpoint(P_list[-1],dr,rho_list[-1],M_r_list[-1],r_list[-1])
        M_r_int=M_r_midpoint(M_r_list[-1],dr,rho_list[-1],r_list[-1])
        x_int=x_midpoint(x_list[-1],P_int,prec)
        rho_int=rho_midpoint(x_int)
        r2=r_midpoint2(r_list[-1],dr)
        P2=P_midpoint2(P_list[-1],dr,rho_int,M_r_int,r_int)
        M_r2=M_r_midpoint2(M_r_list[-1],dr,rho_int,r_int)
        x2=x_midpoint2(x_int,P2,prec)
        rho2=rho_midpoint2(x2)
        r_list.append(r2)
        P_list.append(P2)
        M_r_list.append(M_r2)
        x_list.append(x2)
        rho_list.append(rho2)
    P_array=np.zeros(len(P_list))
    M_r_array=np.zeros(len(M_r_list))
    r_array=np.zeros(len(r_list))
    rho_array=np.zeros(len(rho_list))
    for i in range(0,len(r_list)):#To obtain normalized and adimensional values for the graphs
        P_array[i]=np.log10(P_list[i])
        M_r_array[i]=M_r_list[i]/M_sun
        r_array[i]=r_list[i]/R_sun
        rho_array[i]=np.log10(rho_list[i])
    if save==True:#If the data obtained has to be saved
        data_NT_array=np.array((i_list,r_array,M_r_array,P_array,rho_array,x_list))#"NT": Not Transposed
        data_T_array=data_NT_array.T#"T": Transposed
        variables_list=["i","r","M_r","log (P)","log (rho)","x"]
        with open("data_white_dwarf_numerical_method.csv","w") as file:
            a=csv.writer(file)
            a.writerow(variables_list)
            a.writerows(data_T_array)
    return (r_array,M_r_array,P_array,rho_array,x_list)
#Results + plots of R and M as a function of rho_c and M as a function of R (Numerical method)
j=0
for rho_c in rho_c_values_list:
    R_list.append(main_program(rho_c=rho_c,dr=dr_values_list[j],save=False,prec=10**(-6))[0][-1])
    M_list.append(main_program(rho_c=rho_c,dr=dr_values_list[j],save=False,prec=10**(-6))[1][-1])
    j+=1
plt.figure(figsize=(9,5))
plt.plot(log_rho_c_values_list,R_list,marker="d",c="red")
plt.xlabel(r"$log(\rho_{c})$")
plt.ylabel(r"$R/R_{\odot}$")
plt.title(r"Dependence of $R$ as a function of $\rho_{c}$")
plt.grid()
plt.tight_layout()
plt.show()
plt.figure(figsize=(9,5))
plt.plot(log_rho_c_values_list,M_list,marker="d",c="red")
plt.axhline(M_CH/M_sun,color="black",linestyle="dashed",label=r"Chandrasekhar limit mass ($M_{CH}\approx1.44M_{\odot}$)")
plt.xlabel(r"$log(\rho_{c})$")
plt.ylabel(r"$M/M_{\odot}$")
plt.title(r"Dependence of $M$ as a function of $\rho_{c}$")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#Results + plots of the comparison of rho and P between the numerical and polytropic (with n = 3) methods
rho_aver_value=M_polytrope/((4/3)*np.pi*(R_polytrope)**(3))
rho_c=rho_aver_value*(1/3)*((xi_polytrope**(3))/(2.01))
K_polytrope=lambda_polytrope**(2)*np.pi*G*rho_c**(2/3)
def f_polytrope(u,x):
    return (u[1],-2/x*u[1]-u[0]**(3))
u_polytrope=odeint(f_polytrope,u_0,x_polytrope_list)
u_list=u_polytrope[:,0]
r_polytrope_list=x_polytrope_list*(lambda_polytrope/R_sun)*10**(2)#Factor de 100 añadido para que concuerden las magnitudes de las variables polítropicas con las numéricas.
rho_polytrope_list=rho_c*u_list**(3)
P_polytrope_list=K_polytrope*rho_polytrope_list**(4/3)
num_sol_list=main_program(rho_c,dr=10**(5),save=False,prec=10**(-6))
rho_num_list=np.zeros(len(num_sol_list[3]))
P_num_list=np.zeros(len(num_sol_list[2]))
r_num_list=num_sol_list[0]*10#Factor of 10 added to make the magnitudes of the polytropic variables match the numerical ones.
for j in range(0,len(rho_num_list)):
    rho_num_list[j]=10**(num_sol_list[3][j])
    P_num_list[j]=(10**(num_sol_list[2][j]))/(100)#Factor of 100 added to make the magnitudes of the polytropic variables match the numerical ones.
plt.figure(figsize=(9,5))
plt.plot(r_polytrope_list,rho_polytrope_list,c="red",label=r"$n=3$ polytrope")
plt.plot(r_num_list,rho_num_list,c="green",label="Numerical")
plt.xlabel(r"$r/R_{\odot}$ []")
plt.ylabel(r"$\rho [kg/m^{3}]$")
plt.title(r"Comparison between methods of $\rho$")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
plt.figure(figsize=(9,5))
plt.plot(r_polytrope_list,P_polytrope_list,c="red",label=r"$n=3$ polytrope")
plt.plot(r_num_list,P_num_list,c="green",label="Numerical")
plt.xlabel(r"$r/R_{\odot}$ []")
plt.ylabel(r"$P [(kg\cdot m)/(s^{2})]$")
plt.title(r"Comparison between methods of $P$")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#Results + plots of the comparison of M as a function of R between the numerical and theoretical (relativistic and non-relativistic) methods
def R_non_relativistic(M):
    return ((R_sun)/(74))*((M_sun)/(M))**(1/3)
def R_relativistic(M):
    return 0.0126*R_sun*(((M_sun)/(M))**(1/3))*((1-((M)/(M_CH))**(4/3))**(1/2))
for h in range(0,len(M_list)):
    R_non_relativistic_list.append(R_non_relativistic(M_sun*M_list[h])/(R_sun))#To obtain normalized and adimensional values for the graphs
    R_relativistic_list.append(R_relativistic(M_sun*M_list[h])/(R_sun))#To obtain normalized and adimensional values for the graphs
plt.figure(figsize=(9,5))
plt.plot(R_non_relativistic_list,M_list,marker="d",c="red",label="Theoretical (Non-relativistic)")
plt.plot(R_relativistic_list,M_list,marker="d",c="green",label="Theoretical (Relativistic)")
plt.plot(R_list,M_list,marker="d",c="blue",label="Numerical")
plt.xlabel(r"$R/R_{\odot}$")
plt.ylabel(r"$M/M_{\odot}$")
plt.title(r"Comparison between methods of $M$ as a function of $R$")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
print("Program execution time:",time.process_time()-start_time_program,"seconds.")