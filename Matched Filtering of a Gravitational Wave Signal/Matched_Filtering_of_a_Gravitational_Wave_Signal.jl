# Matched Filtering of a Gravitational Wave (GW) Signal
# Computed in Julia 1.11.3
# author: Lucas Romero Fernández 2025
# Start with "julia Matched_Filtering_of_a_Gravitational_Wave_Signal.jl"

using QuadGK # For the primitives
using LinearAlgebra # For the vector and matrix operations
using PyPlot # For the plots of the results

# Definition of general constants, variables, arrays and functions

parsec = 3.0857e16 # Conversion from 1 parsec (pc) to meters (Units: m)
G = 6.674e-11 # Gravitational constant (Units: m^(3)kg^(-1)s^(-2))
M_sun = 1.98847e30 # Conversion from 1 solar mass to kilograms (Units: kg)
c = 299792458.0 # Speed of light in vaccuum (Units: m/s)
Gpc = 1e9*(parsec/c) # Conversion of gigaparsecs (length) to units of time (pc/c) (Units: s)
m_sun_time = M_sun*(G/(c^3)) # Conversion of solar masses (mass) to units of time (M_sun*G*(c^(-3)) (Units: s)
f_c = 10 # Lower cut-off frequency and lower limit of the SNR integral (Units: Hz)
Masses_array = range(0,stop = 3,length = 2000) .|> x->10^x # Total masses array (Units: Solar masses)
f_ISCO(M) = 1/(M*pi*6^(3/2)) # Mass dependent f_ISCO function (superior limit of the SNR integral) (Units: Hz)
function Sn(f) # Theoretical noise/sensitivity curve (PSD) for Advanced LIGO (Units: Hz^(-1))
    f_0 = 215 # Scaling factor (Units: Hz)
    S_0 = 1e-49 # "Amplitude" (Units: Hz^(-1))
    x = f./f_0
    return ifelse.(f .>= f_c,S_0*(x.^(-4.14) .- 5*x.^(-2) .+ 111*(1 .- x.^2 .+ x.^4/2)./(1 .+ x.^2/2)),Inf)
    end

# Section A: Computation and plotting of the SNR at a distance of 1 Gpc for an equal mass case

# Definition of specific functions

A_rms_INS(chirpM,r) = (chirpM^(5/6))/(r*sqrt(30)*pi^(2/3)) # Mass dependent A^rms_INS function (Units: kg^(5/6)/s)
function SNR(m_1,m_2,r) # Signal-to-Noise Ratio (SNR) function with masses m_1 and m_2 and horizon distance r
    M = (m_1 + m_2) # Total mass (Units: kg)
    mu = (m_1*m_2)/M # Reduced mass (Units: kg)
    eta = mu/M # Symmetric mass ratio (No units)
    chirpM = M^(2/5)*mu^(3/5) # Chirp mass (Units: kg)
    to_integrate(f) = (f^(-7/3))/Sn(f) # Function to integrate
    integral = quadgk(to_integrate,f_c,f_ISCO(M))[1] # Integration computation
    return sqrt(4*A_rms_INS(chirpM,r)^2*integral)
    end

# main_program

SNR_array = [SNR((m/2)*m_sun_time,(m/2)*m_sun_time,Gpc) for m in Masses_array] # r = 1 Gpc

# Plot

figure(figsize = (10,6))
plot(Masses_array,SNR_array,color = "black")
xscale("log")
yscale("log")
xlim(1,1000)
xlabel(L"Total\ mass\ (M)\ [M_\odot]")
ylabel(L"Signal-To-Noise\ Ratio\ (SNR)\ (\rho)")
grid()
title("SNR at a distance of 1 Gpc for the case of equal masses")
tight_layout()
savefig("Results_Section_A.png")
show()

# Section B: Computation and plotting of the horizon distance for an SNR = 8 and an equal mass case

# Definition of specific functions

function r(m_1,m_2,rho) # Horizon distance (r) function, 
    # obtained by rewriting the SNR (rho) function from section A
    M = (m_1+m_2) # Total mass (Units: kg)
    mu = (m_1*m_2)/M # Reduced mass (Units: kg)
    eta = mu/M # Symmetric mass ratio (No units)
    chirpM = M^(2/5)*mu^(3/5) # Chirp mass (Units: kg)
    to_integrate(f) = f^(-7/3)/Sn(f) # Function to integrate
    integral = quadgk(to_integrate,f_c,f_ISCO(M))[1] # Integration computation
    return (((2*chirpM^(5/6))/(rho*sqrt(30)*pi^(2/3)))*sqrt(integral))/Gpc
    end

# main_program

r_array = [r((m/2)*m_sun_time,(m/2)*m_sun_time,8) for m in Masses_array] # SNR = 8

# Plot

figure(figsize = (10,6))
plot(Masses_array,r_array,color = "black")
xscale("log")
yscale("log")
xlim(1,1000)
xlabel(L"Total\ mass\ (M)\ [M_\odot]")
ylabel(L"Horizon\ distance\ (r)\ [Gpc]")
grid()
title("Horizon distance with SNR of 8 for the case of equal masses")
tight_layout()
savefig("Results_Section_B.png")
show()

# Section C: Computation of the overlap and Fitting Factors of the signal with a phase at 2 PN with the one at 1.5
# PN for different mass cases (for 2 PN, all terms have been considered, for 1.5 PN, only up to A_3).
# Mass cases to be studied:
# - m_1 = m_2 = 1.4 M_sun
# - m_1 = m_2 = 10 M_sun
# - m_1 = 10 M_sun and m_2 = 1.4 M_sun

# Definition of specific arrays and functions

x_array1 = range(1.38,stop = 1.41,length = 250) # For the first case
x_array2 = range(10.,stop = 15.,length = 500) # For the second case
x_array3 = range(8.0,stop = 8.5,length = 250) # For the third case
y_array1 = range(1.38,stop = 1.41,length = 250) # For the first case
y_array2 = range(7.,stop = 9.,length = 250) # For the second case
y_array3 = range(1.5,stop = 1.75,length = 250) # For the third case
function varPhase(m_1,m_2,f,PN) # Variable phase function with masses m_1 and m_2 and frequency f, 
# between both phases considered of different Post-Newtonian (PN) orders
    M = m_1 + m_2 # Total mass (Units: kg)
    mu = (m_1*m_2)/M # Reduced mass (Units: kg)
    eta = mu/M # Symmetric mass ratio (No units)
    chirpM = M*eta^(3/5) # Chirp mass (Units: kg)
    u = (pi*chirpM*f)^(1/3) # u term
    A_0 = 1 # Coefficient of u of order 0
    A_1 = 0 # Coefficient of u of order 1
    A_2 = 20/9*(743/336 + (11/4)*eta)*eta^(-2/5) # Coefficient of u of order 2
    A_3 = -16*pi*eta^(-3/5) # Coefficient of u of order 3
    A_4 = 10*(3058673/1016064 + (5429/1008)*eta + (617/144)*eta^2)*eta^(-4/5) # Coefficient of u of order 4
    if PN == 2 # Variable phase for PN = 2 case
        return (3/128)*(A_0/u^5 + A_1/u^4 + A_2/u^3 + A_3/u^2 + A_4/u)
    elseif PN == 1.5 # Variable phase for PN = 1.5 case
        return (3/128)*(A_0/u^5 + A_1/u^4 + A_2/u^3 + A_3/u^2)
    else
        println("Variable phase not specified for this PN")
    end
    end
function overlap(m_1,m_2) # Overlap function with masses m_1 and m_2 with the specifics of this section
    m_1 = m_1*m_sun_time # Conversion of mass units to time units (Units: s)
    m_2 = m_2*m_sun_time # Conversion of mass units to time units (Units: s)
    f_ISCO_sp = f_ISCO(m_1+m_2) # f_ISCO function computation for the specific case of this section (Units: Hz)
    # Numerator computation
    to_integrate_top(f) = (cos(varPhase(m_1,m_2,f,2) - varPhase(m_1,m_2,f,1.5)))*(f^(-7/3)/Sn(f)) # I.Function 
    numerator = quadgk(to_integrate_top,f_c,min(f_ISCO_sp,f_ISCO_sp))[1] # Integration computation
    # Denominator computation
    to_integrate_bottom(f) = f^(-7/3)/Sn(f) # Function to integrate
    integral_bottom1 = quadgk(to_integrate_bottom,f_c,f_ISCO_sp)[1] # Integrations computations
    integral_bottom2 = quadgk(to_integrate_bottom,f_c,f_ISCO_sp)[1] # Integrations computations
    denominator = sqrt(integral_bottom1*integral_bottom2)
    return numerator/denominator
    end
function meshgrid(x_in_array,y_in_array) # Mesh grid function for the 2D sampling, 
# with array x_in_array and array y_in_array
    n_x = length(x_in_array) # Number of points in direction x
    n_y = length(y_in_array) # Number of points in direction y
    x_out_array = zeros(n_y,n_x) # Construction of resulting array in direction x
    y_out_array = zeros(n_y,n_x) # Construction of resulting array in direction y
    for j_x = 1:n_x
        for i_x = 1:n_y
            x_out_array[i_x,j_x] = x_in_array[j_x] # Mesh grid
            y_out_array[i_x,j_x] = y_in_array[i_x] # Mesh grid
            end
        end
    return (x_array = x_out_array,y_array = y_out_array)
    end
function overlap_plots(m_1b,m_2b,m_1a,m_2a) # Overlaps for the plots function, 
# with masses m_1 and m_2 and 2D sampling for both PNs.
    m_1a = m_1a*m_sun_time # Conversion of mass units to time units (Units: s)
    m_2a = m_2a*m_sun_time # Conversion of mass units to time units (Units: s)
    m_1b = m_1b*m_sun_time # Conversion of mass units to time units (Units: s)
    m_2b = m_2b*m_sun_time # Conversion of mass units to time units (Units: s)
    f_ISCO1 = f_ISCO(m_1a + m_2a) # f_ISCO function computation for the specific cases of this section (Units: Hz)
    f_ISCO2 = f_ISCO(m_1b + m_2b) # f_ISCO function computation for the specific cases of this section (Units: Hz)
    # Numerator computation
    to_integrate_top(f) = (cos(varPhase(m_1a,m_2a,f,2) - varPhase(m_1b,m_2b,f,1.5)))*(f^(-7/3)/Sn(f)) # I.Function
    numerator = quadgk(to_integrate_top,f_c,min(f_ISCO1,f_ISCO2))[1] # Integration computation
    # Denominator computation
    to_integrate_bottom(f) = f^(-7/3)/Sn(f) # Function to integrate
    integral_bottom1 = quadgk(to_integrate_bottom,f_c,f_ISCO1)[1] # Integrations computations
    integral_bottom2 = quadgk(to_integrate_bottom,f_c,f_ISCO2)[1] # Integrations computations
    denominator = sqrt(integral_bottom1*integral_bottom2)
    return numerator/denominator
    end

# main_program

# Computation of the overlap for the three cases already indicated
println("Overlap for the m_1 = m_2 = 1.4 M_sun case: ",overlap(1.4,1.4))
println("Overlap for the m_1 = m_2 = 10 M_sun case: ",overlap(10,10))
println("Overlap for the m_1 = 10 M_sun and m_2 = 1.4 M_sun case: ",overlap(10,1.4))
# Final nested sampling for the different cases overlaps
X_array1,Y_array1 = meshgrid(x_array1,y_array1) # For the first case
overlap_array1 = overlap_plots.(X_array1,Y_array1,1.4,1.4) # For the first case
X_array2, Y_array2 = meshgrid(x_array2,y_array2) # For the second case
overlap_array2 = overlap_plots.(X_array2,Y_array2,10,10) # For the second case
X_array3, Y_array3 = meshgrid(x_array3,y_array3) # For the third case
overlap_array3 = overlap_plots.(X_array3,Y_array3,10,1.4) # For the third case

# Plots

# Overlap results for the first case
figure(figsize = (5,7))
pcolormesh(x_array1,y_array1,overlap_array1,cmap =:plasma)
xlabel(L"m^{PN1.5}_1\ [M_\odot]")
ylabel(L"m^{PN1.5}_2\ [M_\odot]")
xlim(minimum(x_array1),maximum(x_array1))
ylim(minimum(y_array1),maximum(y_array1))
grid()
colorbar()
title(L"Overlap\ for\ case\ m^{PN2}_1 = m^{PN2}_2 = 1.4 M_\odot")
tight_layout()
savefig("Results_Section_C_1.png")
show()
# Overlap results for the second case
figure(figsize = (5,7))
pcolormesh(x_array2,y_array2,overlap_array2,cmap =:plasma)
xlabel(L"m^{PN1.5}_1\ [M_\odot]")
ylabel(L"m^{PN1.5}_2\ [M_\odot]")
xlim(minimum(x_array2),maximum(x_array2))
ylim(minimum(y_array2),maximum(y_array2))
grid()
colorbar()
title(L"Overlap\ for\ case\ m^{PN2}_1 = m^{PN2}_2 = 10 M_\odot")
tight_layout()
savefig("Results_Section_C_2.png")
show()
# Overlap results for the third case
figure(figsize = (5,7))
pcolormesh(x_array3,y_array3,overlap_array3,cmap =:plasma)
xlabel(L"m^{PN1.5}_1\ [M_\odot]")
ylabel(L"m^{PN1.5}_2\ [M_\odot]")
xlim(minimum(x_array3),maximum(x_array3))
ylim(minimum(y_array3),maximum(y_array3))
grid()
colorbar()
title(L"Overlap\ for\ case\ m^{PN2}_1 = 10 M_\odot,\ m^{PN2}_2 = 1.4 M_\odot")
tight_layout()
savefig("Results_Section_C_3.png")
show()

# Computation of the searched Fitting Factors as the maximum values in the graphs

# FFs results for the first case
FF1 = maximum(overlap_array1) # Maximum value in the graph
index_max1 = findfirst(x_array1->x_array1 == maximum(overlap_array1),overlap_array1) # Index of the grid with FF1
m_1_index1 = x_array1[index_max1[2]] # Mass m_1 for 1.5 PN with that index
m_2_index1 = y_array1[index_max1[1]] # Mass m_2 for 1.5 PN with that index
println("For the 2 PN, m_1 = m_2 = 1.4 M_sun case: ")
println("Fitting Factor (FF): ",FF1," with masses for 1.5 PN m_1 = ",m_1_index1,"; m_2 = ",m_2_index1)
# FFs results for the second case
FF2 = maximum(overlap_array2) # Maximum value in the graph
index_max2 = findfirst(x_array2->x_array2 == maximum(overlap_array2),overlap_array2) # Index of the grid with FF2
m_1_index2 = x_array2[index_max2[2]] # Mass m_1 for 1.5 PN with that index
m_2_index2 = y_array2[index_max2[1]] # Mass m_2 for 1.5 PN with that index
println("For the 2 PN, m_1 = m_2 = 10 M_sun case: ")
println("Fitting Factor (FF): ",FF2," with masses for 1.5 PN m_1 = ",m_1_index2,"; m_2 = ",m_2_index2)
# FFs result for the third case
FF3 = maximum(overlap_array3) # Maximum value in the graph
index_max3 = findfirst(x_array3->x_array3 == maximum(overlap_array3),overlap_array3) # Index of the grid with FF3
m_1_index3 = x_array3[index_max3[2]] # Mass m_1 for 1.5 PN with that index
m_2_index3 = y_array3[index_max3[1]] # Mass m_2 for 1.5 PN with that index
println("For the 2 PN, m_1 = 10 M_sun and m_2 = 1.4 M_sun case: ")
println("Fitting Factor (FF): ",FF3," with masses for 1.5 PN m_1 = ",m_1_index3,"; m_2 = ",m_2_index3)