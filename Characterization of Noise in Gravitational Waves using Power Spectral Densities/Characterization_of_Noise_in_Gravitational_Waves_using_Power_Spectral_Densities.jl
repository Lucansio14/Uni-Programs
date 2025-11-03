# Characterization of Noise in Gravitational Waves (GWs) using Power Spectral Densities (PSDs) 
# Computed in Julia 1.11.3
# author: Lucas Romero Fernández 2025
# Start with "julia Characterization_of_Noise_in_Gravitational_Waves_using_Power_Spectral_Densities.jl"

using DSP # For the periodograms
using FFTW # For the Discrete Fourier Transforms (DFTs)
using LaTeXStrings # For the plots of the results
using PyPlot # For the plots of the results
using Random # For the generation of random noise.

# Definition of general constants, variables, arrays and functions

T_coh = 300 # Coherent time (in this case, equivalent to the observation time of 5 minutes) (Units: s)
f_Nyq = 2000 # Nyquist frequency (Units: Hz)
f_i = 0 # Initial frequency (Units: Hz)
f_f = f_Nyq # End frequency (Units: Hz)
f_s = 2*(f_f - f_i) # Sampling frequency (Units: Hz)
Deltat = 1/f_s # Time step size (Units: s)
N = T_coh*f_s # Number of sample points
t_array = collect(0:N-1)*Deltat
Random.seed!(1418) # Random noise seed setting (to always have the same random result)
f_c = 10 # Lower cut-off frequency (Units: Hz)
function Sn(f) # Theoretical noise/sensitivity curve (PSD) for Advanced LIGO (Units: Hz^(-1))
    f_0 = 215 # Scaling factor (Units: Hz)
    S_0 = 10^(-49) # "Amplitude" (Units: Hz^(-1))
    x = f./f_0
    return ifelse.(f .>= f_c,S_0*(x.^(-4.14) .- 5*x.^(-2) .+ 111*(1 .- x.^2 .+ x.^4/2)./(1 .+ x.^2/2)),0)
    end

# main_program

f_array = fftshift(fftfreq(N,f_s)) # Frequency array computation
f_pos_array = f_array[findall(f_array .>= 0)] # Only positive frequencies array

# Computing of the positive noise values in the frequency domain
x_k_pos_array = 0.5 .* (randn(length(f_pos_array)) .+ 1im .* randn(length(f_pos_array))
) .* sqrt.(T_coh .* Sn(f_pos_array))

# Filling the noise array with the remaining values corresponding to negative frequencies
x_k_array = [conj.(reverse(x_k_pos_array));0;x_k_pos_array[1:end-1]]

x_k_sort_array = ifftshift(x_k_array) # Sort the noise array from lowest to highest frequencies
x_j_array = ifft(x_k_sort_array)/Deltat # Transform the noise array to the time domain

# Computation of the periodograms

# Periodogram using only the FFT estimation
FFT_array = periodogram(
real(x_j_array), # Detector/noise output
fs = f_s, # Frequency sampling
window = nothing # No window method applied
)
# Welch method for different amount of intervals
# 2 intervals (M = 2)
Welch2_array = welch_pgram(
real(x_j_array), # Detector/noise output
n = div(length(real(x_j_array)),2), # Samples in each interval
fs = f_s, # Frequency sampling
window = hanning , # Hann window method
)
# 4 intervals (M = 4)
Welch4_array = welch_pgram(
real(x_j_array), # Detector/noise output
n = div(length(real(x_j_array)),4), # Samples in each interval
fs = f_s, # Frequency sampling
window = hanning , # Hann window method
)
# 10 intervals (M = 10)
Welch10_array = welch_pgram(
real(x_j_array), # Detector/noise output
n = div(length(real(x_j_array)),10), # Samples in each interval
fs = f_s, # Frequency sampling
window = hanning , # Hann window method
)
# Frequency axis for the different periodograms
f_FFT_array = range(f_i,f_f,length(FFT_array.power))
f_Welch2_array = range(f_i,f_f,length(Welch2_array.power))
f_Welch4_array = range(f_i,f_f,length(Welch4_array.power))
f_Welch10_array = range(f_i,f_f,length(Welch10_array.power))

# Plots

# Theoretical noise/sensitivity curve (PSD)
figure(figsize = (10,6))
plot(f_pos_array,log10.(Sn(f_pos_array)),c = "red")
xlim(f_c,f_f)
xlabel(L"f\ [Hz]")
ylabel(L"log(S_{n})\ [Hz^{-1}]")
grid()
title("Advanced LIGO single-sided Power Spectral Density")
tight_layout()
savefig("Theoretical noise curve")
show()

# Generated noise in the frequency domain
plt.figure(figsize = (10,6))
xlim(-f_f,f_f)
plt.plot(f_array,abs.(x_k_array),c = "cyan")
plt.xlabel(L"f\ [Hz]")
plt.ylabel(L"|n(f)|\ [Hz^{-1/2}]")
plt.grid()
plt.title("Generated noise in the frequency domain")
plt.tight_layout()
plt.savefig("Generated noise (Frequency domain)")
plt.show()

# Generated noise in the time domain
plt.figure(figsize = (10,6))
plt.xlim(0,T_coh)
plt.plot(t_array,real(x_j_array),c = "cyan")
plt.xlabel(L"t\ [s]")
plt.ylabel(L"n(t)\ [s^{-1}]")
plt.grid()
plt.title("Generated noise in the time domain")
plt.tight_layout()
plt.savefig("Generated noise (Time domain)")
plt.show()

# Comparison between PSDs
plt.figure(figsize = (10,6))
plt.xlim(f_c,f_f)
plt.ylim(-55,-40)
plt.plot(f_FFT_array,log10.(FFT_array.power),c = "blue",label = "FFT only",alpha = 0.8)
plt.plot(f_Welch2_array,log10.(Welch2_array.power),c = "cyan",label = "Welch method (M = 2)",alpha = 0.8)
plt.plot(f_Welch4_array,log10.(Welch4_array.power),c = "grey",label = "Welch method (M = 4)",alpha = 0.8)
plt.plot(f_Welch10_array,log10.(Welch10_array.power),c = "orange",label = "Welch method (M = 10)",alpha = 0.8)
plt.plot(f_pos_array,log10.(Sn(f_pos_array)),c = "red",label = L"S_{n,theo}")
plt.xlabel(L"f\ [Hz]")
plt.ylabel(L"log(S_{n})\ [Hz^{-1}]")
plt.grid()
plt.title("Comparison between Power Spectral Densities periodograms")
plt.legend()
plt.tight_layout()
plt.savefig("Comparison between PSDs")
plt.show()