# Wave equation in 3d: ParallelStencil version
# author: Lucas Romero Fernández 2025
# Multithreaded option:
# Start by using the desired number of threads "t" with "julia -t 4 Solver_for_3D_wave_equation_ps.jl"
# or by using the environment var export JULIA_NUM_THREADS=4 (in both cases, 4 threads selected)
# GPU option:
# Start with "julia Solver_for_3D_wave_equation_ps.jl gpu"

using Plots
using ParallelStencil

@static if "gpu" in ARGS
    @init_parallel_stencil(CUDA,Float64,3)
else
    @init_parallel_stencil(Threads,Float64,3)
end

Gaussian(x,y,z,alpha,dt,μ=0.5,σ=0.1) = exp(-((x + alpha*dt)-μ)^2/2σ^2 -((y + alpha*dt)-μ)^2/2σ^2 -((z + alpha*dt)-μ)^2/2σ^2)

function Initialcond!(nx,ny,nz,u₀,u₁,x,y,z,alpha,dt)
    u₀ .= Data.Array([Gaussian(x[i],y[j],z[k],alpha,-dt) for i=1:nx,j=1:ny,k=1:nz])
    u₁ .= Data.Array([Gaussian(x[i],y[j],z[k],alpha,0) for i=1:nx,j=1:ny,k=1:nz])
    return
end

function Introduction()
    println("Wave Equation 3D using ParallelStencil") 
    if "gpu" in ARGS println("and using CUDA GPU") end
end

@parallel_indices (i,j,k) function ExpEulerSecOrd!(tau,u,u₀,u₁)
    # 2nd order difference quotient (time) and discrete Laplace operator of 2th order (spatial) or, equivalently, an explicit Euler method of second order
    # From: https://beltoforion.de/en/recreational_mathematics/3d-wave-equation.php
    u[i,j,k] = tau*(u₁[i-1,j,k] + u₁[i+1,j,k] + u₁[i,j-1,k] + u₁[i,j+1,k]
    + u₁[i,j,k-1] + u₁[i,j,k+1] - 6*u₁[i,j,k]) + 2*u₁[i,j,k] - u₀[i,j,k]
    return
end

function Waveevo!(tau,nx,ny,nz,u,u₀,u₁)
    @parallel (2:nx-1,2:ny-1,2:nz-1) ExpEulerSecOrd!(tau,u,u₀,u₁)
    return
end

@parallel_indices (j,k) function RadiativeBoundcond!_x!(kappa,A::Data.Array,A₁::Data.Array)
    # Radiative/absorbing boundary conditions for x
    A[1,j,k] = A₁[2,j,k] + ((kappa - 1)/(kappa + 1))*(A[2,j,k] - A₁[1,j,k])
    A[end,j,k] = A₁[end-1,j,k] + ((kappa - 1)/(kappa + 1))*(A[end-1,j,k] - A₁[end,j,k])
    return
end

@parallel_indices (i,k) function RadiativeBoundcond!_y!(kappa,A::Data.Array,A₁::Data.Array)
    # Radiative/absorbing boundary conditions for y
    A[i,1,k] = A₁[i,2,k] + ((kappa - 1)/(kappa + 1))*(A[i,2,k] - A₁[i,1,k])
    A[i,end,k] = A₁[i,end-1,k] + ((kappa - 1)/(kappa + 1))*(A[i,end-1,k] - A₁[i,end,k])
    return
end

@parallel_indices (i,j) function RadiativeBoundcond!_z!(kappa,A::Data.Array,A₁::Data.Array)
    # Radiative/absorbing boundary conditions for z
    A[i,j,1] = A₁[i,j,2] + ((kappa - 1)/(kappa + 1))*(A[i,j,2] - A₁[i,j,1])
    A[i,j,end] = A₁[i,j,end-1] + ((kappa - 1)/(kappa + 1))*(A[i,j,end-1] - A₁[i,j,end])
    return
end

function Boundcond!(kappa,nx,ny,nz,u,u₁)
    @parallel (1:ny-1,1:nz-1) RadiativeBoundcond!_x!(kappa,u,u₁)
    @parallel (1:nx-1,1:nz-1) RadiativeBoundcond!_y!(kappa,u,u₁)
    @parallel (1:nx-1,1:ny-1) RadiativeBoundcond!_z!(kappa,u,u₁)
    return
end

function Simulation(N,nplot=20,tend=1,cfl=0.5,alpha=1)
    # Initial grid, deltas and constants
    nx = ny = nz = N #dx=dy=dz
    dx,dy,dz = 1/(nx-1), 1/(ny-1), 1/(nz-1)
    x,y,z = 0:dx:1, 0:dy:1, 0:dz:1
    dt = cfl*dx # Courant condition
    tau = ((alpha*dt)/dx)^2
    kappa = dt*alpha/dx
    
    # Start code
    Introduction()
    
    # Initial conditions
    t = 0.
    iter = 0
    u = @zeros(nx,ny,nz)
    u₀ = @zeros(nx,ny,nz)
    u₁ = @zeros(nx,ny,nz)
    Initialcond!(nx,ny,nz,u₀,u₁,x,y,z,alpha,dt)

    # viz (movie animation)
    domovie = nplot > 0
    if domovie
        gr(); ENV["GKSwstype"]="nul"        
        anim = Animation()
        uv = Array(u₀)
        println("Movie frame at iteration $iter"," umax=",maximum(u₁))
        surface(x,y,uv[:,:,nz÷2],title="Time $t") # To plot the initial condition
        frame(anim)
    end

    walltime = 0
    # Main evolution loop
    while t <= tend 
        iter += 1
        t += dt
        t0 = time()
        Waveevo!(tau,nx,ny,nz,u,u₀,u₁)
        Boundcond!(kappa,nx,ny,nz,u,u₁)
        walltime += time() - t0
        if domovie && iter%nplot==0 
            println("Movie frame at iteration $iter"," umax=",maximum(u))
            uv = Array(u)
            surface(x,y,uv[:,:,nz÷2],title="Time "*string(round(t,digits=2)))
            frame(anim) 
        end
        # Update the time positions of the solutions for new iterations
        u₀, u₁ = u₁, u
    end
    println("Computation time: ",walltime)
    if domovie mp4(anim,"wave3d_$nx.mp4") end
end

# main_program
Simulation(401,20)