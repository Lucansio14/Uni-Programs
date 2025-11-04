# Wave equation in 3d: Multi-threaded version
# author: Lucas Romero Fernández 2025
# Start by using the desired number of threads "t" with "julia -t 4 Solver_for_3D_wave_equation_t.jl"
# or by using the environment var export JULIA_NUM_THREADS=4 (in both cases, 4 threads selected)

using Plots

Gaussian(x,y,z,alpha,dt,μ=0.5,σ=0.1) = exp(-((x + alpha*dt)-μ)^2/2σ^2 -((y + alpha*dt)-μ)^2/2σ^2 -((z + alpha*dt)-μ)^2/2σ^2)

function Initialcond!(nx,ny,nz,u₀,u₁,x,y,z,alpha,dt)
    for k=1:nz
        for j=1:ny
            for i=1:nx
                u₀[i,j,k] = Gaussian(x[i],y[j],z[k],alpha,-dt)
                u₁[i,j,k] = Gaussian(x[i],y[j],z[k],alpha,0)
            end
        end
    end
end

function Introduction()
    println("Wave Equation 3D: Multithreaded version") 
end

function Waveevo!(tau,nx,ny,nz,u,u₀,u₁)
    # 2nd order difference quotient (time) and discrete Laplace operator of 2th order (spatial) or, equivalently, an explicit Euler method of second order
    # From: https://beltoforion.de/en/recreational_mathematics/3d-wave-equation.php
    Threads.@threads for k=2:nz-1
        for j=2:ny-1
            for i=2:nx-1
                u[i,j,k] = tau*(u₁[i-1,j,k] + u₁[i+1,j,k] + u₁[i,j-1,k] + u₁[i,j+1,k]
                            + u₁[i,j,k-1] + u₁[i,j,k+1] - 6*u₁[i,j,k]) + 2*u₁[i,j,k] - u₀[i,j,k]
            end
        end
    end
end

function Boundcond!(kappa,nx,ny,nz,u,u₁)
    # Radiative/absorbing boundary conditions
    Threads.@threads for j=1:ny-1 for k=1:nz-1 u[1,j,k] = u₁[2,j,k] + ((kappa - 1)/(kappa + 1))*(u[2,j,k] - u₁[1,j,k]) end end
    Threads.@threads for j=1:ny-1 for k=1:nz-1 u[nx,j,k] = u₁[nx-1,j,k] + ((kappa - 1)/(kappa + 1))*(u[nx-1,j,k] - u₁[nx,j,k]) end end
    Threads.@threads for i=1:nx-1 for k=1:nz-1 u[i,1,k] = u₁[i,2,k] + ((kappa - 1)/(kappa + 1))*(u[i,2,k] - u₁[i,1,k]) end end
    Threads.@threads for i=1:nx-1 for k=1:nz-1 u[i,ny,k] = u₁[i,ny-1,k] + ((kappa - 1)/(kappa + 1))*(u[i,ny-1,k] - u₁[i,ny,k]) end end
    Threads.@threads for i=1:nx-1 for j=1:ny-1 u[i,j,1] = u₁[i,j,2] + ((kappa - 1)/(kappa + 1))*(u[i,j,2] - u₁[i,j,1]) end end
    Threads.@threads for i=1:nx-1 for j=1:ny-1 u[i,j,nz] = u₁[i,j,nz-1] + ((kappa - 1)/(kappa + 1))*(u[i,j,nz-1] - u₁[i,j,nz]) end end
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
    u = zeros(nx,ny,nz)
    u₀ = zeros(nx,ny,nz)
    u₁ = zeros(nx,ny,nz)
    Initialcond!(nx,ny,nz,u₀,u₁,x,y,z,alpha,dt)

    # viz (movie animation)
    domovie = nplot > 0
    if domovie
        gr(); ENV["GKSwstype"]="nul"        
        anim = Animation()
        println("Movie frame at iteration $iter"," umax=",maximum(u₁))
        surface(x,y,u₁[:,:,nz÷2],title="Time $t") # To plot the initial condition
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
            surface(x,y,u[:,:,nz÷2],title="Time "*string(round(t,digits=2)))
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