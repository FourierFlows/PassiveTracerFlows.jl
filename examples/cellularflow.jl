# # Advection-diffusion of tracer by cellular flow
#
#md # This example can be viewed as a Jupyter notebook via [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/cellularflow.ipynb).
# 
# An example demonstrating the advection-diffusion of a tracer by a cellular flow.
#
# ## Install dependencies
#
# First let's make sure we have all required packages installed.

# ```julia
# using Pkg
# pkg"add PassiveTracerFlows, Plots, Printf"
# ```

# ## Let's begin
# Let's load `PassiveTracerFlows.jl` and some other needed packages.
#
using PassiveTracerFlows, Plots, Printf


# ## Choosing a device: CPU or GPU

dev = CPU()     # Device (CPU/GPU)
nothing # hide


# ## Numerical parameters and time-stepping parameters

      n = 128            # 2D resolution = n²
stepper = "RK4"          # timestepper
     dt = 0.02           # timestep
 nsteps = 800            # total number of time-steps
 nsubs  = 25             # number of time-steps for intermediate logging/plotting (nsteps must be multiple of nsubs)
nothing # hide


# ## Numerical parameters and time-stepping parameters

L = 2π        # domain size
κ = 0.002     # diffusivity
nothing # hide


# ## Set up cellular flow
# We create a two-dimensional grid to construct the cellular flow. Our cellular flow is derived
# from a streamfunction ``ψ(x, y) = ψ₀ \cos(x) \cos(y)`` as ``(u, v) = (-∂_y ψ, ∂_x ψ)``.

grid = TwoDGrid(n, L)
x, y = gridpoints(grid)

ψ₀ = 0.2
mx, my = 1, 1

ψ = @. ψ₀ * cos(mx * x) * cos(my * y)

uvel(x, y) =  ψ₀ * mx * cos(mx * x) * sin(my * y)
vvel(x, y) = -ψ₀ * my * sin(mx * x) * cos(my * y)
nothing # hide


# ## Problem setup
# We initialize a `Problem` by providing a set of keyword arguments.
prob = TracerAdvectionDiffusion.Problem(dev; nx=n, Lx=L, κ=κ, steadyflow=true, u=uvel, v=vvel,
                                          dt=dt, stepper=stepper)
nothing # hide

# and define some shortcuts
sol, clock, vars, params, grid = prob.sol, prob.clock, prob.vars, prob.params, prob.grid
x, y = grid.x, grid.y
nothing # hide


# ## Setting initial conditions

# Our initial condition for the tracer ``c`` is a gaussian centered at ``(x, y) = (L_x/5, 0)``.

gaussian(x, y, σ) = exp(-(x^2 + y^2) / (2σ^2))

amplitude, spread = 0.5, 0.15
c₀ = [amplitude * gaussian(x[i] - 0.2 * grid.Lx, y[j], spread) for i=1:grid.nx, j=1:grid.ny]

TracerAdvectionDiffusion.set_c!(prob, c₀)
nothing # hide

# Let's plot the initial tracer concentration and streamlines. Note that when plotting, we decorate 
# the variable to be plotted with `Array()` to make sure it is brought back on the CPU when 
# `vars` live on the GPU.

function plot_output(prob)
  c = prob.vars.c
  
  p = heatmap(x, y, Array(vars.c'),
         aspectratio = 1,
              c = :balance,
         legend = :false,
           clim = (-0.2, 0.2),
          xlims = (-grid.Lx/2, grid.Lx/2),
          ylims = (-grid.Ly/2, grid.Ly/2),
         xticks = -3:3,
         yticks = -3:3,
         xlabel = "x",
         ylabel = "y",
          title = "initial tracer concentration (shading) + streamlines",
     framestyle = :box)

  contour!(p, x, y, Array(ψ'),
     levels=0.0125:0.025:0.2,
     lw=2, c=:black, ls=:solid, alpha=0.7)

  contour!(p, x, y, Array(ψ'),
     levels=-0.1875:0.025:-0.0125,
     lw=2, c=:black, ls=:dash, alpha=0.7)

  return p
end
nothing # hide


# ## Time-stepping the `Problem` forward

# We time-step the `Problem` forward in time.

startwalltime = time()

p = plot_output(prob)

anim = @animate for j = 0:round(Int, nsteps/nsubs)

 if j % (200 / nsubs) == 0
    log = @sprintf("step: %04d, t: %d, walltime: %.2f min", 
                   clock.step, clock.t, (time()-startwalltime)/60)
    
    println(log)
  end  

  p[1][1][:z] = Array(vars.c)
  p[1][:title] = "concentration, t=" * @sprintf("%.2f", clock.t)

  stepforward!(prob, nsubs)
  TracerAdvectionDiffusion.updatevars!(prob)
end

mp4(anim, "cellularflow.mp4", fps=12)
