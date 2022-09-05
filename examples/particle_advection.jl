# # Particle advection in two dimensions
#
# This is an example demonstrating the advection of particles in
# two dimensions. The flow field used is a *cat's eye flow*[^1].
#
# ## Install dependencies
#
# First let's make sure we have all the required packages installed

# ```julia
# using Pkg
# pkg.add(["PassiveTracerFlows", "CairoMakie", "JLD2"])
# ```
#
# ## Let's begin
# First load packages needed to run this example.
using PassiveTracerFlows, CairoMakie, Printf, JLD2, LinearAlgebra, Random

# ## Choosing a device: CPU or GPU

dev = CPU()     # Device (CPU/GPU)
nothing # hide

# ## Numerical parameters and time-stepping parameters

     nx = 128            # 2D resolution = nx²
stepper = "RK4"          # timestepper
     dt = 0.02           # timestep
 nsteps = 800            # total number of time-steps
 nsubs  = 25             # number of time-steps for intermediate logging/plotting (nsteps must be multiple of nsubs)
nothing # hide


# ## Numerical parameters and time-stepping parameters

Lx = 2π       # domain size
 κ = 0        # diffusivity, makes this an advection only problem
nothing # hide


# ## Set up a cats eye flow
# We create a two-dimensional grid to construct the cats eye flow. The cats eye flow is derived
# from a streamfunction ``ψ(x, y) = \sin(x)\sin(y) + \eps \cos(x) \cos(y)`` as ``(u, v) = (-∂_y ψ, ∂_x ψ)``.
# The cats eye flow is then passed into the `TwoDAdvectingFlow` constructor with `steadyflow = true`
# to indicate that the flow is not time dependent.
grid = TwoDGrid(dev; nx, Lx)

ϵ = 0.1

ψ(x, y) = sin(x) * sin(y) + ϵ * cos(x) * cos(y)

u(x, y) = -sin(x) * cos(y) + ϵ * cos(x) * sin(y)
v(x, y) =  cos(x) * sin(y) - ϵ * sin(y) * cos(x)

advecting_flow = TwoDAdvectingFlow(; u, v, steadyflow = true)
nothing # hide

# ## Problem setup
# We initialize a `Problem` by providing a set of keyword arguments. Again note that this will be an 
# advection only problem as we have set the diffusivity, ``\kappa``, to zero
prob = TracerAdvectionDiffusion.Problem(dev, advecting_flow; nx, Lx, κ, dt, stepper)

# and define some shortcuts
sol, clock, vars, params, grid = prob.sol, prob.clock, prob.vars, prob.params, prob.grid
x, y = grid.x, grid.y
xgrid, ygrid = gridpoints(grid)
nothing # hide

# ## Setting initial conditions
# We choose some random locations on the grid to place particles that will be advected by the flow.

pᵢ = rand(1:grid.nx, (5, 2))
p₀ = zeros(grid.nx, grid.ny)
for j ∈ 1:lastindex(pᵢ[:, 1])
  p₀[pᵢ[j, 1], pᵢ[j, 2]] = 1.0
end

TracerAdvectionDiffusion.set_c!(prob, p₀)
nothing # hide

# ## Time-stepping the `Problem` forward

# We want to step the `Problem` forward in time and, whilst doing so, we'd like
# to produce an animation of the particles being advected.
#
# First we create a figure using [`Observable`](https://makie.juliaplots.org/stable/documentation/nodes/)s.

c_anim = Observable(Array(vars.c))
title = Observable(@sprintf("particle position, t = %.2f", clock.t))

Lx, Ly = grid.Lx, grid.Ly

fig = Figure(resolution = (600, 600))

ax = Axis(fig[1, 1], 
          xlabel = "x",
          ylabel = "y",
          aspect = 1,
          title = title,
          limits = ((-Lx/2, Lx/2), (-Ly/2, Ly/2)))

hm = heatmap!(ax, x, y, c_anim;
              colormap = :balance)

contour!(ax, x, y, ψ.(xgrid, ygrid);
         color = :grey, linestyle = :solid)

fig

# Now we time-step `Problem` and update the `c_anim` and `title` observables as we go
# to create an animation.

startwalltime = time()

frames = 0:round(Int, nsteps/nsubs)
record(fig, "particle-advection.mp4", frames, framerate = 12) do j
   if j % (200 / nsubs) == 0
      log = @sprintf("step: %04d, t: %d, walltime: %.2f min", 
                     clock.step, clock.t, (time()-startwalltime)/60)
      
      println(log)
    end

  c_anim[] = vars.c
  title[] = @sprintf("particle position, t = %.2f", clock.t)

  stepforward!(prob, nsubs)
  TracerAdvectionDiffusion.updatevars!(prob)
end

# ![](particle-advection.mp4)

## References

# [^1] S. Childress and A. M. Soward, [Scalar transport and alpha-effect for a family of cat’s-eye flows](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/scalar-transport-and-alphaeffect-for-a-family-of-catseye-flows/A649BCF133BB25DA3B80239024945C77), J. Fluid Mech. 205, 99 (1989).