# # Advection-diffusion of tracer in one dimension
#
# This is an example demonstrating the advection-diffusion of a passive tracer in 
# one dimension. 
#
# ## Install dependencies
#
# First let's make sure we have all the required packages installed

# ```julia
# using Pkg
# pkg.add(["PassiveTracerFlows", "Printf", "Plots", "JLD2"])
# ```
#
# ## Let's begin
# First load packages needed to run this example.
using PassiveTracerFlows, Plots, Printf, JLD2, LinearAlgebra

# ## Choosing a device: CPU or GPU

dev = CPU()     # Device (CPU/GPU)
nothing # hide


# ## Numerical parameters and time-stepping parameters

      n = 128            # 2D resolution = n²
stepper = "RK4"          # timestepper
     dt = 0.02           # timestep
 nsteps = 5000           # total number of time-steps
nothing # hide


# ## Physical parameters

L = 2π       # domain size
κ = 0.01     # diffusivity
nothing # hide

# ## Flow
# We set a constant background flow and pass this to `OneDAdvectingFlow` with `steadyflow = true` to indicate the flow is not time dependent.
u(x) = 0.05
advecting_flow = OneDAdvectingFlow(; u, steadyflow = true)

# ## Problem setup
# We initialize a `Problem` by providing a set of keyword arguments.
prob = TracerAdvectionDiffusion.Problem(dev, advecting_flow; nx=n, Lx=L, κ, dt, stepper)
nothing # hide

# and define some shortcuts
sol, clock, vars, params, grid = prob.sol, prob.clock, prob.vars, prob.params, prob.grid
x = grid.x

# ## Initial condition
# We advect-diffuse a concentration field that has an initial concentration set to Gaussian.
gaussian(x, σ) = exp(-x^2 / 2σ^2)

amplitude, spread = 1, 0.15
c₀ = [amplitude * gaussian(x[i], spread) for i in 1:grid.nx]

TracerAdvectionDiffusion.set_c!(prob, c₀)
nothing #hide

# ## Saving output
# We create the saved output using the `Output` function from `FourierFlows.jl` then
# save the concentration field using the `get_concentration` function every 50 timesteps.
function get_concentration(prob)
  ldiv!(prob.vars.c, prob.grid.rfftplan, deepcopy(prob.sol))

  return prob.vars.c
end
  
output = Output(prob, "advection-diffusion1D.jld2",
                (:concentration, get_concentration))

# By calling `saveproblem(output)` we save information that we will use for plotting later on.
saveproblem(output)

# ## Stepping the problem forward
# Now we step the problem forward and save output every 50 timesteps.

save_frequency = 50 # frequency at which output is saved

startwalltime = time()
while clock.step <= nsteps
  if clock.step % save_frequency == 0
    saveoutput(output)
    log = @sprintf("Output saved, step: %04d, t: %.2f, walltime: %.2f min",
                   clock.step, clock.t, (time()-startwalltime) / 60)

    println(log)
  end

  stepforward!(prob)
end

# ## Visualising the output
# We load the `.jld2` file and create a timeseries of the concentration field
file = jldopen(output.path)

iterations = parse.(Int, keys(file["snapshots/t"]))

t = [file["snapshots/t/$i"] for i ∈ iterations]
c = [file["snapshots/concentration/$i"] for i ∈ iterations]
nothing # hide

# Set up the plotting arguments and look at the initial concentration.
x, Lx  = file["grid/x"], file["grid/Lx"]

plot_args = (xlabel = "x",
             ylabel = "c",
             framestyle = :box,
             xlims = (-Lx/2, Lx/2),
             ylims = (0, maximum(c[1])),
             legend = :false,
             color = :balance)

p = plot(x, Array(c[1]);
         title = "concentration, t = " * @sprintf("%.2f", t[1]),
         plot_args...)
nothing # hide

# Now, we create a movie of the tracer concentration being advected and  diffused.

anim = @animate for i ∈ 1:length(t)
  p[1][:title] = "concentration, t = " * @sprintf("%.2f", t[i])
  p[1][1][:y] = Array(c[i])
end

mp4(anim, "1D_advection-diffusion.mp4", fps = 12)
