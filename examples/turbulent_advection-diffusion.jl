# # Advection-diffusion of tracer by a turbulent flow
#
# This is an example demonstrating the advection-diffusion of a tracer using a 
# turbulent flow generated by the `GeophysicalFlows.jl` package.
#
# ## Install dependencies
#
# First let's make sure we have all the required packages installed

# ```julia
# using Pkg
# pkg.add(["PassiveTracerFlows", "Printf", "Plots", "JLD2"])
#
# ## Let's begin
# First load `PassiveTracerFlows.jl` and the other packages needed to run this example.
using PassiveTracerFlows, Printf, Plots, JLD2, LinearAlgebra
using Random: seed!

# ## Choosing a device: CPU or GPU
dev = CPU()
nothing # hide

# ## Setting up a `MultiLayerQG.Problem` to generate a turbulent flow
#
# The tubulent flow we will use to advect the passive tracer is generated using the 
# [`MultiLayerQG`](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/modules/multilayerqg/) module from the [`GeophysicalFlows.jl`](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/) package. A more detailed setup  
# of this two layer system can be found [here](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/literated/multilayerqg_2layer/).
#
# ### Numerical and time stepping parameters for the flow

      n = 128            # 2D resolution = n²
stepper = "FilteredRK4"  # timestepper
     dt = 2.5e-3         # timestep

# Physical parameters 
L = 2π                   # domain size
μ = 5e-2                 # bottom drag
β = 5                    # the y-gradient of planetary PV

nlayers = 2              # number of layers
f₀, g = 1, 1             # Coriolis parameter and gravitational constant
 H = [0.2, 0.8]          # the rest depths of each layer
 ρ = [4.0, 5.0]          # the density of each layer

 U = zeros(nlayers) # the imposed mean zonal flow in each layer
 U[1] = 1.0
 U[2] = 0.0

# ### `MultiLayerQG.Problem` setup, shortcuts and initial conditions
MQGprob = MultiLayerQG.Problem(nlayers, dev;
                        nx=n, Lx=L, f₀=f₀, g=g, H=H, ρ=ρ, U=U, μ=μ, β=β,
                        dt=dt, stepper=stepper, aliased_fraction=0)
grid = MQGprob.grid
x, y = grid.x, grid.y

# Initial conditions                        
seed!(1234) # reset of the random number generator for reproducibility
q₀  = 1e-2 * ArrayType(dev)(randn((grid.nx, grid.ny, nlayers)))
q₀h = MQGprob.timestepper.filter .* rfft(q₀, (1, 2)) # apply rfft  only in dims=1, 2
q₀  = irfft(q₀h, grid.nx, (1, 2))                 # apply irfft only in dims=1, 2

MultiLayerQG.set_q!(MQGprob, q₀)
nothing 

# ## Tracer advection-diffusion setup
# Now that we have a `MultiLayerQG.Problem` setup to generate our turbulent flow, we
# setup an advection-diffusion simulation. This is done by passing the `MultiLayerQG.Problem`
# as an argument to `TracerAdvectionDiffusion.Problem` which sets up an advection-diffusion problem
# with same parameters where applicable. We also need to pass a value for the constant diffusivity `κ`,
# the `stepper` used to step the problem forward and when we want the tracer released into the flow.
# We will let the flow run until it reaches a statistical equilibrium and then advect-diffuse the tracer.

κ = 0.002
nsteps = 4000               # total number of time-steps
nsubs = 1                   # number of steps the simulation takes at each iteration 
tracer_release = dt * 8000  # run flow for some time before releasing tracer

ADprob = TracerAdvectionDiffusion.Problem(dev, MQGprob; κ, stepper, tracer_release)
nothing

# ## Initial condition for concentration in both layers
# We have a two layer system so we will advect-diffuse the tracer in both layers.
# To do this we set the initial condition for tracer concetration as a Gaussian centred at the origin.
# Then we create some shortcuts for the `TracerAdvectionDiffusion.Problem`.
gaussian(x, y, σ) = exp(-(x^2 + y^2) / (2σ^2))

amplitude, spread = 10, 0.15
c₀ = [amplitude * gaussian(x[i], y[j], spread) for i=1:grid.nx, j=1:grid.ny]

TracerAdvectionDiffusion.set_c!(ADprob, c₀, nlayers)

# Shortcuts for advection-diffusion problem
sol, clock, vars, params, grid = ADprob.sol, ADprob.clock, ADprob.vars, ADprob.params, ADprob.grid
x, y = grid.x, grid.y

# ## Saving output
# The parent package `FourierFlows.jl` provides the functionality to save the output from our simulation.
# To do this we write a function `get_concentration` and pass this to the `Output` function along 
# with the `TracerAdvectionDiffusion.Problem` and the name of the output file.

function get_concentration(prob)
    ldiv!(prob.vars.c, prob.grid.rfftplan, deepcopy(prob.sol))
    return prob.vars.c
end
output = Output(ADprob, "advection-diffusion.jld2", (:concentration, get_concentration))
# This saves information that we will use for plotting later on
saveproblem(output)

# ## Step the problem forward and save the output
# We specify that we would like to save the concentration every 50 timesteps using `save_frequency`
# then step the problem forward.
save_frequency = 50 # Freqeuncy at which output is saved

startwalltime = time()
while clock.step <= nsteps

    if clock.step % save_frequency == 0

       saveoutput(output)
       log = @sprintf("Output saved, step: %04d, t: %d, walltime: %.2f min",
                      clock.step, clock.t, (time()-startwalltime)/60)
   
       println(log)
     end
   
     stepforward!(ADprob)
     TracerAdvectionDiffusion.MQGupdatevars!(ADprob)
end

# Append this information to our saved data for plotting later on
jldopen(output.path, "a+") do path
    path["save_frequency"] = save_frequency
    path["final_step"] = ADprob.clock.step - 1
end

# ## Visualising the output
# We now have output from our simulation saved in `advection-diffusion.jld2`.
# From this we can create a time series for the tracer that has been advected-diffused
# in the lower layer of our turbulent flow (the simulation from the upper layer can be obtained in a similar manner).

# Create time series for the concentration in the upper layer
conc_data = load("advection-diffusion.jld2")

saved_data = 0:conc_data["save_frequency"]:conc_data["final_step"]
t = [conc_data["snapshots/t/"*string(i)] for i ∈ saved_data]
# Concentration time series in the lower layer
cₗ = [abs.(conc_data["snapshots/concentration/"*string(i)][:, :, 2]) for i ∈ saved_data]

x, y,  = conc_data["grid/x"], conc_data["grid/y"]
Lx, Ly = conc_data["grid/Lx"], conc_data["grid/Ly"]
plot_args = (xlabel = "x",
             ylabel = "y",
             aspectratio = 1,
             framestyle = :box,
             xlims = (-Lx/2, Lx/2),
             ylims = (-Ly/2, Ly/2),
             colorbar = true,
             colorbar_title = " \nConcentration",
             color = :deep)

p = heatmap(x, y, cₗ[1]', title = "Concentration, t = " * @sprintf("%.2f", t[1]); plot_args...)
conc_anim = @animate for i ∈ 2:length(t)

    heatmap!(p, x, y, cₗ[i]', title = "Concentration, t = $(t[i])"; plot_args...)

end

# Create a movie of the tracer
mp4(conc_anim, "conc_adv-diff.mp4", fps = 12)
