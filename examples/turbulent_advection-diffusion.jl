using PassiveTracerFlows, Printf, JLD2, Plots

using Random: seed!

# Choosing a device: CPU or GPU
dev = CPU()

# Numerical and time stepping parameters
n = 128                  # 2D resolution = n²
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

# Problem setup and shortcuts
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

# Tracer advection-diffusion setup
κ = 0.002
nsteps = 4000               # total number of time-steps
nsubs = 1                   # number of steps the simulation takes at each iteration 
tracer_release = dt * 8000  # run flow for some time before releasing tracer

ADprob = TracerAdvectionDiffusion.Problem(dev, MQGprob; κ = κ, stepper = stepper, tracer_release = tracer_release)

# Initial condition for concentration in both layers
gaussian(x, y, σ) = exp(-(x^2 + y^2) / (2σ^2))

amplitude, spread = 1, 0.15
c₀ = [gaussian(x[i], y[j], spread) for i=1:grid.nx, j=1:grid.ny]

TracerAdvectionDiffusion.set_c!(ADprob, c₀, nlayers)

# Shortcuts for advection-diffusion problem
sol, clock, vars, params, grid = ADprob.sol, ADprob.clock, ADprob.vars, ADprob.params, ADprob.grid
x, y = grid.x, grid.y

# Function to extract concentration field and create file to save output
function GetConcentration(prob)
    Concentration = @. prob.vars.c
    return Concentration
end
output = Output(ADprob, "advection-diffusion.jld2", (:Concentration, GetConcentration))
saveproblem(output)

# Step the problem forward and save the output

save_frequency = 100 # Freqeuncy at which output is saved

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

jldopen(output.path, "a+") do path
    path["save_frequency"] = save_frequency
    path["final_step"] = ADprob.clock.step - 1
end

# Create time series for the concentration in the upper layer
conc_data = load("advection-diffusion.jld2")

saved_data = 0:conc_data["save_frequency"]:conc_data["final_step"]
t = [conc_data["snapshots/t/"*string(i)] for i ∈ saved_data]
Cᵤ = [abs.(conc_data["snapshots/Concentration/"*string(i)][:, :, 2]) for i ∈ saved_data]

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
climits = (minimum(Cᵤ[1]), maximum(Cᵤ[1]))
p = heatmap(x, y, Cᵤ[1]', title = "Concentration, t = $(t[1])", clims = climits; plot_args...)
conc_anim = @animate for i ∈ 2:length(t)

    climits = (minimum(Cᵤ[i]), maximum(Cᵤ[i]))
    heatmap!(p, x, y, Cᵤ[i]', title = "Concentration, t = $(t[i])", clims = climits; plot_args...)

end

mp4(conc_anim, "conc_adv-diff.mp4", fps = 12)