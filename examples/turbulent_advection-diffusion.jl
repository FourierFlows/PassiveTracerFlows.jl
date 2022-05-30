using PassiveTracerFlows

using Random: seed!

# Choosing a device: CPU or GPU
dev = CPU()

# Numerical and time stepping parameters
n = 128                  # 2D resolution = n²
stepper = "FilteredRK4"  # timestepper
     dt = 2.5e-3         # timestep
 nsteps = 20000          # total number of time-steps
 nsubs  = 50             # number of time-steps for plotting (nsteps must be multiple of nsubs)

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
sol, clock, params, vars, grid = MQGprob.sol, MQGprob.clock, MQGprob.params, MQGprob.vars, MQGprob.grid
x, y = grid.x, grid.y

# Initial conditions                        
seed!(1234) # reset of the random number generator for reproducibility
q₀  = 1e-2 * ArrayType(dev)(randn((grid.nx, grid.ny, nlayers)))
q₀h = MQGprob.timestepper.filter .* rfft(q₀, (1, 2)) # apply rfft  only in dims=1, 2
q₀  = irfft(q₀h, grid.nx, (1, 2))                 # apply irfft only in dims=1, 2

MultiLayerQG.set_q!(MQGprob, q₀)

# Tracer advection-diffusion setup
κ = 0.002

ADprob = TracerAdvectionDiffusion.Problem(dev, MQGprob; κ = κ, stepper = stepper)

gaussian(x, y, σ) = exp(-(x^2 + y^2) / (2σ^2))

amplitude, spread = 0.5, 0.15
c₀ = [amplitude * gaussian(x[i], y[j], spread) for i=1:grid.nx, j=1:grid.ny]

TracerAdvectionDiffusion.set_c!(ADprob, nlayers, c₀)