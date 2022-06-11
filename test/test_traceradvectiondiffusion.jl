"""
    test_constvel1D(; kwargs...)

Advect a gaussian concentration `c0(x, t)` with a constant velocity flow
`u(x, y) = uvel` and compare the final state with
`cfinal = c0(x - uvel * tfinal)`.
"""
function test_constvel1D(stepper, dt, nsteps, dev::Device=CPU())

  nx, Lx = 128, 2π
  uvel = 0.05
  u(x) = uvel

  prob = TracerAdvectionDiffusion.Problem(dev, true; nx, Lx, κ=0.0, u, dt, stepper, steadyflow=true)
  sol, cl, vs, pr, gr = prob.sol, prob.clock, prob.vars, prob.params, prob.grid

  x = gr.x

  σ = 0.1
  c0ampl = 0.1
  c0func(x) = c0ampl * exp(-(x^2) / (2σ^2))

  c0 = c0func.(x)
  tfinal = nsteps * dt
  cfinal = @. c0func(x - uvel * tfinal)

  TracerAdvectionDiffusion.set_c!(prob, c0)

  stepforward!(prob, nsteps)
  TracerAdvectionDiffusion.updatevars!(prob)

  return isapprox(cfinal, vs.c, rtol = gr.nx*nsteps*1e-12)
end

"""
    test_timedependenttvel1D(; kwargs...)

Advect a gaussian concentration `c0(x, t)` with a time-varying velocity flow
`u(x, t) = uvel * sign(-t + tfinal/2)` and compares the final
state with `cfinal = c0(x - uvel * tfinal)`.
"""
function test_timedependentvel1D(stepper, dt, tfinal, dev::Device=CPU(); uvel = 0.05)
  
  nx, Lx = 128, 2π
  nsteps = round(Int, tfinal/dt)
  
  if !isapprox(tfinal, nsteps*dt, rtol=rtol_traceradvectiondiffusion)
    error("tfinal is not multiple of dt")
  end
  
  u(x, t) = uvel * t + uvel * dt/2

  prob = TracerAdvectionDiffusion.Problem(dev, true; nx, Lx, κ=0.0, u, dt, stepper)
  sol, cl, vs, pr, gr = prob.sol, prob.clock, prob.vars, prob.params, prob.grid
  x = gr.x

  σ = 0.2
  c0func(x) = 0.1 * exp(-(x^2) / (2σ^2))

  c0 = @. c0func(x)
  tfinal = nsteps * dt
  cfinal = @. c0func(x -  0.5uvel * tfinal^2)

  TracerAdvectionDiffusion.set_c!(prob, c0)

  stepforward!(prob, nsteps)
  TracerAdvectionDiffusion.updatevars!(prob)

  return isapprox(cfinal, vs.c, rtol = gr.nx*nsteps*1e-12)
end

"""
    test_constvel(; kwargs...)

Advect a gaussian concentration `c0(x, y, t)` with a constant velocity flow
`u(x, y) = uvel` and `v(x, y) = vvel` and compares the final state with
`cfinal = c0(x - uvel * tfinal, y - vvel * tfinal)`.
"""
function test_constvel(stepper, dt, nsteps, dev::Device=CPU())

  nx, Lx = 128, 2π
  uvel, vvel = 0.2, 0.1
  u(x, y) = uvel
  v(x, y) = vvel

  prob = TracerAdvectionDiffusion.Problem(dev; nx, Lx, κ=0.0, u, v, dt, stepper, steadyflow=true)
  sol, cl, vs, pr, gr = prob.sol, prob.clock, prob.vars, prob.params, prob.grid

  x, y = gridpoints(gr)

  σ = 0.1
  c0ampl = 0.1
  c0func(x, y) = c0ampl * exp(-(x^2 + y^2) / (2σ^2))

  c0 = c0func.(x, y)
  tfinal = nsteps * dt
  cfinal = @. c0func(x - uvel * tfinal, y - vvel * tfinal)

  TracerAdvectionDiffusion.set_c!(prob, c0)

  stepforward!(prob, nsteps)
  TracerAdvectionDiffusion.updatevars!(prob)

  return isapprox(cfinal, vs.c, rtol = gr.nx*gr.ny*nsteps*1e-12)
end


"""
    test_timedependenttvel(; kwargs...)

Advect a gaussian concentration `c0(x, y, t)` with a time-varying velocity flow
`u(x, y, t) = uvel` and `v(x, y, t) = vvel * sign(-t + tfinal/2)` and compares the final
state with `cfinal = c0(x - uvel * tfinal, y)`.
"""
function test_timedependentvel(stepper, dt, tfinal, dev::Device=CPU(); uvel = 0.5, αv = 0.5)
  
  nx, Lx = 128, 2π
  nsteps = round(Int, tfinal/dt)
  
  if !isapprox(tfinal, nsteps*dt, rtol=rtol_traceradvectiondiffusion)
    error("tfinal is not multiple of dt")
  end
  
  u(x, y, t) = uvel
  v(x, y, t) = αv * t + αv * dt/2

  prob = TracerAdvectionDiffusion.Problem(dev; nx, Lx, κ=0.0, u, v, dt, stepper)
  sol, cl, vs, pr, gr = prob.sol, prob.clock, prob.vars, prob.params, prob.grid
  x, y = gridpoints(gr)

  σ = 0.2
  c0func(x, y) = 0.1 * exp(-(x^2 + y^2) / (2σ^2))

  c0 = @. c0func(x, y)
  tfinal = nsteps * dt
  cfinal = @. c0func(x - uvel * tfinal, y - 0.5αv * tfinal^2)

  TracerAdvectionDiffusion.set_c!(prob, c0)

  stepforward!(prob, nsteps)
  TracerAdvectionDiffusion.updatevars!(prob)

  return isapprox(cfinal, vs.c, rtol = gr.nx*gr.ny*nsteps*1e-12)
end

"""
    test_diffusion1D(; kwargs...)

Diffuses a gaussian concentration c0(x, t) and compares the final state with
the analytic solution of the heat equation, cfinal
"""
function test_diffusion1D(stepper, dt, tfinal, dev::Device=CPU(); steadyflow = true)

  nx = 128
  Lx = 2π
   κ = 0.01
  nsteps = round(Int, tfinal/dt)

  if !isapprox(tfinal, nsteps*dt, rtol=rtol_traceradvectiondiffusion)
    error("tfinal is not multiple of dt")
  end

  prob = TracerAdvectionDiffusion.Problem(dev, true; steadyflow=steadyflow, nx=nx,
    Lx=Lx, κ=κ, dt=dt, stepper=stepper)
  sol, cl, vs, pr, gr = prob.sol, prob.clock, prob.vars, prob.params, prob.grid
  x = gr.x

  
  c0ampl, σ₀ = 0.1, 0.1
  σ(t) = sqrt(2κ * t + σ₀)
  c0func(x, t) = (c0ampl / σ(t)) * exp(-(x^2) / (2 * σ(t)^2))

  c0 = @. c0func(x, 0)
  tfinal = nsteps * dt
  cfinal = @. c0func(x, tfinal)

  TracerAdvectionDiffusion.set_c!(prob, ArrayType(dev)(c0))

  stepforward!(prob, nsteps)
  TracerAdvectionDiffusion.updatevars!(prob)

  return isapprox(cfinal, vs.c, rtol=gr.nx*nsteps*1e-12)
end

"""
    test_diffusion(; kwargs...)

Diffuses a gaussian concentration c0(x, y, t) and compares the final state with
the analytic solution of the heat equation, cfinal
"""
function test_diffusion(stepper, dt, tfinal, dev::Device=CPU(); steadyflow = true)

  nx = 128
  Lx = 2π
   κ = 0.01
  nsteps = round(Int, tfinal/dt)

  if !isapprox(tfinal, nsteps*dt, rtol=rtol_traceradvectiondiffusion)
    error("tfinal is not multiple of dt")
  end

  prob = TracerAdvectionDiffusion.Problem(dev; steadyflow=steadyflow, nx=nx,
    Lx=Lx, κ=κ, dt=dt, stepper=stepper)
  sol, cl, vs, pr, gr = prob.sol, prob.clock, prob.vars, prob.params, prob.grid
  x, y = gridpoints(gr)

  c0ampl, σ = 0.1, 0.1
  c0func(x, y) = c0ampl * exp(-(x^2 + y^2) / (2σ^2))

  c0 = @. c0func(x, y)
  tfinal = nsteps * dt
  σt = sqrt(2κ * tfinal + σ^2)
  cfinal = @. c0ampl * (σ^2 / σt^2) * exp(-(x^2 + y^2) / (2σt^2))

  TracerAdvectionDiffusion.set_c!(prob, c0)

  stepforward!(prob, nsteps)
  TracerAdvectionDiffusion.updatevars!(prob)

  return isapprox(cfinal, vs.c, rtol=gr.nx*gr.ny*nsteps*1e-12)
end

function test_diffusion_multilayerqg(stepper, dt, tfinal, dev::Device=CPU())
  # Set up MQGprob to generate zero flow and diffuse concentration field

  nx = 128
  Lx = 2π

  μ = 0                 
  β = 0                    

  nlayers = 2              
  f₀, g = 1, 1             
  H = [0.2, 0.8]          
  ρ = [4.0, 5.0]          

  U = zeros(nlayers) 
  U[1] = 0.0
  U[2] = 0.0

  MQGprob = MultiLayerQG.Problem(nlayers, dev;
                                 nx, Lx, f₀, g, H, ρ, U, μ, β, dt,
                                 stepper="FilteredRK4", aliased_fraction=0)
  grid = MQGprob.grid
  q₀ = zeros(dev, eltype(grid), (grid.nx, grid.ny, nlayers))
  
  MultiLayerQG.set_q!(MQGprob, q₀)
  
  κ = 0.01
  tracer_release_time = dt * 50

  nsteps = round(Int, tfinal/dt)
  if !isapprox(tfinal, nsteps*dt, rtol=rtol_traceradvectiondiffusion)
    error("tfinal is not multiple of dt")
  end

  ADprob = TracerAdvectionDiffusion.Problem(dev, MQGprob; κ, stepper, tracer_release_time)
  sol, cl, vs, pr, gr = ADprob.sol, ADprob.clock, ADprob.vars, ADprob.params, ADprob.grid
  x, y = gridpoints(gr)

  c0ampl, σ = 0.1, 0.1
  c0func(x, y) = c0ampl * exp(-(x^2 + y^2) / (2σ^2))
  c0 = @. c0func(x, y)
  tfinal = nsteps * dt
  σt = sqrt(2κ * tfinal + σ^2)
  cfinal = @. c0ampl * (σ^2 / σt^2) * exp(-(x^2 + y^2) / (2σt^2))

  TracerAdvectionDiffusion.set_c!(ADprob, c0)

  stepforward!(ADprob, nsteps)
  TracerAdvectionDiffusion.updatevars!(ADprob)

  # Compare to analytic solution
  return isapprox(cfinal, vs.c[:, :, 1], rtol = gr.nx*gr.ny*nsteps*1e-12)  &&
         isapprox(cfinal, vs.c[:, :, 2], rtol = gr.nx*gr.ny*nsteps*1e-12)
end

"""
    test_hyperdiffusion(; kwargs...)

Diffuses a gaussian concentration c0(x, y, t) using hyperdiffusivity and
compares the final state with the analytic solution of the heat equation, cfinal
"""
function test_hyperdiffusion(stepper, dt, tfinal, dev::Device=CPU(); steadyflow = true)

   nx = 128
   Lx = 2π
    κ = 0.0   # no diffusivity
    η = κ     # no diffusivity
   κh = 0.01  # hyperdiffusivity coeff
  nκh = 1     # nκh=1 converts hyperdiffusivity to plain diffusivity
              # so we can compare with the analytic solution of heat equation

  nsteps = round(Int, tfinal/dt)

  if !isapprox(tfinal, nsteps*dt, rtol=rtol_traceradvectiondiffusion)
    error("tfinal is not multiple of dt")
  end

  gr = TwoDGrid(dev, nx, Lx)
  x, y = gridpoints(gr)

  #u, v = zero(x), zero(x) #0*x, 0*x
  u(x, y) = 0.0
  v(x, y) = 0.0

  vs = TracerAdvectionDiffusion.Vars(dev, gr)
  pr = TracerAdvectionDiffusion.ConstDiffSteadyFlowParams(η, κ, κh, nκh, u, v, gr)
  eq = TracerAdvectionDiffusion.Equation(dev, pr, gr)
  prob = FourierFlows.Problem(eq, stepper, dt, gr, vs, pr, dev)

  c0ampl, σ = 0.1, 0.1
  c0func(x, y) = c0ampl * exp(-(x^2 + y^2) / (2σ^2))

  c0 = @. c0func(x, y)
  tfinal = nsteps * dt
  σt = sqrt(2κh * tfinal + σ^2)
  cfinal = @. c0ampl * σ^2 / σt^2 * exp(-(x^2 + y^2) / (2σt^2))

  TracerAdvectionDiffusion.set_c!(prob, c0)

  stepforward!(prob, nsteps)
  TracerAdvectionDiffusion.updatevars!(prob)

  return isapprox(cfinal, vs.c, rtol = gr.nx*gr.ny*nsteps*1e-12)
end
