"""
    test_constvel(; kwargs...)

Advects a gaussian concentration c0(x, y, t) with a constant velocity flow
u(x, y) = uvel and v(x, y) = vvel and compares the final state with
cfinal = c0(x-uvel*tfinal, y-vvel*tfinal)
"""
function test_constvel(stepper, dt, nsteps, dev::Device=CPU())

  nx, Lx = 128, 2π
  uvel, vvel = 0.2, 0.1
  u(x, y) = uvel
  v(x, y) = vvel

  prob = TracerAdvectionDiffusion.Problem(dev; nx=nx, Lx=Lx, κ=0.0, u=u, v=v, dt=dt, stepper=stepper, steadyflow=true)
  sol, cl, vs, pr, gr = prob.sol, prob.clock, prob.vars, prob.params, prob.grid

  x, y = gridpoints(gr)

  σ = 0.1
  c0func(x, y) = 0.1*exp(-(x^2+y^2)/(2σ^2))

  c0 = c0func.(x, y)
  tfinal = nsteps*dt
  cfinal = @. c0func(x - uvel*tfinal, y - vvel*tfinal)

  TracerAdvectionDiffusion.set_c!(prob, c0)

  stepforward!(prob, nsteps)
  TracerAdvectionDiffusion.updatevars!(prob)

  isapprox(cfinal, vs.c, rtol=gr.nx*gr.ny*nsteps*1e-12)
end


"""
    test_timedependenttvel(; kwargs...)

Advects a gaussian concentration c0(x, y, t) with a time-varying velocity flow
u(x, y, t) = uvel and v(x, y, t) = vvel*sign(-t+tfinal/2) and compares the final
state with cfinal = c0(x-uvel*tfinal, y)
"""
function test_timedependentvel(stepper, dt, tfinal, dev::Device=CPU(); uvel=0.5, αv=0.5)
  
  nx, Lx = 128, 2π
  nsteps = round(Int, tfinal/dt)
  
  if !isapprox(tfinal, nsteps*dt, rtol=rtol_traceradvectiondiffusion)
    error("tfinal is not multiple of dt")
  end
  
  u(x, y, t) = uvel
  v(x, y, t) = αv*t + αv*dt/2

  prob = TracerAdvectionDiffusion.Problem(dev; nx=nx, Lx=Lx, κ=0.0, u=u, v=v, dt=dt, stepper=stepper)
  sol, cl, vs, pr, gr = prob.sol, prob.clock, prob.vars, prob.params, prob.grid
  x, y = gridpoints(gr)

  σ = 0.2
  c0func(x, y) = 0.1*exp(-(x^2+y^2)/(2σ^2))

  c0 = @. c0func(x, y)
  tfinal = nsteps*dt
  cfinal = @. c0func(x - uvel*tfinal, y - 0.5αv*tfinal^2)

  TracerAdvectionDiffusion.set_c!(prob, c0)

  stepforward!(prob, nsteps)
  TracerAdvectionDiffusion.updatevars!(prob)
  isapprox(cfinal, vs.c, rtol=gr.nx*gr.ny*nsteps*1e-12)
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
  c0func(x, y) = c0ampl*exp(-(x^2+y^2)/(2σ^2))

  c0 = @. c0func(x, y)
  tfinal = nsteps*dt
  σt = sqrt(2*κ*tfinal + σ^2)
  cfinal = @. c0ampl*(σ^2/σt^2)*exp(-(x^2+y^2)/(2*σt^2))

  TracerAdvectionDiffusion.set_c!(prob, c0)

  stepforward!(prob, nsteps)
  TracerAdvectionDiffusion.updatevars!(prob)

  isapprox(cfinal, vs.c, rtol=gr.nx*gr.ny*nsteps*1e-12)
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

  u, v = zero(x), zero(x) #0*x, 0*x

  vs = TracerAdvectionDiffusion.Vars(dev, gr)
  pr = TracerAdvectionDiffusion.ConstDiffSteadyFlowParams(η, κ, κh, nκh, u, v)
  eq = TracerAdvectionDiffusion.Equation(pr, gr)
  prob = FourierFlows.Problem(eq, stepper, dt, gr, vs, pr, dev)

  c0ampl, σ = 0.1, 0.1
  c0func(x, y) = c0ampl*exp(-(x^2+y^2)/(2σ^2))

  c0 = @. c0func(x, y)
  tfinal = nsteps*dt
  σt = sqrt(2*κh*tfinal + σ^2)
  cfinal = @. c0ampl*σ^2/σt^2 * exp(-(x^2+y^2)/(2*σt^2))

  TracerAdvectionDiffusion.set_c!(prob, c0)

  stepforward!(prob, nsteps)
  TracerAdvectionDiffusion.updatevars!(prob)

  isapprox(cfinal, vs.c, rtol=gr.nx*gr.ny*nsteps*1e-12)
end
