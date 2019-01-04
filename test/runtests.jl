#!/usr/bin/env julia

using
  FourierFlows,
  Test,
  Statistics,
  Random,
  FFTW

import # use 'import' rather than 'using' for submodules to keep namespace clean
  PassiveTracerFlows.TracerAdvDiff

const rtol_traceradvdiff = 1e-12 # tolerance for rtol_traceradvdiff tests

# Run tests
testtime = @elapsed begin

@testset "TracerAdvDiff" begin
  include("test_traceradvdiff.jl")

  stepper = "RK4"
  dt, nsteps  = 1e-2, 40
  @test test_constvel(stepper, dt, nsteps)
  dt, tfinal  = 0.002, 0.1
  @test test_timedependentvel(stepper, dt, tfinal)
  dt, tfinal  = 0.005, 0.1
  @test test_diffusion(stepper, dt, tfinal; steadyflow=true)
  dt, tfinal  = 0.005, 0.1
  @test test_diffusion(stepper, dt, tfinal; steadyflow=false)
  dt, tfinal  = 0.005, 0.1
  @test test_hyperdiffusion(stepper, dt, tfinal)
end

end
println("Total test time: ", testtime)
