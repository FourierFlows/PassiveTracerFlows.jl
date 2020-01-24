#!/usr/bin/env julia

using
  FourierFlows,
  Test,
  Statistics,
  Random,
  FFTW

import # use 'import' rather than 'using' for submodules to keep namespace clean
  PassiveTracerFlows.TracerAdvDiff

# the devices on which tests will run
devices = (CPU(),)
@has_cuda devices = (CPU(), GPU())
@has_cuda using CuArrays

const rtol_traceradvdiff = 1e-12 # tolerance for rtol_traceradvdiff tests

# Run tests
testtime = @elapsed begin
  
for dev in devices
  
  println("testing on "*string(typeof(dev)))

  @testset "TracerAdvDiff" begin
    include("test_traceradvdiff.jl")

    @test test_noflow()
    
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

end #time

println("Total test time: ", testtime)
