using
  PassiveTracerFlows,
  Test,
  Statistics,
  CUDA,
  Random

import # use 'import' rather than 'using' for submodules to keep namespace clean
  PassiveTracerFlows.TracerAdvectionDiffusion

# the devices on which tests will run
devices = CUDA.has_cuda() ? (CPU(), GPU()) : (CPU(),)

const rtol_traceradvectiondiffusion = 1e-12 # tolerance for rtol_traceradvdiff tests

# Run tests
testtime = @elapsed begin
  
for dev in devices
  
  println("testing on "*string(typeof(dev)))

  @testset "TracerAdvectionDiffusion" begin
    include("test_traceradvectiondiffusion.jl")

    stepper = "RK4"
    dt, nsteps  = 1e-2, 40
    @test test_constvel1D(stepper, dt, nsteps, dev)
    dt, tfinal  = 0.002, 0.1
    @test test_timedependentvel1D(stepper, dt, tfinal, dev)
    dt, nsteps  = 1e-2, 40
    @test test_constvel2D(stepper, dt, nsteps, dev)
    dt, tfinal  = 0.002, 0.1
    @test test_timedependentvel2D(stepper, dt, tfinal, dev)
    dt, nsteps  = 1e-2, 40
    @test test_constvel3D(stepper, dt, nsteps, dev)
    dt, tfinal  = 0.002, 0.1
    @test test_timedependentvel3D(stepper, dt, tfinal, dev)
    dt, tfinal  = 0.005, 0.1
    @test test_diffusion1D(stepper, dt, tfinal, dev; steadyflow=true)
    dt, tfinal  = 0.005, 0.1
    @test test_diffusion1D(stepper, dt, tfinal, dev; steadyflow=false)
    dt, tfinal  = 0.005, 0.1
    @test test_diffusion2D(stepper, dt, tfinal, dev; steadyflow=true)
    dt, tfinal  = 0.005, 0.1
    @test test_diffusion2D(stepper, dt, tfinal, dev; steadyflow=false)
    dt, tfinal  = 0.005, 0.1
    @test test_diffusion3D(stepper, dt, tfinal, dev; steadyflow=true)
    dt, tfinal  = 0.005, 0.1
    @test test_diffusion3D(stepper, dt, tfinal, dev; steadyflow=false)
    dt, tfinal  = 0.005, 0.1
    @test test_diffusion_multilayerqg(stepper, dt, tfinal, dev)
    dt, tfinal  = 0.005, 0.1
    @test test_hyperdiffusion(stepper, dt, tfinal, dev)
    
    @test TracerAdvectionDiffusion.noflow(π) == 0
  end
    
end

end #time

println("Total test time: ", testtime)
