module PassiveTracerFlows

using
  CUDA,
  Reexport
  
@reexport using FourierFlows

include("traceradvectiondiffusion.jl")

@reexport using PassiveTracerFlows.TracerAdvectionDiffusion

end # module
