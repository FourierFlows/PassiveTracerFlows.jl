module TracerAdvectionDiffusion

export
   Problem,
   set_c!,
   updatevars!

using
   DocStringExtensions,
   Reexport

@reexport using FourierFlows, GeophysicalFlows.MultiLayerQG

using GeophysicalFlows.MultiLayerQG: SingleLayerParams, TwoLayerParams, numberoflayers

import LinearAlgebra: mul!, ldiv!
import GeophysicalFlows.MultiLayerQG

# --
# Problems
# --

"""
    Problem(; parameters...)

Construct a constant diffusivity problem with steady or time-varying flow.
"""
noflow(args...) = 0.0 # used as defaults for u, v functions in Problem()

function Problem(dev;
          nx = 128,
          Lx = 2π,
          ny = nx,
          Ly = Lx,
           κ = 0.1,
           η = κ,
           u = noflow,
           v = noflow,
          dt = 0.01,
     stepper = "RK4",
  steadyflow = false,
           T = Float64
  )
  
  grid = TwoDGrid(dev, nx, Lx, ny, Ly; T=T)
  params = steadyflow==true ? ConstDiffSteadyFlowParams(η, κ, u, v, grid) : ConstDiffParams(η, κ, u, v)
  vars = Vars(dev, grid)
  equation = Equation(params, grid)

  return FourierFlows.Problem(equation, stepper, dt, grid, vars, params, dev)
end

"""
    Problem(dev, MQGprob; parameters...)

Construct a constant diffusivity problem on device `dev` using the flow from a
`GeophysicalFlows.MultiLayerQG` problem.
"""
function Problem(dev, MQGprob::FourierFlows.Problem;
                     κ = 0.1,
                     η = κ,
               stepper = "FilteredRK4",
   tracer_release_time = 0
  )
  
  grid = MQGprob.grid
  
  tracer_release_time < 0 && throw(ArgumentError("tracer_release_time must be non-negative!"))

  if tracer_release_time > 0
    @info "Stepping the flow forward until t = tracer_release_time = $tracer_release_time"
    step_until!(MQGprob, tracer_release_time)
  end

  params = TurbulentFlowParams(η, κ, tracer_release_time, MQGprob)
  vars = Vars(dev, grid, MQGprob)
  equation = Equation(params, grid)

  dt = MQGprob.clock.dt

  return FourierFlows.Problem(equation, stepper, dt, grid, vars, params, dev)
end

# --
# Params
# --

abstract type AbstractConstDiffParams <: AbstractParams end
abstract type AbstractSteadyFlowParams <: AbstractParams end
abstract type AbstractTurbulentFlowParams <: AbstractParams end

"""
    struct ConstDiffParams{T} <: AbstractConstDiffParams

A struct containing the parameters for a constant diffusivity problem with time-varying flow.
Included are:

$(TYPEDFIELDS)
"""
struct ConstDiffParams{T} <: AbstractConstDiffParams
    "isotropic horizontal diffusivity coefficient"
       η :: T
    "isotropic vertical diffusivity coefficient"
       κ :: T
    "isotropic hyperdiffusivity coefficient"
      κh :: T
    "isotropic hyperdiffusivity order"  
     nκh :: Int
    "function returning the x-component of advecting flow"
       u :: Function
    "function returning the y-component of advecting flow"
       v :: Function
end

"""
    ConstDiffParams(η, κ, u, v)

The constructor for the `params` struct for constant diffusivity problem and time-varying flow.
"""
ConstDiffParams(η, κ, u, v) = ConstDiffParams(η, κ, 0η, 0, u, v)

"""
    struct ConstDiffParams{T} <: AbstractConstDiffParams

A struct containing the parameters for a constant diffusivity problem with steady flow.
Included are:

$(TYPEDFIELDS)
"""
struct ConstDiffSteadyFlowParams{T, A} <: AbstractSteadyFlowParams
  "isotropic horizontal diffusivity coefficient"
     η :: T
  "isotropic vertical diffusivity coefficient"
     κ :: T
  "isotropic hyperdiffusivity coefficient"
    κh :: T
  "isotropic hyperdiffusivity order"  
   nκh :: Int
   "x-component of advecting flow"
     u :: A
   "y-component of advecting flow"
     v :: A
end

"""
    ConstDiffSteadyFlowParams(η, κ, κh, nκh, u::Function, v::Function, grid)
    ConstDiffSteadyFlowParams(η, κ, u, v, grid)

The constructor for the `params` struct for constant diffusivity problem and steady flow.
"""
function ConstDiffSteadyFlowParams(η, κ, κh, nκh, u::Function, v::Function, grid)
   x, y = gridpoints(grid)
  ugrid = u.(x, y)
  vgrid = v.(x, y)
  
  return ConstDiffSteadyFlowParams(η, κ, κh, nκh, ugrid, vgrid)
end

ConstDiffSteadyFlowParams(η, κ, u, v, grid) = ConstDiffSteadyFlowParams(η, κ, 0η, 0, u, v, grid)

"""
    struct TurbulentFlowParams{T} <: AbstractTurbulentFlowParams

A struct containing the parameters for a constant diffusivity problem with flow obtained
from a `GeophysicalFlows.MultiLayerQG` problem.

$(TYPEDFIELDS)
"""
struct TurbulentFlowParams{T} <: AbstractTurbulentFlowParams
  "isotropic horizontal diffusivity coefficient"
                    η :: T
  "isotropic vertical diffusivity coefficient"
                    κ :: T
  "isotropic hyperdiffusivity coefficient"
                   κh :: T
  "isotropic hyperdiffusivity order"  
                  nκh :: Int
  "number of layers in which the tracer is advected-diffused"
              nlayers :: Int 
  "flow time prior to releasing tracer"
  tracer_release_time :: T
  "`MultiLayerQG.Problem` to generate the advecting flow"
              MQGprob :: FourierFlows.Problem            
end

"""
    TurbulentFlowParams(η, κ, tracer_release_time, MQGprob)

The constructor for the `params` for a constant diffusivity problem with flow obtained
from a `GeophysicalFlows.MultiLayerQG` problem.
"""
function TurbulentFlowParams(η, κ, tracer_release_time, MQGprob)
  nlayers = numberoflayers(MQGprob.params)
  
  MultiLayerQG.updatevars!(MQGprob)

  return TurbulentFlowParams(η, κ, 0η, 0, nlayers, tracer_release_time, MQGprob)
end

# --
# Equations
# --

"""
    Equation(params, grid)

Return the equation for constant diffusivity problem with params p and grid g.
"""
function Equation(params::ConstDiffParams, grid)
  L = @. -params.η * grid.kr^2 - params.κ * grid.l^2 - params.κh * grid.Krsq^params.nκh
  
  return FourierFlows.Equation(L, calcN!, grid)
end

function Equation(params::ConstDiffSteadyFlowParams, grid)
  L = @. -params.η * grid.kr^2 - params.κ * grid.l^2 - params.κh * grid.Krsq^params.nκh
  
  return FourierFlows.Equation(L, calcN_steadyflow!, grid)
end

function Equation(params::TurbulentFlowParams, grid)
  L = zeros(grid.nkr, grid.nl, params.nlayers)

  for j in 1:params.nlayers
      @. L[:, :, j] = -params.η * grid.kr^2 - params.κ * grid.l^2 - params.κh * grid.Krsq^params.nκh
  end

  return FourierFlows.Equation(L, calcN_turbulentflow!, grid)
end

# --
# Vars
# --

"""
    struct Vars{Aphys, Atrans} <: AbstractVars

The variables for TracerAdvectionDiffussion problems.

$(FIELDS)
"""
struct Vars{Aphys, Atrans} <: AbstractVars
    "tracer concentration"
       c :: Aphys
    "tracer concentration x-derivative, ∂c/∂x"
      cx :: Aphys
    "tracer concentration y-derivative, ∂c/∂y"
      cy :: Aphys
    "Fourier transform of tracer concentration"
      ch :: Atrans
    "Fourier transform of tracer concentration x-derivative, ∂c/∂x"
     cxh :: Atrans
    "Fourier transform of tracer concentration y-derivative, ∂c/∂y"
     cyh :: Atrans
end

"""
    Vars(dev, grid)

Return the variables `vars` for a constant diffusivity problem on `grid` and device `dev`.
"""
function Vars(::Dev, grid::AbstractGrid{T}) where {Dev, T}

  @devzeros Dev T (grid.nx, grid.ny) c cx cy
  @devzeros Dev Complex{T} (grid.nkr, grid.nl) ch cxh cyh
  
  return Vars(c, cx, cy, ch, cxh, cyh)
end

function Vars(dev::Dev, grid::AbstractGrid{T}, MQGprob::FourierFlows.Problem) where {Dev, T}
  nlayers = numberoflayers(MQGprob.params)

  if nlayers == 1
    return Vars(dev, grid)
  else
    @devzeros Dev T (grid.nx, grid.ny, nlayers) c cx cy
    @devzeros Dev Complex{T} (grid.nkr, grid.nl, nlayers) ch cxh cyh
    
    return Vars(c, cx, cy, ch, cxh, cyh)
  end
end



# --
# Solvers
# --

"""
    calcN!(N, sol, t, clock, vars, params, grid)

Calculate the advective terms for a tracer equation with constant diffusivity and time-varying flow.
"""
function calcN!(N, sol, t, clock, vars, params::AbstractConstDiffParams, grid)
  @. vars.cxh = im * grid.kr * sol
  @. vars.cyh = im * grid.l  * sol

  ldiv!(vars.cx, grid.rfftplan, vars.cxh) # destroys vars.cxh when using fftw
  ldiv!(vars.cy, grid.rfftplan, vars.cyh) # destroys vars.cyh when using fftw

  x, y = gridpoints(grid)
  @. vars.cx = -params.u(x, y, clock.t) * vars.cx - params.v(x, y, clock.t) * vars.cy # copies over vars.cx so vars.cx = N in physical space
  mul!(N, grid.rfftplan, vars.cx)
  
  return nothing
end


"""
    calcN_steadyflow!(N, sol, t, clock, vars, params, grid)

Calculate the advective terms for a tracer equation with constant diffusivity and time-constant flow.
"""
function calcN_steadyflow!(N, sol, t, clock, vars, params::AbstractSteadyFlowParams, grid)
  @. vars.cxh = im * grid.kr * sol
  @. vars.cyh = im * grid.l  * sol

  ldiv!(vars.cx, grid.rfftplan, vars.cxh) # destroys vars.cxh when using fftw
  ldiv!(vars.cy, grid.rfftplan, vars.cyh) # destroys vars.cyh when using fftw

  @. vars.cx = -params.u * vars.cx - params.v * vars.cy # copies over vars.cx so vars.cx = N in physical space
  mul!(N, grid.rfftplan, vars.cx)
  
  return nothing
end

"""
    calcN_turbulentflow!(N, solt, t, clock, vars, params::AbstractTurbulentFlowParams, grid)

Calculate the advective terms for a tracer equation with constant diffusivity and turbulent `MultiLayerQG` flow.
"""
function calcN_turbulentflow!(N, sol, t, clock, vars, params::AbstractTurbulentFlowParams, grid)
  @. vars.cxh = im * grid.kr * sol
  @. vars.cyh = im * grid.l  * sol

  invtransform!(vars.cx, vars.cxh, params.MQGprob.params)
  invtransform!(vars.cy, vars.cyh, params.MQGprob.params)

  u = @. params.MQGprob.vars.u + params.MQGprob.params.U
  v = params.MQGprob.vars.v

  @. vars.cx = - u * vars.cx - v * vars.cy # copies over vars.cx so vars.cx = N in physical space
  fwdtransform!(N, vars.cx, params.MQGprob.params)

  return nothing
end

# --
# Helper functions
# --

"""
    updatevars!(prob)

Update the `prob.vars` in problem `prob` using the solution `prob.sol`.
"""
function updatevars!(params, vars, grid, sol)
  @. vars.ch = sol
  
  ldiv!(vars.c, grid.rfftplan, deepcopy(vars.ch))
  
  return nothing
end

"""
    updatevars!(params::AbstractTurbulentFlowParams, vars, grid, sol)

Update the `vars`` on the `grid` with the solution in `sol` for a problem `prob`
that is being advected by a turbulent flow.     
"""
function updatevars!(params::AbstractTurbulentFlowParams, vars, grid, sol)  
  @. vars.ch = sol
  
  invtransform!(vars.c, deepcopy(vars.ch), params.MQGprob.params)

  return nothing
end

updatevars!(prob) = updatevars!(prob.params, prob.vars, prob.grid, prob.sol)

"""
    set_c!(sol, params::Union{AbstractConstDiffParams, AbstractSteadyFlowParams}, grid, c)

Set the solution `sol` as the transform of `c` and update variables `vars`.
"""
function set_c!(sol, params::Union{AbstractConstDiffParams, AbstractSteadyFlowParams}, vars, grid, c)
  mul!(sol, grid.rfftplan, c)
  
  updatevars!(params, vars, grid, sol)
  
  return nothing
end

"""
    set_c!(sol, params::AbstractTurbulentFlowParams, grid, c)

Set the initial condition for tracer concentration in all layers of a
`TracerAdvectionDiffusion.Problem` that uses a `MultiLayerQG` flow to 
advect the tracer.
"""
function set_c!(sol, params::AbstractTurbulentFlowParams, vars, grid::AbstractGrid{T}, c) where T
  nlayers = numberoflayers(params.MQGprob.params)
  
  fwdtransform!(sol, repeat(c, 1, 1, nlayers), params.MQGprob.params)
  updatevars!(params, vars, grid, sol)

  return nothing
end

set_c!(prob, c) = set_c!(prob.sol, prob.params, prob.vars, prob.grid, c)

end # module
