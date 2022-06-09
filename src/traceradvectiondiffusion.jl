module TracerAdvectionDiffusion

export
  Problem,
  set_c!,
  updatevars!

using
  CUDA,
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

Construct a constant diffusivity problem with steady or time-varying flow in one or two
dimensions.

One dimension
"""
noflow(args...) = 0.0 # used as defaults for u, v functions in Problem()

function Problem(dev, onedim::Bool;
    nx = 128,
    Lx = 2π,
     κ = 0.1,
     u = noflow,
    dt = 0.01,
stepper = "RK4",
steadyflow = false,
     T = Float64
)

grid = OneDGrid(dev, nx, Lx; T=T)
params = steadyflow==true ? ConstDiffSteadyFlowParams(κ, u, grid::OneDGrid) : ConstDiffTimeVaryingFlowParams(κ, u)
vars = Vars(dev, grid; T=T)
equation = Equation(dev, params, grid)

return FourierFlows.Problem(equation, stepper, dt, grid, vars, params, dev)
end
"""
Two dimensions
"""
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
  params = steadyflow==true ? ConstDiffSteadyFlowParams(η, κ, u, v, grid::TwoDGrid) : ConstDiffTimeVaryingFlowParams(η, κ, u, v)
  vars = Vars(dev, grid; T=T)
  equation = Equation(dev, params, grid)

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

  params = ConstDiffTurbulentFlowParams(η, κ, tracer_release_time, MQGprob)
  vars = Vars(dev, grid, MQGprob)
  equation = Equation(dev, params, grid)

  dt = MQGprob.clock.dt

  return FourierFlows.Problem(equation, stepper, dt, grid, vars, params, dev)
end

# --
# Params
# --

abstract type AbstractTimeVaryingFlowParams <: AbstractParams end
abstract type AbstractSteadyFlowParams <: AbstractParams end
abstract type AbstractTurbulentFlowParams <: AbstractParams end

"""
    struct ConstDiffTimeVaryingFlowParams1D{T} <: AbstractTimeVaryingFlowParams

A struct containing the parameters for a constant diffusivity problem with time-varying flow in one dimension.
Included are:

$(TYPEDFIELDS)
"""
struct ConstDiffTimeVaryingFlowParams1D{T} <: AbstractTimeVaryingFlowParams
    "isotropic horizontal diffusivity coefficient"
       κ :: T
    "isotropic hyperdiffusivity coefficient"
      κh :: T
    "isotropic hyperdiffusivity order"  
     nκh :: Int
    "function returning the x-component of advecting flow"
       u :: Function
end

"""
    struct ConstDiffTimeVaryingFlowParams2D{T} <: AbstractTimeVaryingFlowParams

A struct containing the parameters for a constant diffusivity problem with time-varying flow in two dimensions.
Included are:

$(TYPEDFIELDS)
"""
struct ConstDiffTimeVaryingFlowParams2D{T} <: AbstractTimeVaryingFlowParams
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
    ConstDiffTimeVaryingFlowParams(κ, u)
    ConstDiffTimeVaryingFlowParams(η, κ, u, v)

The constructor for the `params` struct for constant diffusivity problem and time-varying flow.
"""
ConstDiffTimeVaryingFlowParams(κ, u) = ConstDiffTimeVaryingFlowParams1D(κ, 0κ, 0, u)
ConstDiffTimeVaryingFlowParams(η, κ, u, v) = ConstDiffTimeVaryingFlowParams2D(η, κ, 0η, 0, u, v)

"""
    struct ConstDiffSteadyFlowParams1D{T} <: AbstractSteadyFlowParams

A struct containing the parameters for a constant diffusivity problem with steady flow in one dimensions.
Included are:

$(TYPEDFIELDS)
"""
struct ConstDiffSteadyFlowParams1D{T, A} <: AbstractSteadyFlowParams
  "isotropic horizontal diffusivity coefficient"
     κ :: T
  "isotropic hyperdiffusivity coefficient"
    κh :: T
  "isotropic hyperdiffusivity order"  
   nκh :: Int
   "x-component of advecting flow"
     u :: A
end

"""
    struct ConstDiffSteadyFlowParams2D{T} <: AbstractSteadyFlowParams

A struct containing the parameters for a constant diffusivity problem with steady flow in two dimensions.
Included are:

$(TYPEDFIELDS)
"""
struct ConstDiffSteadyFlowParams2D{T, A} <: AbstractSteadyFlowParams
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
    ConstDiffSteadyFlowParams(κ, κh, nκh, u::Function, grid::OneDGrid)
    ConstDiffSteadyFlowParams(κ, u, grid::OneDGrid)
    ConstDiffSteadyFlowParams(η, κ, κh, nκh, u::Function, v::Function, grid::TwoDGrid)
    ConstDiffSteadyFlowParams(η, κ, u, v, grid::TwoDGrid)

The constructor for the `params` struct for constant diffusivity problem and steady flow.
"""
function ConstDiffSteadyFlowParams(κ, κh, nκh, u::Function, grid::OneDGrid)
    x = grid.x
   ugrid = u.(x)
   
   return ConstDiffSteadyFlowParams1D(κ, κh, nκh, ugrid)
 end
 
 ConstDiffSteadyFlowParams(κ, u, grid::OneDGrid) = ConstDiffSteadyFlowParams(κ, 0κ, 0, u, grid)

function ConstDiffSteadyFlowParams(η, κ, κh, nκh, u::Function, v::Function, grid::TwoDGrid)
   x, y = gridpoints(grid)
  ugrid = u.(x, y)
  vgrid = v.(x, y)
  
  return ConstDiffSteadyFlowParams2D(η, κ, κh, nκh, ugrid, vgrid)
end

ConstDiffSteadyFlowParams(η, κ, u, v, grid::OneDGrid) = ConstDiffSteadyFlowParams(η, κ, 0η, 0, u, v, grid)

"""
    struct ConstDiffTurbulentFlowParams{T} <: AbstractTurbulentFlowParams

A struct containing the parameters for a constant diffusivity problem with flow obtained
from a `GeophysicalFlows.MultiLayerQG` problem.

$(TYPEDFIELDS)
"""
struct ConstDiffTurbulentFlowParams{T} <: AbstractTurbulentFlowParams
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
    ConstDiffTurbulentFlowParams(η, κ, tracer_release_time, MQGprob)

The constructor for the `params` for a constant diffusivity problem with flow obtained
from a `GeophysicalFlows.MultiLayerQG` problem.
"""
function ConstDiffTurbulentFlowParams(η, κ, tracer_release_time, MQGprob)
  nlayers = numberoflayers(MQGprob.params)
  
  MultiLayerQG.updatevars!(MQGprob)

  return ConstDiffTurbulentFlowParams(η, κ, 0η, 0, nlayers, tracer_release_time, MQGprob)
end

# --
# Equations
# --

"""
    Equation(dev, params, grid)

Return the equation for constant diffusivity problem with `params` and `grid` on device `dev`.
"""
function Equation(dev, params::ConstDiffTimeVaryingFlowParams1D, grid)
    L = zeros(dev, eltype(grid), (grid.nkr))
    @. L = - params.κ * grid.kr^2 - params.κh * (grid.kr^2)^params.nκh
    
    return FourierFlows.Equation(L, calcN!, grid)
end

function Equation(dev, params::ConstDiffTimeVaryingFlowParams2D, grid)
  L = zeros(dev, eltype(grid), (grid.nkr, grid.nl))
  @. L = - params.η * grid.kr^2 - params.κ * grid.l^2 - params.κh * grid.Krsq^params.nκh
  
  return FourierFlows.Equation(L, calcN!, grid)
end

function Equation(dev, params::ConstDiffSteadyFlowParams1D, grid)
    L = zeros(dev, eltype(grid), (grid.nkr))
    @. L = - params.κ * grid.kr^2 - params.κh * (grid.kr^2)^params.nκh
    
    return FourierFlows.Equation(L, calcN_steadyflow!, grid)
end

function Equation(dev, params::ConstDiffSteadyFlowParams2D, grid)
  L = zeros(dev, eltype(grid), (grid.nkr, grid.nl))
  @. L = - params.η * grid.kr^2 - params.κ * grid.l^2 - params.κh * grid.Krsq^params.nκh
  
  return FourierFlows.Equation(L, calcN_steadyflow!, grid)
end

function Equation(dev, params::ConstDiffTurbulentFlowParams, grid)
  L = zeros(dev, eltype(grid), (grid.nkr, grid.nl, params.nlayers))

  for j in 1:params.nlayers
      @. L[:, :, j] = - params.η * grid.kr^2 - params.κ * grid.l^2 - params.κh * grid.Krsq^params.nκh
  end

  return FourierFlows.Equation(L, calcN_turbulentflow!, grid)
end

# --
# Vars
# --
"""
    struct Vars{Aphys, Atrans} <: AbstractVars

The variables for TracerAdvectionDiffussion problems in one dimension.

$(FIELDS)
"""
struct Vars1D{Aphys, Atrans} <: AbstractVars
    "tracer concentration"
       c :: Aphys
    "tracer concentration x-derivative, ∂c/∂x"
      cx :: Aphys
    "Fourier transform of tracer concentration"
      ch :: Atrans
    "Fourier transform of tracer concentration x-derivative, ∂c/∂x"
     cxh :: Atrans
end

"""
    struct Vars2D{Aphys, Atrans} <: AbstractVars

The variables for TracerAdvectionDiffussion problems in two dimensions.

$(FIELDS)
"""
struct Vars2D{Aphys, Atrans} <: AbstractVars
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
function Vars(::Dev, grid::OneDGrid; T=T) where {Dev}

    @devzeros Dev T (grid.nx) c cx
    @devzeros Dev Complex{T} (grid.nkr) ch cxh
    
    return Vars1D(c, cx, ch, cxh)
end

function Vars(::Dev, grid::TwoDGrid; T=T) where {Dev}

  @devzeros Dev T (grid.nx, grid.ny) c cx cy
  @devzeros Dev Complex{T} (grid.nkr, grid.nl) ch cxh cyh
  
  return Vars2D(c, cx, cy, ch, cxh, cyh)
end

function Vars(dev::Dev, grid::AbstractGrid{T}, MQGprob::FourierFlows.Problem) where {Dev, T}
  nlayers = numberoflayers(MQGprob.params)

  if nlayers == 1
    return Vars(dev, grid; T=T)
  else
    @devzeros Dev T (grid.nx, grid.ny, nlayers) c cx cy
    @devzeros Dev Complex{T} (grid.nkr, grid.nl, nlayers) ch cxh cyh
    
    return Vars2D(c, cx, cy, ch, cxh, cyh)
  end
end



# --
# Solvers
# --

"""
    calcN!(N, sol, t, clock, vars, params, grid)

Calculate the advective terms for a tracer equation with constant diffusivity and time-varying flow.
"""
function calcN!(N, sol, t, clock, vars, params::AbstractTimeVaryingFlowParams, grid::OneDGrid)
    @. vars.cxh = im * grid.kr * sol
  
    ldiv!(vars.cx, grid.rfftplan, vars.cxh) # destroys vars.cxh when using fftw
  
    x = grid.x
    @. vars.cx = -params.u(x, clock.t) * vars.cx # copies over vars.cx so vars.cx = N in physical space
    mul!(N, grid.rfftplan, vars.cx)
    
    return nothing
end

function calcN!(N, sol, t, clock, vars, params::AbstractTimeVaryingFlowParams, grid::TwoDGrid)
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
function calcN_steadyflow!(N, sol, t, clock, vars, params::AbstractSteadyFlowParams, grid::OneDGrid)
    @. vars.cxh = im * grid.kr * sol
  
    ldiv!(vars.cx, grid.rfftplan, vars.cxh) # destroys vars.cxh when using fftw
  
    @. vars.cx = -params.u * vars.cx # copies over vars.cx so vars.cx = N in physical space
    mul!(N, grid.rfftplan, vars.cx)
    
    return nothing
end

function calcN_steadyflow!(N, sol, t, clock, vars, params::AbstractSteadyFlowParams, grid::TwoDGrid)
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
    set_c!(sol, params::Union{AbstractTimeVaryingFlowParams, AbstractSteadyFlowParams}, grid, c)

Set the solution `sol` as the transform of `c` and update variables `vars`.
"""
function set_c!(sol, params::Union{AbstractTimeVaryingFlowParams, AbstractSteadyFlowParams}, vars, grid, c)
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
  
  C = @CUDA.allowscalar repeat(c, 1, 1, nlayers)
  fwdtransform!(sol, C, params.MQGprob.params)
  updatevars!(params, vars, grid, sol)

  return nothing
end

set_c!(prob, c) = set_c!(prob.sol, prob.params, prob.vars, prob.grid, c)

end # module
