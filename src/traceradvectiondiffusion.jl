module TracerAdvectionDiffusion

export
   Problem,
   set_c!,
   updatevars!

using
   DocStringExtensions,
   Reexport

@reexport using FourierFlows, GeophysicalFlows.MultiLayerQG

using GeophysicalFlows.MultiLayerQG: SingleLayerParams, TwoLayerParams

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
    Problem(MQGprob; parameters...)

Construct a constant diffusivity problem with a turbulent flow from the `GeophysicalFlows.jl` package.
"""
function Problem(dev, MQGprob::FourierFlows.Problem;
                     κ = 0.1,
                     η = κ,
               stepper = "FilteredRK4",
   tracer_release_time = 0
  )
  
  nlayers = typeof(MQGprob.params) <: SingleLayerParams ? 1 : 
            typeof(MQGprob.params) <: TwoLayerParams ? 2 : MQGprob.params.nlayers

  grid = MQGprob.grid
  
  if tracer_release_time > 0
    @info "Stepping the flow forward until `t = tracer_release_time`"
    step_until!(MQGprob, tracer_release_time)
  end

  params = TurbulentFlowParams(η, κ, nlayers, tracer_release_time, MQGprob)
  vars = Vars(dev, grid, nlayers)
  equation = Equation(params, grid)

  dt = MQGprob.clock.dt

  return FourierFlows.Problem(equation, stepper, dt, grid, vars, params, dev)
end

# --
# Params
# --

abstract type AbstractTracerParams <: AbstractParams end
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
struct ConstDiffSteadyFlowParams{T,A} <: AbstractSteadyFlowParams
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

A struct containing the parameters for a constant diffusivity problem and turbulent 
`MultiLayerQG` flow.

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
    TurbulentFlowParams(η, κ, nlayers, tracer_release_time)

The constructor for the `params` struct for a constant diffusivity and turbulent flow.    
"""
TurbulentFlowParams(η, κ, nlayers, tracer_release_time, MQGprob) =
  TurbulentFlowParams(η, κ, 0η, 0, nlayers, tracer_release_time, MQGprob)

# --
# Equations
# --

"""
    Equation(params, grid)

Returns the equation for constant diffusivity problem with params p and grid g.
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

    for n in 1:params.nlayers
        @. L[:, :, n] = -params.η * grid.kr^2 - params.κ * grid.l^2 - params.κh * grid.Krsq^params.nκh
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

Returns the variables `vars` for a constant diffusivity problem on `grid` and device `dev`.
"""
function Vars(::Dev, grid::AbstractGrid{T}) where {Dev, T}

  @devzeros Dev T (grid.nx, grid.ny) c cx cy
  @devzeros Dev Complex{T} (grid.nkr, grid.nl) ch cxh cyh
  
  return Vars(c, cx, cy, ch, cxh, cyh)
end

function Vars(::Dev, grid::AbstractGrid{T}, nlayers::Int) where {Dev, T}
  if nlayers == 1
    @devzeros Dev T (grid.nx, grid.ny) c cx cy
    @devzeros Dev Complex{T} (grid.nkr, grid.nl) ch cxh cyh
  else
    @devzeros Dev T (grid.nx, grid.ny, nlayers) c cx cy
    @devzeros Dev Complex{T} (grid.nkr, grid.nl, nlayers) ch cxh cyh
  end
  
  return Vars(c, cx, cy, ch, cxh, cyh)
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
  
    MultiLayerQG.invtransform!(vars.cx, vars.cxh, params.MQGprob.params)
    MultiLayerQG.invtransform!(vars.cy, vars.cyh, params.MQGprob.params)
  
    u = @. params.MQGprob.vars.u + params.MQGprob.params.U
    v = params.MQGprob.vars.u

    # Step the flow forward for next iteration
    MultiLayerQG.stepforward!(params.MQGprob)
    MultiLayerQG.updatevars!(params.MQGprob)
    @. vars.cx = -u * vars.cx - v * vars.cy # copies over vars.cx so vars.cx = N in physical space

    MultiLayerQG.fwdtransform!(N, vars.cx, params.MQGprob.params)
    
    return nothing
  end

# --
# Helper functions
# --

"""
    updatevars!(prob)

Update the `vars` on the `grid` with the solution in `sol`.
"""
function updatevars!(prob)
  vars, grid, sol = prob.vars, prob.grid, prob.sol
  
  @. vars.ch = sol
  
  ldiv!(vars.c, grid.rfftplan, deepcopy(vars.ch))
  
  return nothing
end

"""
    MQGupdatevars!(prob)

Update the `vars`` on the `grid` with the solution in `sol` for a problem `prob`
that is being advected by a turbulent flow.     
"""
function MQGupdatevars!(prob)
  vars, grid, sol = prob.vars, prob.grid, prob.sol
  
  @. vars.ch = sol
  
  MultiLayerQG.invtransform!(vars.c, deepcopy(vars.ch), prob.params.MQGprob.params)

  return nothing
end

"""
    set_c!(prob, c)

Set the solution sol as the transform of c and update variables v
on the grid g.
"""
function set_c!(prob, c)
  sol, vars, grid = prob.sol, prob.vars, prob.grid

  mul!(sol, grid.rfftplan, c)
  
  updatevars!(prob)
  
  return nothing
end

"""
    function set_c!(prob, c, nlayers)

Set the initial condition for tracer concentration in all layers of a
`TracerAdvectionDiffusion.Problem` that uses a `MultiLayerQG` flow to 
advect the tracer.
"""
function set_c!(prob, c, nlayers)
  sol, vars, grid = prob.sol, prob.vars, prob.grid

  C = Array{Float64}(undef, grid.nx, grid.ny, nlayers)

  if size(c) == size(C)
    @. C = c
  else
    for n in 1:nlayers
        C[:, :, n] = c
    end
  end
  
  MultiLayerQG.fwdtransform!(sol, C, prob.params.MQGprob.params)
  MQGupdatevars!(prob)

  return nothing
end

end # module
