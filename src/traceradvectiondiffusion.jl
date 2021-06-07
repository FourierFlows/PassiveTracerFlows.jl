module TracerAdvectionDiffusion

export
   Problem,
   set_c!,
   updatevars!

using
   DocStringExtensions,
   Reexport

@reexport using FourierFlows

import LinearAlgebra: mul!, ldiv!

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


# --
# Params
# --

abstract type AbstractTracerParams <: AbstractParams end
abstract type AbstractConstDiffParams <: AbstractParams end
abstract type AbstractSteadyFlowParams <: AbstractParams end

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
   "x-component of advecting flow"
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


# --
# Helper functions
# --

"""
    updatevars!(prob)

Update the vars in v on the grid g with the solution in sol.
"""
function updatevars!(prob)
  vars, grid, sol = prob.vars, prob.grid, prob.sol
  
  @. vars.ch = sol
  
  ldiv!(vars.c, grid.rfftplan, deepcopy(vars.ch))
  
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

end # module
