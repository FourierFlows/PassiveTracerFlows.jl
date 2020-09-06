module TracerAdvDiff

export
   Problem,
   set_c!,
   updatevars!

using
  FFTW,
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

function Problem(;
          nx = 128,
          Lx = 2π,
          ny = nx,
          Ly = Lx,
         kap = 0.1,
         eta = kap,
           u = noflow,
           v = noflow,
          dt = 0.01,
     stepper = "RK4",
  steadyflow = false,
           T = Float64,
         dev = CPU()
  )
  
  grid = TwoDGrid(dev, nx, Lx, ny, Ly; T=T)
  params = steadyflow==true ? ConstDiffSteadyFlowParams(eta, kap, u, v, grid) : ConstDiffParams(eta, kap, u, v)
  vars = Vars(dev, grid)
  equation = Equation(params, grid)

  FourierFlows.Problem(equation, stepper, dt, grid, vars, params, dev)
end


# --
# Params
# --

abstract type AbstractTracerParams <: AbstractParams end
abstract type AbstractConstDiffParams <: AbstractParams end
abstract type AbstractSteadyFlowParams <: AbstractParams end

"""
    ConstDiffParams(eta, kap, kaph, nkaph, u, v)
    ConstDiffParams(eta, kap, u, v)

Returns the params for constant diffusivity problem with time-varying flow.
"""
struct ConstDiffParams{T} <: AbstractConstDiffParams
  eta :: T           # Constant isotropic horizontal diffusivity
  kap :: T           # Constant isotropic vertical diffusivity
 kaph :: T           # Constant isotropic hyperdiffusivity
nkaph :: Int         # Constant isotropic hyperdiffusivity order
    u :: Function    # Advecting x-velocity
    v :: Function    # Advecting y-velocity
end
ConstDiffParams(eta, kap, u, v) = ConstDiffParams(eta, kap, 0eta, 0, u, v)

"""
    ConstDiffSteadyFlowParams(eta, kap, kaph, nkaph, u, v, g)
    ConstDiffSteadyFlowParams(eta, kap, u, v, g)

Returns the params for constant diffusivity problem with time-steady flow.
"""
struct ConstDiffSteadyFlowParams{T,A} <: AbstractSteadyFlowParams
  eta :: T           # Constant horizontal diffusivity
  kap :: T           # Constant vertical diffusivity
 kaph :: T           # Constant isotropic hyperdiffusivity
nkaph :: Int         # Constant isotropic hyperdiffusivity order
    u :: A           # Advecting x-velocity
    v :: A           # Advecting y-velocity
end

function ConstDiffSteadyFlowParams(eta, kap, kaph, nkaph, u::Function, v::Function, grid)
   x, y = gridpoints(grid)
  ugrid = u.(x, y)
  vgrid = v.(x, y)
  ConstDiffSteadyFlowParams(eta, kap, kaph, nkaph, ugrid, vgrid)
end

ConstDiffSteadyFlowParams(eta, kap, u, v, grid) = ConstDiffSteadyFlowParams(eta, kap, 0eta, 0, u, v, grid)


# --
# Equations
# --

"""
    Equation(p, g)

Returns the equation for constant diffusivity problem with params p and grid g.
"""
function Equation(params::ConstDiffParams, grid)
  L = @. -params.eta * grid.kr^2 - params.kap * grid.l^2 - params.kaph * grid.Krsq^params.nkaph
  FourierFlows.Equation(L, calcN!, grid)
end

function Equation(params::ConstDiffSteadyFlowParams, grid)
  L = @. -params.eta * grid.kr^2 - params.kap * grid.l^2 - params.kaph * grid.Krsq^params.nkaph
  FourierFlows.Equation(L, calcN_steadyflow!, grid)
end


# --
# Vars
# --

struct Vars{Aphys, Atrans} <: AbstractVars
    c :: Aphys
   cx :: Aphys
   cy :: Aphys
   ch :: Atrans
  cxh :: Atrans
  cyh :: Atrans
end

"""
    Vars(g)

Returns the vars for constant diffusivity problem on grid g.
"""
function Vars(::Dev, grid::AbstractGrid{T}) where {Dev, T}
  @devzeros Dev T (grid.nx, grid.ny) c cx cy
  @devzeros Dev Complex{T} (grid.nkr, grid.nl) ch cxh cyh
  Vars(c, cx, cy, ch, cxh, cyh)
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
