# TracerAdvectionDiffusion Module

### Basic Equations

This module solves the advection-diffusion equation for a passive tracer concentration in
1D, 2D, or 3D domains. 

For 1D problems the tracer concentration ``c(x, t)`` evolves under:

```math
\partial_t c + u \partial_x c = \underbrace{\kappa \partial_x^2 c}_{\textrm{diffusivity}} + \underbrace{\kappa_h (-1)^{n_h + 1} \partial_x^{2 n_h} c}_{\textrm{hyper-diffusivity}} \ ,
```

where ``u(x, t)`` is the advecting flow and ``\kappa`` the diffusivity, ``\kappa_h`` is the hyper-diffusivity
coefficient and ``n_h`` the hyper-diffusivity order.

For 2D problems, ``\boldsymbol{x} = (x, y)``, the tracer concentration ``c(\boldsymbol{x}, t)`` evolves under:

```math
\partial_t c + \bm{u} \bm{\cdot} \bm{\nabla} c = \underbrace{\kappa \partial_x^2 c + \eta \partial_y^2 c}_{\textrm{diffusivity}} + \underbrace{\kappa_h (-1)^{n_h + 1} \nabla^{2 n_h} c}_{\textrm{hyper-diffusivity}} \ ,
```

where ``\bm{u} = (u, v)`` is the two-dimensional advecting flow, ``\kappa`` the ``x``-diffusivity and ``\eta``
is the ``y``-diffusivity. If ``\eta`` is not defined then by default it is set to have the same value as
``\kappa``. See [`TracerAdvectionDiffusion.Problem`]

For 3D problems, ``\boldsymbol{x} = (x, y, z)``, the tracer concentration ``c(\boldsymbol{x}, t)`` evolves via:

```math
\partial_t c + \bm{u} \bm{\cdot} \bm{\nabla} c = \underbrace{\kappa \partial_x^2 c + \eta \partial_y^2 c + \ell \partial_z^2 c}_{\textrm{diffusivity}} + \underbrace{\kappa_h (-1)^{n_h + 1} \nabla^{2 n_h} c}_{\textrm{hyper-diffusivity}} \ ,
```

where ``\bm{u} = (u, v, w)`` is the three-dimensional advecting flow, ``\kappa`` the ``x``-diffusivity,
``\eta`` is the ``y``-diffusivity, and ``\ell`` the ``z``-diffusivity. If ``\eta`` or ``\ell`` are not
defined then by default are set to have the same value as ``\kappa``.

The advecting flow can be either compressible or incompressible. 


### Implementation

The equations are time-stepped forward in Fourier space. For example, for 2D problems:

```math
\partial_t \widehat{c} = - \widehat{\bm{u} \bm{\cdot} \bm{\nabla} c} - \left ( \kappa k_x^2 + \eta k_y^2 + \kappa_h |\bm{k}|^{2 n_h} \right) \widehat{c} \ ,
```
where ``\bm{k} = (k_x, k_y)``.

Thus:

```math
\begin{aligned}
L & = - \kappa k_x^2 - \eta k_y^2 - \kappa_h |\bm{k}|^{2 n_h} \ , \\
N(\widehat{c}) &= - \mathrm{FFT}(u \partial_x c + v \partial_y c) \ .
\end{aligned}
```
