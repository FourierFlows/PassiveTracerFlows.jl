# TracerAdvectionDiffusion Module

### Basic Equations

This module solves the advection-diffusion equation for a passive tracer concentration in
1D or 2D domains. 

For 1D problems the tracer concentration ``c(x, t)`` evolves under:

```math
\partial_t c + u \partial_x c = \underbrace{\kappa \partial_x^2 c}_{\textrm{diffusivity}} + \underbrace{\kappa_h (-1)^{n_{h}} \partial_x^{2n_{h}}c}_{\textrm{hyper-diffusivity}}\ ,
```

where ``u(x, t)`` is the advecting flow and ``\kappa`` the diffusivity. The advecting flow could be either compressible or incompressible. 

For 2D problems the tracer concentration ``c(x, y, t)`` evolves under:

```math
\partial_t c + \bm{u} \bm{\cdot} \bm{\nabla} c = \underbrace{\eta \partial_x^2 c + \kappa \partial_y^2 c}_{\textrm{diffusivity}} + \underbrace{\kappa_h (-1)^{n_{h}} \nabla^{2n_{h}}c}_{\textrm{hyper-diffusivity}}\ ,
```

where ``\bm{u} = (u, v)`` is the two-dimensional advecting flow, ``\eta`` the ``x``-diffusivity and ``\kappa`` is the ``y``-diffusivity. If ``\eta`` is not defined then the code uses isotropic diffusivity, i.e., ``\eta \partial_x^2 c + \kappa \partial_y^2 c \mapsto \kappa \nabla^2``. The advecting flow could be either compressible or incompressible. 


### Implementation

The equation is time-stepped forward in Fourier space. For example, for 2D problems:

```math
\partial_t \widehat{c} = - \widehat{\bm{u} \bm{\cdot} \bm{\nabla} c} - \left[ (\eta k_x^2 + \kappa k_y^2) + \kappa_h |\bm{k}|^{2\nu_h} \right] \widehat{c}\ ,
```
where ``\bm{k} = (k_x, k_y)``.

Thus:

```math
\begin{aligned}
L & = -\eta k_x^2 - \kappa k_y^2 - \kappa_h |\bm{k}|^{2\nu_h} , \\
N(\widehat{c}) &= - \mathrm{FFT}(u \partial_x c + v \partial_y c) .
\end{aligned}
```
