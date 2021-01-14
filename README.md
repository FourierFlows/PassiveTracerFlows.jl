# PassiveTracerFlows.jl

 <p align="left">
 <p align="center">
     <a href="https://buildkite.com/julialang/passivetracerflows-dot-jl">
         <img alt="Buildkite CPU+GPU build status" src="https://img.shields.io/buildkite/4d921fc17b95341ea5477fb62df0e6d9364b61b154e050a123/master?logo=buildkite&label=Buildkite%20CPU%2BGPU">
     </a>
     <a href="https://ci.appveyor.com/project/navidcy/passivetracerflows-jl">
         <img alt="Build Status for Window" src="https://img.shields.io/appveyor/ci/navidcy/passivetracerflows-jl/master?label=Window&logo=appveyor&logoColor=white&style=flat-square">
     </a>
     <a href="https://fourierflows.github.io/PassiveTracerFlowsDocumentation/stable/">
         <img src="https://img.shields.io/badge/docs-stable-blue.svg">
     </a>
     <a href="https://fourierflows.github.io/PassiveTracerFlowsDocumentation/latest/">
         <img src="https://img.shields.io/badge/docs-dev-blue.svg">
     </a>
     <a href="https://codecov.io/gh/FourierFlows/PassiveTracerFlows.jl">
         <img src="https://codecov.io/gh/FourierFlows/PassiveTracerFlows.jl/branch/master/graph/badge.svg" title="codecov">
     </a>
     <a href="https://doi.org/10.5281/zenodo.2535983">
        <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.2535983.svg" alt="DOI">
    </a>
 </p>

This package leverages the [FourierFlows.jl]() framework to provide modules for solving passive tracer advection-diffusion problems on periodic domains using Fourier-based pseudospectral methods.

 ## Installation

 To install, do
 ```julia
 ] add PassiveTracerFlows
 ```

 See `examples/` for example scripts.

 ## Modules

 * `TracerAdvectionDiffusion`: advection-diffusion of a passive tracer in a 2D domain.


 ## Cite

 The code is citable via [zenodo](https://zenodo.org). Please cite as:

 > Navid C. Constantinou. and Gregory L. Wagner (2021). FourierFlows/PassiveTracerFlows.jl: PassiveTracerFlows v0.4.1 (Version v0.4.1). Zenodo.  [https://doi.org/10.5281/zenodo.2535983](https://doi.org/10.5281/zenodo.2535983)

 [FourierFlows.jl]: https://github.com/FourierFlows/FourierFlows.jl
