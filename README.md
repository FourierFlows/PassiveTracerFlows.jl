# PassiveTracerFlows.jl

 <p align="left">
     <a href="https://travis-ci.com/FourierFlows/PassiveTracerFlows.jl">
         <img alt="Build Status for CPU" src="https://img.shields.io/travis/com/FourierFlows/PassiveTracerFlows.jl/master?label=CPU&logo=travis&logoColor=white&style=flat-square">
     </a>
     <a href="https://gitlab.com/JuliaGPU/PassiveTracerFlows-jl/commits/master">
       <img alt="Build Status for GPU" src="https://img.shields.io/gitlab/pipeline/JuliaGPU/PassiveTracerFlows-jl/master?label=GPU&logo=gitlab&logoColor=white&style=flat-square">
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
