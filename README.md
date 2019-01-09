# PassiveTracerFlows.jl

 <p align="left">
     <a href="https://travis-ci.org/FourierFlows/PassiveTracerFlows.jl">
          <img src="https://travis-ci.org/FourierFlows/PassiveTracerFlows.jl.svg?branch=master" title="Build Status">
     </a>
     <a href="https://ci.appveyor.com/project/navidcy/passivetracerflows-jl">
          <img src="https://ci.appveyor.com/api/projects/status/yb5ywk9bof6u3nyg?svg=true" title="Build Status">
     </a>
     <a href="https://fourierflows.github.io/PassiveTracerFlows.jl/stable/">
         <img src="https://img.shields.io/badge/docs-stable-blue.svg">
     </a>
     <a href="https://fourierflows.github.io/PassiveTracerFlows.jl/latest/">
         <img src="https://img.shields.io/badge/docs-latest-blue.svg">
     </a>
     <a href='https://coveralls.io/github/FourierFlows/PassiveTracerFlows.jl?branch=master'><img src='https://coveralls.io/repos/github/FourierFlows/PassiveTracerFlows.jl/badge.svg?branch=master' alt='Coverage Status' />
     </a>
     <a href="https://doi.org/10.5281/zenodo.2535984">
        <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.2535984.svg" alt="DOI">
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

 * `TracerAdvDiff`: advection-diffusion of a passive tracer in 2D domain.


 [FourierFlows.jl]: https://github.com/FourierFlows/FourierFlows.jl
