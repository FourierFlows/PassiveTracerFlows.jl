using Documenter, FourierFlows

makedocs(
   modules = [PassiveTracerFlows],
   doctest = false, clean = true,
 checkdocs = :all,
    format = :html,
   authors = "Gregory L. Wagner and Navid C. Constantinou",
  sitename = "PassiveTracerFlows.jl",
     pages = Any[
              "Home" => "index.md",
              "Modules" => Any[
                "modules/traceradvdiff.md",
              ],
              "DocStrings" => Any[
              "man/types.md",
              "man/functions.md"]
             ]
)

deploydocs(
       repo = "github.com/FourierFlows/PassiveTracerFlows.jl.git",
     target = "build",
      julia = "0.7",
     osname = "linux",
       deps = nothing,
       make = nothing
 )
