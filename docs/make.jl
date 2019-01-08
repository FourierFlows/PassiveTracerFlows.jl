# Workaround for JuliaLang/julia/pull/28625
if Base.HOME_PROJECT[] !== nothing
  Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

using Documenter, PassiveTracerFlows

makedocs(
   modules = [PassiveTracerFlows],
   doctest = false,
     clean = true,
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

deploydocs( repo = "github.com/FourierFlows/PassiveTracerFlows.jl.git" )
