using
  Documenter,
  PassiveTracerFlows

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
      canonical = "https://fourierflows.github.io/PassiveTracerFlowsDocumentation/dev/"
)

makedocs(
   modules = [PassiveTracerFlows],
   doctest = false,
     clean = true,
 checkdocs = :all,
    format = format,
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

withenv("GITHUB_REPOSITORY" => "FourierFlows/PassiveTracerFlowsDocumentation") do
  deploydocs(        repo = "github.com/FourierFlows/PassiveTracerFlowsDocumentation.git",
                versions = ["stable" => "v^", "v#.#"],
            push_preview = true
            )
end
