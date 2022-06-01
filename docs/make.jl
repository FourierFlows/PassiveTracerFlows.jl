using
  Documenter,
  Literate,
  Plots,   # to not capture precompilation output
  PassiveTracerFlows

# Gotta set this environment variable when using the GR run-time on Travis CI.
# This happens as examples will use Plots.jl to make plots and movies.
# See: https://github.com/jheinen/GR.jl/issues/278
ENV["GKSwstype"] = "100"

#####
##### Generate examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")

examples = [
    "cellularflow.jl",
    "turbulent_advection-diffusion.jl"
]

for example in examples
  withenv("GITHUB_REPOSITORY" => "FourierFlows/PassiveTracerFlowsDocumentation") do
    example_filepath = joinpath(EXAMPLES_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor())
    Literate.notebook(example_filepath, OUTPUT_DIR)
    Literate.script(example_filepath, OUTPUT_DIR)
  end
end

#####
##### Build and deploy docs
#####

# Set up a timer to print a space ' ' every 240 seconds. This is to avoid Travis CI
# timing out when building demanding Literate.jl examples.
Timer(t -> println(" "), 0, interval=240)

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
      canonical = "https://fourierflows.github.io/PassiveTracerFlowsDocumentation/stable/"
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
           "Examples" => Any[
             "TracerAdvectionDiffusion" => Any[
               "literated/cellularflow.md",
               "literated/turbulent_advection-diffusion.md",
               ]
           ],
           "Modules" => Any[
             "modules/traceradvectiondiffusion.md",
           ],
           "DocStrings" => Any[
             "man/types.md",
             "man/functions.md"
             ]
           ]
)

withenv("GITHUB_REPOSITORY" => "FourierFlows/PassiveTracerFlowsDocumentation") do
  deploydocs(       repo = "github.com/FourierFlows/PassiveTracerFlowsDocumentation.git",
                versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
            push_preview = false,
               devbranch = "main"
            )
end
