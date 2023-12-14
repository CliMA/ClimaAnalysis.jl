using ClimaAnalysis
using Documenter

makedocs(;
    modules = [ClimaAnalysis],
    authors = "Climate Modelling Alliance",
    repo = "https://github.com/Sbozzolo/ClimaAnalysis.jl",
    sitename = "ClimaAnalysis.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://sbozzolo.github.io/ClimaAnalysis.jl",
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/Sbozzolo/ClimaAnalysis.jl")
