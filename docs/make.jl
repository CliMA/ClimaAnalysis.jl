using ClimaAnalysis
using Documenter
import CairoMakie

DocMeta.setdocmeta!(
    ClimaAnalysis,
    :DocTestSetup,
    :(using ClimaAnalysis.Utils);
    recursive = true,
)

makedocs(;
    modules = [
        ClimaAnalysis,
        Base.get_extension(ClimaAnalysis, :CairoMakieExt),
    ],
    authors = "Climate Modelling Alliance",
    repo = "https://github.com/CliMA/ClimaAnalysis.jl",
    sitename = "ClimaAnalysis.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://CliMA.github.io/ClimaAnalysis.jl",
    ),
    checkdocs = :exports,
    pages = ["Home" => "index.md", "APIs" => "api.md"],
)

deploydocs(; repo = "github.com/CliMA/ClimaAnalysis.jl")
