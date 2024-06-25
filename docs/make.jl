using ClimaAnalysis
using Documenter
import CairoMakie
import GeoMakie

DocMeta.setdocmeta!(
    ClimaAnalysis,
    :DocTestSetup,
    :(using ClimaAnalysis;
    using ClimaAnalysis.Utils;
    using ClimaAnalysis.Var;
    using ClimaAnalysis.Sim);
    recursive = true,
)

makedocs(;
    modules = [
        ClimaAnalysis,
        Base.get_extension(ClimaAnalysis, :CairoMakieExt),
        Base.get_extension(ClimaAnalysis, :GeoMakieExt),
    ],
    authors = "Climate Modelling Alliance",
    repo = "https://github.com/CliMA/ClimaAnalysis.jl",
    sitename = "ClimaAnalysis.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://CliMA.github.io/ClimaAnalysis.jl",
    ),
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "APIs" => "api.md",
        "How do I?" => "howdoi.md",
    ],
)

deploydocs(; repo = "github.com/CliMA/ClimaAnalysis.jl")
