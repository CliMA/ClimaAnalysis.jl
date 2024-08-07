using ClimaAnalysis
using Documenter
import GeoMakie

DocMeta.setdocmeta!(
    ClimaAnalysis,
    :DocTestSetup,
    :(using ClimaAnalysis;
    using ClimaAnalysis.Utils;
    using ClimaAnalysis.Var;
    using ClimaAnalysis.Atmos;
    using ClimaAnalysis.Sim);
    recursive = true,
)

makedocs(;
    modules = [
        ClimaAnalysis,
        Base.get_extension(ClimaAnalysis, :MakieExt),
        Base.get_extension(ClimaAnalysis, :GeoMakieExt),
    ],
    authors = "Climate Modelling Alliance",
    sitename = "ClimaAnalysis.jl",
    format = Documenter.HTML(;
        prettyurls = !isempty(get(ENV, "CI", "")),
        collapselevel = 1,
    ),
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "APIs" => "api.md",
        "How do I?" => "howdoi.md",
    ],
)

deploydocs(; repo = "github.com/CliMA/ClimaAnalysis.jl")
