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
        Base.get_extension(ClimaAnalysis, :ClimaAnalysisMakieExt),
        Base.get_extension(ClimaAnalysis, :ClimaAnalysisGeoMakieExt),
    ],
    authors = "Climate Modelling Alliance",
    sitename = "ClimaAnalysis.jl",
    format = Documenter.HTML(;
        prettyurls = !isempty(get(ENV, "CI", "")),
        collapselevel = 1,
        size_threshold_ignore = ["api.md"],
    ),
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "OutputVars" => "var.md",
        "Visualizing OutputVars" => "visualize.md",
        "RMSEVariables" => "rmse_var.md",
        "Visualizing RMSEVariables" => "visualize_rmse_var.md",
        "FlatVar" => "flat.md",
        "APIs" => "api.md",
        "How do I?" => "howdoi.md",
        "Developer Documentation" => "developer.md",
    ],
)

deploydocs(; repo = "github.com/CliMA/ClimaAnalysis.jl", push_preview = true)
