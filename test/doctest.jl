using Documenter
import ClimaAnalysis

@testset "Test docstrings" begin

    DocMeta.setdocmeta!(
        ClimaAnalysis,
        :DocTestSetup,
        :(using ClimaAnalysis;
        using ClimaAnalysis.Utils;
        using ClimaAnalysis.Var;
        using ClimaAnalysis.Sim);
        recursive = true,
    )

    doctest(ClimaAnalysis; manual = false)
end
