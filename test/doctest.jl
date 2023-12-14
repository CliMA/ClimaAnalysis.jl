using Documenter
import ClimaAnalysis

@testset "Test docstrings" begin

    DocMeta.setdocmeta!(
        ClimaAnalysis,
        :DocTestSetup,
        :(using ClimaAnalysis.Utils);
        recursive = true,
    )

    doctest(ClimaAnalysis; manual = false)
end
