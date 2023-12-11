using Test
using ClimaAnalysis
using Aqua

@testset "Aqua tests" begin
    Aqua.test_undefined_exports(ClimaAnalysis)
    Aqua.test_stale_deps(ClimaAnalysis)
    Aqua.test_deps_compat(ClimaAnalysis)
    Aqua.detect_ambiguities(ClimaAnalysis; recursive = true)
    Aqua.test_piracies(ClimaAnalysis)
end
