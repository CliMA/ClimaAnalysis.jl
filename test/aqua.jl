using Test
using ClimaAnalysis
using Aqua
import ExplicitImports

@testset "Aqua tests" begin
    Aqua.test_undefined_exports(ClimaAnalysis)
    Aqua.test_stale_deps(ClimaAnalysis)
    Aqua.test_deps_compat(ClimaAnalysis)
    Aqua.detect_ambiguities(ClimaAnalysis; recursive = true)
    Aqua.test_piracies(ClimaAnalysis)
end

@testset "Explicit Imports" begin
    @test isnothing(ExplicitImports.check_no_implicit_imports(ClimaAnalysis))
    @test isnothing(
        ExplicitImports.check_no_stale_explicit_imports(ClimaAnalysis),
    )
    @test isnothing(
        ExplicitImports.check_all_qualified_accesses_via_owners(ClimaAnalysis),
    )
    @test isnothing(
        ExplicitImports.check_no_self_qualified_accesses(ClimaAnalysis),
    )
end
