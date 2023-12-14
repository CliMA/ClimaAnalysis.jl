using Test

import JuliaFormatter
import ClimaAnalysis

@testset "Formatting" begin
    @test JuliaFormatter.format(
        ClimaAnalysis;
        verbose = false,
        overwrite = false,
    )
end
