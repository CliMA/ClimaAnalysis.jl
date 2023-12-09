using Test
import ClimaAnalysis
import CairoMakie

@testset "CairoMakieExt" begin
    @test length(methods(ClimaAnalysis.Visualize.contour_plot!)) == 1
end
