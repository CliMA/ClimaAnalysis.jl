using Test
import ClimaAnalysis: Utils

@testset "Regexp" begin

    @test Utils.match_nc_filename("bob") === nothing

    @test Utils.match_nc_filename("ta_1d_average.nc") ==
          Tuple(["ta", "1d", "average"])

    @test Utils.match_nc_filename("pfull_6.0min_max.nc") ==
          Tuple(["pfull", "6.0min", "max"])

    @test Utils.match_nc_filename("hu_inst.nc") ==
          Tuple(["hu", nothing, "inst"])
end

@testset "Squeeze" begin
    @test Utils.squeeze([[1 2] [3 4]]) == [1, 2, 3, 4]
end

@testset "nearest_index" begin
    @test Utils.nearest_index([10, 20, 30], 0) == 1
    @test Utils.nearest_index([10, 20, 30], 100) == 3
    @test Utils.nearest_index([10, 20, 30], 20) == 2
    @test Utils.nearest_index([10, 20, 30], typemax(Int)) == 3
    @test Utils.nearest_index([10, 20, 30], -typemax(Int)) == 1
end

@testset "kwargs" begin
    @test values(Utils.kwargs(a = 1, b = 2)) == (; a = 1, b = 2)
end
