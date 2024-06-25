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

@testset "seconds_to_prettystr" begin
    @test Utils.seconds_to_prettystr(1) == "1s"
    @test Utils.seconds_to_prettystr(2) == "2s"
    @test Utils.seconds_to_prettystr(60) == "1m"
    @test Utils.seconds_to_prettystr(61) == "1m 1s"
    @test Utils.seconds_to_prettystr(3600) == "1h"
    @test Utils.seconds_to_prettystr(3666) == "1h 1m 6s"
    @test Utils.seconds_to_prettystr(86400) == "1d"
    @test Utils.seconds_to_prettystr(86401) == "1d 1s"
    @test Utils.seconds_to_prettystr(86522) == "1d 2m 2s"
    @test Utils.seconds_to_prettystr(24 * 60 * 60 * 365) == "1y"
    @test Utils.seconds_to_prettystr(24 * 60 * 60 * 365 + 60) == "1y 1m"
    @test Utils.seconds_to_prettystr(24 * 60 * 60 * 365 + 68) == "1y 1m 8s"

    @test Utils.seconds_to_prettystr(24 * 60 * 60 * 365 + 68.5) == "1y 1m 8.5s"
end
