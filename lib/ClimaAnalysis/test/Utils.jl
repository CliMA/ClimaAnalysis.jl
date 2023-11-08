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
