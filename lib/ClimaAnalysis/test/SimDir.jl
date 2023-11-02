using Test
import ClimaAnalysis

@testset "SimDir" begin

    simulation_path = joinpath(@__DIR__, "sample_data")

    simdir = ClimaAnalysis.SimDir(simulation_path)

    @test ClimaAnalysis.available_vars(simdir) ==
          Set(["va", "ua", "orog", "ta", "ts"])

    @test simdir.variables_paths["orog"]["inst"] ==
          joinpath(simulation_path, "orog_inst.nc")
    @test simdir.variables_paths["ta"]["max"] == Dict{Any, Any}(
        "3.0h" => joinpath(simulation_path, "ta_3.0h_max.nc"),
        "4.0h" => joinpath(simulation_path, "ta_4.0h_max.nc"),
    )

    expected_files = Set(
        joinpath(simulation_path, f) for
        f in readdir(simulation_path) if endswith(f, ".nc")
    )

    @test simdir.allfiles == expected_files
end
