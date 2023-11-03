using Test
import ClimaAnalysis

@testset "SimDir" begin

    simulation_path = joinpath(@__DIR__, "sample_data")

    simdir = ClimaAnalysis.SimDir(simulation_path)

    @test ClimaAnalysis.available_vars(simdir) ==
          Set(["va", "ua", "orog", "ta", "ts"])

    @test ClimaAnalysis.available_reductions(simdir, short_name = "ta") ==
          Set(["average", "max", "min"])
    @test_throws ErrorException ClimaAnalysis.available_reductions(
        simdir,
        short_name = "bob",
    )

    @test ClimaAnalysis.available_periods(
        simdir,
        short_name = "ta",
        reduction = "max",
    ) == Set(["3.0h", "4.0h"])
    @test_throws ErrorException ClimaAnalysis.available_periods(
        simdir,
        short_name = "ta",
        reduction = "bob",
    )

    @test simdir.variable_paths["orog"]["inst"][nothing] ==
          joinpath(simulation_path, "orog_inst.nc")
    @test simdir.variable_paths["ta"]["max"] == Dict{Any, Any}(
        "3.0h" => joinpath(simulation_path, "ta_3.0h_max.nc"),
        "4.0h" => joinpath(simulation_path, "ta_4.0h_max.nc"),
    )

    @test simdir.vars["ta"]["max"] ==
          Dict{Any, Any}("3.0h" => nothing, "4.0h" => nothing)

    expected_files = Set(
        joinpath(simulation_path, f) for
        f in readdir(simulation_path) if endswith(f, ".nc")
    )

    @test simdir.allfiles == expected_files
end

@testset "OutputVar" begin

    simulation_path = joinpath(@__DIR__, "sample_data")

    simdir = ClimaAnalysis.SimDir(simulation_path)

    ta_max = get(simdir, short_name = "ta", reduction = "max", period = "4.0h")

    @test ta_max.attributes["units"] == "K"
    @test size(ta_max.var) == (3, 180, 80, 10)

    @test_throws ErrorException get(simdir, short_name = "bob")

    @test_throws ErrorException get(
        simdir,
        short_name = "ta",
        reduction = "max",
        period = "10.0h",
    )
end
