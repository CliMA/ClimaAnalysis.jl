using Test
import ClimaAnalysis

@testset "SimDir" begin

    simulation_path = joinpath(@__DIR__, "sample_data")

    simdir = ClimaAnalysis.SimDir(simulation_path)

    @test ClimaAnalysis.available_vars(simdir) ==
          Set(["va", "ua", "orog", "ta", "ts", "pfull"])

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
    ts_max = get(simdir, short_name = "ts", reduction = "max", period = "1.0h")

    @test ta_max.attributes["units"] == "K"

    @test collect(keys(ta_max.dims)) == Array(["time", "lon", "lat", "z"])

    @test ta_max.dim2index ==
          Dict(["time" => 1, "lon" => 2, "lat" => 3, "z" => 4])

    @test ta_max.index2dim == ["time", "lon", "lat", "z"]

    @test ta_max.dim_attributes["lat"]["units"] == "degrees_north"

    @test size(ta_max.data) == (3, 180, 80, 10)

    @test_throws ErrorException get(simdir, short_name = "bob")

    @test_throws ErrorException get(
        simdir,
        short_name = "ta",
        reduction = "max",
        period = "10.0h",
    )

    # Check if reduction = inst and passing a period
    orog = get(simdir, short_name = "orog", reduction = "inst", period = "4.0h")
    @test orog == get(simdir, short_name = "orog", reduction = "inst")

    # Check is_z_1d
    pfull = get(simdir, short_name = "pfull", reduction = "inst")
    @test ClimaAnalysis.is_z_1D(pfull) == false

    # Check the shorter gets
    @test_throws ErrorException get(simdir, "ta")
    @test_throws ErrorException get(simdir; short_name = "ta")
    @test_throws ErrorException get(
        simdir;
        short_name = "ta",
        reduction = "max",
    )

    @test orog == get(simdir, "orog")
    @test orog ==
          get(simdir; short_name = "orog", reduction = "inst", period = "4.0h")
    @test ts_max == get(simdir; short_name = "ts", reduction = "max")

    # short_name, long_name, units
    @test ClimaAnalysis.short_name(orog) == "orog"
    @test ClimaAnalysis.long_name(orog) == "Surface Altitude, Instantaneous"
    @test ClimaAnalysis.units(orog) == "m"

    # The ts_max file bundled with ClimaAnalysis does not have the short_name
    # attribute because it was generated before that was a feature. We use that
    # to check the empty tring
    @test ClimaAnalysis.short_name(ts_max) == ""
end
