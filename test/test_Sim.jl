using Test
import ClimaAnalysis

@testset "SimDir" begin

    simulation_path = joinpath(@__DIR__, "sample_data")

    simdir = ClimaAnalysis.SimDir(simulation_path)

    @test ClimaAnalysis.available_vars(simdir) ==
          Set(["va", "ua", "orog", "ta", "ts", "pfull", "thetaa"])

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

    @test ClimaAnalysis.available_coord_types(
        simdir,
        short_name = "thetaa",
        reduction = "average",
        period = "1d",
    ) == Set([nothing, "pressure"])
    @test_throws ErrorException ClimaAnalysis.available_coord_types(
        simdir,
        short_name = "thetaa",
        reduction = "average",
        period = "2d",
    )

    @test simdir.variable_paths["orog"]["inst"][nothing][nothing] ==
          [joinpath(simulation_path, "orog_inst.nc")]
    @test simdir.variable_paths["ta"]["max"] == Dict{Any, Any}(
        "3.0h" =>
            Dict(nothing => [joinpath(simulation_path, "ta_3.0h_max.nc")]),
        "4.0h" =>
            Dict(nothing => [joinpath(simulation_path, "ta_4.0h_max.nc")]),
    )

    @test simdir.vars["ta"]["max"] == Dict{Any, Any}(
        "3.0h" => Dict{Any, Any}(nothing => nothing),
        "4.0h" => Dict{Any, Any}(nothing => nothing),
    )

    expected_files = Set(
        joinpath(simulation_path, f) for
        f in readdir(simulation_path) if endswith(f, ".nc")
    )

    @test simdir.allfiles == expected_files

    @test !isempty(simdir)
    @test isempty(ClimaAnalysis.SimDir(mktempdir()))
end

@testset "Show and summary for Simdir" begin
    simulation_path = joinpath(@__DIR__, "sample_data")

    simdir = ClimaAnalysis.SimDir(simulation_path)

    @test sprint(show, simdir) ==
          "Output directory: " *
          simulation_path *
          "\nVariables:\n- va\n    average (2.0h)\n- pfull\n    inst (1.0d)\n- ua\n    average (6.0h)\n- orog\n    inst (nothing)\n- thetaa\n    average (1d)\n- ta\n    average (3.0h)\n    max (4.0h, 3.0h)\n    min (3.0h)\n- ts\n    max (1.0h)"
    @test sprint(summary, simdir) == sprint(show, simdir)

    # Empty simdir
    simdir = ClimaAnalysis.SimDir{Dict{Any, Any}, Dict{Any, Any}}(
        "path",
        Dict(),
        Dict(),
        Set(),
    )
    @test sprint(summary, simdir) == "Output directory: path\nVariables:"
end

@testset "OutputVar" begin

    simulation_path = joinpath(@__DIR__, "sample_data")

    simdir = ClimaAnalysis.SimDir(simulation_path)

    ta_max = get(simdir, short_name = "ta", reduction = "max", period = "4.0h")
    ts_max = get(simdir, short_name = "ts", reduction = "max", period = "1.0h")

    @test ClimaAnalysis.units(ta_max) == "K"

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
    @test_throws ErrorException orog =
        get(simdir, short_name = "orog", reduction = "inst", period = "4.0h")
    orog = get(simdir, short_name = "orog", reduction = "inst")

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
          get(simdir; short_name = "orog", reduction = "inst", period = nothing)
    @test ts_max == get(simdir; short_name = "ts", reduction = "max")

    # Test `get` with multiple arguments
    orog2, pfull2 = get(simdir, "orog", "pfull")
    @test orog == orog2
    @test pfull == pfull2

    # short_name, long_name, units
    @test ClimaAnalysis.short_name(orog) == "orog"
    @test ClimaAnalysis.long_name(orog) == "Surface Altitude, Instantaneous"
    @test ClimaAnalysis.units(orog) == "m"

    # The ts_max file bundled with ClimaAnalysis does not have the short_name
    # attribute because it was generated before that was a feature. We use that
    # to check the empty tring
    @test ClimaAnalysis.short_name(ts_max) == ""

    # Test coord_type keyword argument
    @test_throws ErrorException get(
        simdir,
        short_name = "thetaa",
        reduction = "average",
    )
    @test_throws ErrorException get(
        simdir,
        short_name = "thetaa",
        reduction = "average",
        coord_type = "not_pfull",
    )
    thetaa_var = get(
        simdir,
        short_name = "thetaa",
        reduction = "average",
        coord_type = nothing,
    )
    thetaa_var_pressure_coords = get(
        simdir,
        short_name = "thetaa",
        reduction = "average",
        coord_type = "pressure",
    )
    @test !ClimaAnalysis.has_pressure(thetaa_var)
    @test ClimaAnalysis.has_pressure(thetaa_var_pressure_coords)
end

@testset "SimDir with stitching" begin
    simulation_path = joinpath(@__DIR__, "sample_outputs")

    simdir = ClimaAnalysis.SimDir(simulation_path)

    @test ClimaAnalysis.available_vars(simdir) ==
          Set(["pfull", "orog", "ts", "test"])

    @test ClimaAnalysis.available_reductions(simdir, short_name = "pfull") ==
          Set(["inst"])
    @test_throws ErrorException ClimaAnalysis.available_reductions(
        simdir,
        short_name = "bob",
    )

    @test ClimaAnalysis.available_periods(
        simdir,
        short_name = "pfull",
        reduction = "inst",
    ) == Set(["2.0d"])
    @test_throws ErrorException ClimaAnalysis.available_periods(
        simdir,
        short_name = "pfull",
        reduction = "bob",
    )

    @test simdir.variable_paths["ts"]["max"]["1.0h"][nothing] == [
        joinpath(simulation_path, "output_0001", "ts_1.0h_max.nc"),
        joinpath(simulation_path, "output_0002", "ts_1.0h_max.nc"),
    ]
    @test simdir.variable_paths["pfull"]["inst"]["2.0d"][nothing] == [
        joinpath(simulation_path, "output_0001", "pfull_2.0d_inst.nc"),
        joinpath(simulation_path, "output_0002", "pfull_2.0d_inst.nc"),
        joinpath(simulation_path, "output_0003", "pfull_2.0d_inst.nc"),
    ]

    @test simdir.vars["pfull"]["inst"] ==
          Dict{Any, Any}("2.0d" => Dict(nothing => nothing))

    expected_files = Set()
    for (root, _, files) in walkdir(simulation_path)
        for file in files
            if endswith(file, ".nc")
                push!(expected_files, joinpath(root, file))
            end
        end
    end

    @test simdir.allfiles == expected_files

    @test !isempty(simdir)
end

@testset "OutputVar with stitching" begin
    simulation_path = joinpath(@__DIR__, "sample_outputs")

    simdir = ClimaAnalysis.SimDir(simulation_path)

    # We do not care about what these variable mean, but only that the stitching works
    # Normal case: points in time are in sequential order
    orog_var = get(simdir, short_name = "orog", reduction = "inst")

    # Wrong names for time dimension
    @test_throws ErrorException get(
        simdir,
        short_name = "ts",
        reduction = "max",
        period = "1.0h",
    )

    # Time dimension is not in increasing order
    @test_throws ErrorException get(
        simdir,
        short_name = "pfull",
        reduction = "inst",
        period = "2.0d",
    )

    # Start dates do not match
    @test_throws ErrorException get(
        simdir,
        short_name = "ts",
        reduction = "average",
        period = "10d",
    )

    # Time dimension does not exists
    @test_throws ErrorException get(
        simdir,
        short_name = "test",
        reduction = "average",
        period = "10d",
    )

    @test orog_var |> ClimaAnalysis.times == [0.0, 1.0, 2.0]
    @test orog_var.data |> size == (3, 180, 90)
end
