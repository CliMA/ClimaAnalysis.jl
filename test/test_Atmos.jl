using Test
import ClimaAnalysis
import OrderedCollections: OrderedDict

@testset "To pressure coordinates" begin

    # Let start with testing a single column
    z_alt = 0:100.0 |> collect
    data = copy(z_alt)

    zvar = ClimaAnalysis.OutputVar(Dict("z" => z_alt), data)

    # Fake pressure, linearly decreasing, so that we can check precisely
    pressure = 300.0:-2.0:100.0 |> collect
    pdata = copy(pressure)

    attribs = Dict("short_name" => "pfull")
    dim_attribs = Dict{String, Any}()
    pressure_var =
        ClimaAnalysis.OutputVar(attribs, Dict("z" => z_alt), dim_attribs, pdata)

    pressure_in_pressure_coordinates =
        ClimaAnalysis.Atmos.to_pressure_coordinates(pressure_var, pressure_var)

    @test collect(keys(pressure_in_pressure_coordinates.dims)) == ["pfull"]
    # reverse because we go from min to max for pressure (to have it increasing
    # for Interpolations.jl)
    @test pressure_in_pressure_coordinates.dims["pfull"] == reverse(pdata)
    @test pressure_in_pressure_coordinates.data == reverse(pdata)

    # Fake var, again linear. When everything is linear we should obtain the
    # input variable back

    var_func(z) = 500 + 10 * z
    myvardata = var_func.(z_alt)
    mydata = copy(myvardata)

    attribs = Dict("short_name" => "myvar")
    dim_attribs = Dict{String, Any}()
    myvar = ClimaAnalysis.OutputVar(
        attribs,
        Dict("z" => z_alt),
        dim_attribs,
        mydata,
    )

    myvar_in_pressure_coordinates =
        ClimaAnalysis.Atmos.to_pressure_coordinates(myvar, pressure_var)

    @test collect(keys(myvar_in_pressure_coordinates.dims)) == ["pfull"]
    @test myvar_in_pressure_coordinates.dims["pfull"] == reverse(pdata)
    @test myvar_in_pressure_coordinates.data == reverse(mydata)

    exp_pressure = exp.(-0.01 .* z_alt |> collect)
    exp_pdata = copy(exp_pressure)

    attribs = Dict("short_name" => "pfull")
    dim_attribs = Dict{String, Any}()
    exp_pressure_var = ClimaAnalysis.OutputVar(
        attribs,
        Dict("z" => z_alt),
        dim_attribs,
        exp_pdata,
    )

    myvar_in_exp_pressure_coordinates =
        ClimaAnalysis.Atmos.to_pressure_coordinates(myvar, exp_pressure_var)

    # Linear range from min to max
    expected_range = collect(
        range(
            minimum(exp_pressure),
            maximum(exp_pressure),
            length = length(exp_pdata),
        ),
    )
    @test myvar_in_exp_pressure_coordinates.dims["pfull"] == expected_range

    # Specify pressure levels
    @test ClimaAnalysis.Atmos.to_pressure_coordinates(
        myvar,
        exp_pressure_var,
        target_pressure = [1.0],
    ).data[1] == mydata[1]

    # From the pressure range we can compute the corresponding altitudes
    corresponding_z = -100 .* log.(expected_range)
    @test myvar_in_exp_pressure_coordinates.data ≈ var_func.(corresponding_z) rtol =
        1e-5

    # 3D test
    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect
    zzz = 0.0:10.0 |> collect
    var3D_func(lg, lt, zz) = 2lg + 3lt + 4zz
    data = [var3D_func(lg, lt, zz) for lg in long, lt in lat, zz in zzz]
    dims = OrderedDict(["lon" => long, "lat" => lat, "z" => zzz])
    dim_attribs = Dict{String, Any}()
    var3D = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    altitude = [zz for lg in long, lt in lat, zz in zzz]
    pdata3D = exp.(-0.01 .* altitude)

    attribs = Dict("short_name" => "pfull")
    dim_attribs = Dict{String, Any}()
    pressure3D = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, pdata3D)

    my3Dvar_in_exp_pressure_coordinates =
        ClimaAnalysis.Atmos.to_pressure_coordinates(var3D, pressure3D)

    # Linear range from min to max
    overall_range =
        collect(range(minimum(pdata3D), maximum(pdata3D), length = length(zzz)))
    @test my3Dvar_in_exp_pressure_coordinates.dims["pfull"] == overall_range

    # From the pressure range we can compute the corresponding altitudes
    expected_output = [
        var3D_func.(lg, lt, -100 .* log(p)) for lg in long, lt in lat,
        p in overall_range
    ]
    @test my3Dvar_in_exp_pressure_coordinates.data ≈ expected_output rtol = 1e-5

    # Error checking
    @test_throws ErrorException ClimaAnalysis.Atmos.to_pressure_coordinates(
        var3D,
        pressure_var,
    )
    # This will fail because there is no recognized altitude dimension
    Zvar = ClimaAnalysis.OutputVar(Dict("Z" => z_alt), z_alt)
    @test_throws ErrorException ClimaAnalysis.Atmos.to_pressure_coordinates(
        Zvar,
        Zvar,
    )
end

@testset "RMSE of pressure coordinates" begin
    lon = 0.5:1.0:359.5 |> collect
    lat = -89.5:1.0:89.5 |> collect
    zzz = 1.0:10.0 |> collect
    var3D_func(lg, lt, zz) = zz
    data = [var3D_func(lg, lt, zz) for lg in lon, lt in lat, zz in zzz]
    dims = OrderedDict(["lon" => lon, "lat" => lat, "z" => zzz])
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "degrees_east"),
        "lat" => Dict("units" => "degrees_north"),
        "z" => Dict("units" => "m"),
    ])
    attribs = Dict("short_name" => "T", "units" => "K")
    var3D = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    altitude = [zz for lg in lon, lt in lat, zz in zzz]
    pdata3D = -1.0 .* altitude

    pfull_attribs = Dict("short_name" => "pfull")
    pfull_dim_attribs = OrderedDict([
        "lon" => Dict("units" => "degrees_east"),
        "lat" => Dict("units" => "degrees_north"),
        "z" => Dict("units" => "m"),
    ])
    pfull_attribs = Dict("short_name" => "pfull", "units" => "kg m^-2 s^-2")
    pressure3D =
        ClimaAnalysis.OutputVar(pfull_attribs, dims, pfull_dim_attribs, pdata3D)

    global_rmse_pfull = ClimaAnalysis.Atmos.global_rmse_pfull(
        var3D,
        var3D,
        sim_pressure = pressure3D,
        obs_pressure = pressure3D,
    )
    @test global_rmse_pfull == 0.0

    # Test if the computation is the same as a manual computation
    zero_data = zeros(size(data))
    zero_var3D = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, zero_data)
    pressure3D =
        ClimaAnalysis.OutputVar(pfull_attribs, dims, pfull_dim_attribs, pdata3D)
    zero_var =
        ClimaAnalysis.Atmos.to_pressure_coordinates(zero_var3D, pressure3D)

    global_rmse_pfull = ClimaAnalysis.Atmos.global_rmse_pfull(
        var3D,
        zero_var,
        sim_pressure = pressure3D,
    )
    min_pfull, max_pfull = extrema(zero_var.dims["pfull"])

    # (4pi * (1:10 .|> x -> ^(x, 2))) is from finding the squared error and integrating with
    # respect to lon and lat for each pressure level
    # Multiply by (1 / 4pi) to normalize with respect to lon and lat
    # Summing is integrating with respect to pressure
    # Divide by (max_pfull - min_pfull + 1.0) to normalize with respect to pressure which gives
    # MSE; the 1.0 comes from the pressure levels being equispaced
    # Take square root to get RMSE
    global_rmse_test = sqrt(
        sum(((1 / 4pi) * 4pi * (1:10 .|> x -> ^(x, 2)))) /
        (max_pfull - min_pfull + 1.0),
    )
    @test global_rmse_pfull ≈ global_rmse_test

    # Test with data from era5; checking if we get no error
    ncpath = joinpath(@__DIR__, "sample_nc/test_pfull.nc")
    pfull_obs_var = ClimaAnalysis.OutputVar(ncpath, "t")
    ClimaAnalysis.set_dim_units!(
        pfull_obs_var,
        "pressure_level",
        "kg m^-2 s^-2",
    )
    pfull_obs_var =
        ClimaAnalysis.window(pfull_obs_var, "pressure_level", left = -10)
    global_rmse_pfull = ClimaAnalysis.Atmos.global_rmse_pfull(
        var3D,
        pfull_obs_var,
        sim_pressure = pressure3D,
    )
    @test global_rmse_pfull >= 0.0

    # Error handling
    # OutputVars are not in pressure coordinates and no pressure OutputVar is supplied
    @test_throws ErrorException ClimaAnalysis.Atmos.global_rmse_pfull(
        var3D,
        var3D,
        sim_pressure = nothing,
        obs_pressure = pressure3D,
    )
    @test_throws ErrorException ClimaAnalysis.Atmos.global_rmse_pfull(
        var3D,
        var3D,
        sim_pressure = pressure3D,
        obs_pressure = nothing,
    )

    # No pressure dimension
    lat = collect(range(-89.5, 89.5, 180))
    lon = collect(range(-179.5, 179.5, 360))
    data = ones(length(lat), length(lon))
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "lon" => Dict("units" => "deg"),
    ])
    var_latlon = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test_throws ErrorException ClimaAnalysis.Atmos.global_rmse_pfull(
        var_latlon,
        var3D,
        obs_pressure = pressure3D,
    )
    @test_throws ErrorException ClimaAnalysis.Atmos.global_rmse_pfull(
        var3D,
        var_latlon,
        sim_pressure = pressure3D,
    )
end
