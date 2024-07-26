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
