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
    myvardata = 500.0:10.0:1500.0 |> collect
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

    exp_pressure = exp.(1:-0.01:0 |> collect)
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
    @test myvar_in_exp_pressure_coordinates.dims["pfull"] == collect(range(1, exp(1), length = length(exp_pdata)))
    # TODO: Add tests here

    # 3D test
    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect
    zzz = 0.0:10.0 |> collect
    data = reshape(1.0:(91 * 181 * 11), (181, 91, 11))
    dims = OrderedDict(["lon" => long, "lat" => lat, "z" => zzz])
    dim_attribs = Dict{String, Any}()
    var3D = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    pdata3D = data.^-1

    attribs = Dict("short_name" => "pfull")
    dim_attribs = Dict{String, Any}()
    pressure3D = ClimaAnalysis.OutputVar(
        attribs,
        dims,
        dim_attribs,
        pdata3D,
    )

    ClimaAnalysis.Atmos.to_pressure_coordinates(var3D, pressure3D)

end
