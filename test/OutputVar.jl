using Test
import ClimaAnalysis

import Statistics: mean
import OrderedCollections: OrderedDict

@testset "General" begin
    # Add test for short constructor
    long = -180.0:180.0 |> collect
    data = copy(long)

    longvar = ClimaAnalysis.OutputVar(Dict("long" => long), data)

    @test longvar.dims["long"] == long


    # center_longitude!
    #
    # Check var without long
    dims = Dict("z" => long)
    var_error = ClimaAnalysis.OutputVar(
        Dict{String, Any}(),
        dims,
        Dict{String, Any}(),
        data,
    )
    @test_throws ErrorException ClimaAnalysis.center_longitude!(
        var_error,
        180.0,
    )

    time = 0:10.0 |> collect
    dims = Dict("lon" => long, "time" => time)
    data = collect(reshape(1:(361 * 11), (361, 11)))
    var_good = ClimaAnalysis.OutputVar(
        Dict{String, Any}(),
        dims,
        Dict{String, Any}(),
        data,
    )
    ClimaAnalysis.center_longitude!(var_good, 90.0)
    # We are shifting by 91
    @test var_good.dims["lon"][180] == 90
    @test var_good.data[3, :] == data[3, :]
    @test var_good.data[180, 1] == 271
end

@testset "Arithmetic operations" begin
    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect
    time = 0.0:10.0 |> collect

    data1 = reshape(1.0:(91 * 181 * 10), (10, 181, 91))

    dims = OrderedDict(["time" => time, "lon" => long, "lat" => lat])
    dim_attributes = OrderedDict([
        "time" => Dict(),
        "lon" => Dict("b" => 2),
        "lat" => Dict("a" => 1),
    ])
    attribs = Dict("short_name" => "bob", "long_name" => "hi")
    var1 = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data1)

    dim_attributes2 = OrderedDict([
        "time" => Dict(),
        "lon" => Dict("lol" => 2),
        "lat" => Dict("a" => 1),
    ])

    var2 = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes2, data1)

    data3 = 5.0 .+ reshape(1.0:(91 * 181 * 10), (10, 181, 91))
    attribs3 = Dict("long_name" => "bob", "short_name" => "bula")
    var3 = ClimaAnalysis.OutputVar(attribs3, dims, dim_attributes, data3)

    # Check arecompatible
    @test !ClimaAnalysis.arecompatible(var1, var2)
    @test ClimaAnalysis.arecompatible(var1, var3)

    var1plus10 = var1 + 10

    @test var1plus10.data == data1 .+ 10
    @test ClimaAnalysis.short_name(var1plus10) == "bob + 10"
    @test ClimaAnalysis.long_name(var1plus10) == "hi + 10"

    tenplusvar1 = 10 + var1

    @test tenplusvar1.data == data1 .+ 10
    @test ClimaAnalysis.short_name(tenplusvar1) == "10 + bob"
    @test ClimaAnalysis.long_name(tenplusvar1) == "10 + hi"

    var1plusvar3 = var1 + var3

    @test var1plusvar3.data == data1 .+ data3
    @test ClimaAnalysis.short_name(var1plusvar3) == "bob + bula"
    @test ClimaAnalysis.long_name(var1plusvar3) == "hi + bob"
end

@testset "Reductions" begin
    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect
    time = 0.0:10.0 |> collect

    data = reshape(1.0:(91 * 181 * 10), (10, 181, 91))

    dims = OrderedDict(["time" => time, "lon" => long, "lat" => lat])
    dim_attributes = OrderedDict([
        "time" => Dict(),
        "lon" => Dict("b" => 2),
        "lat" => Dict("a" => 1),
    ])
    attribs = Dict("long_name" => "hi")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)

    # Test copy
    var_copied = copy(var)
    fields = fieldnames(ClimaAnalysis.OutputVar)
    for field in fields
        @test getfield(var, field) == getfield(var_copied, field)
        @test getfield(var, field) !== getfield(var_copied, field)
    end

    # Test reduction
    lat_avg = ClimaAnalysis.average_lat(var)
    @test lat_avg.dims == OrderedDict(["lon" => long, "time" => time])
    @test lat_avg.dim_attributes ==
          OrderedDict(["lon" => Dict("b" => 2), "time" => Dict()])
    @test lat_avg.data == dropdims(mean(data, dims = 3), dims = 3)

    lat_lon_avg = ClimaAnalysis.average_lon(lat_avg)
    @test lat_lon_avg.dims == OrderedDict(["time" => time])
    @test lat_lon_avg.dim_attributes == OrderedDict(["time" => Dict()])

    @test lat_lon_avg.data == dropdims(mean(lat_avg.data, dims = 2), dims = 2)

    lat_lon_time_avg = ClimaAnalysis.average_time(lat_lon_avg)
    @test lat_lon_time_avg.dims == OrderedDict()
    @test lat_lon_time_avg.dim_attributes == OrderedDict()

    @test lat_lon_time_avg.data[] == mean(data)

    @test lat_lon_time_avg.attributes["long_name"] ==
          "hi averaged over latitudes averaged over longitudes averaged over time"
end

@testset "Slicing" begin
    z = 0.0:20.0 |> collect
    time = 100.0:110.0 |> collect

    data = reshape(1.0:(11 * 21), (11, 21))

    dims = OrderedDict(["time" => time, "z" => z])
    dim_attributes =
        OrderedDict(["time" => Dict("units" => "s"), "z" => Dict("b" => 2)])
    attribs = Dict("long_name" => "hi")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)

    z_sliced = ClimaAnalysis.slice(var, z = 1.0)
    # 1.0 is the second index
    z_expected_data = data[:, 2]
    @test z_sliced.dims == OrderedDict(["time" => time])
    @test z_sliced.dim_attributes ==
          OrderedDict(["time" => Dict("units" => "s")])
    @test z_sliced.data == z_expected_data

    t_sliced = ClimaAnalysis.slice(var, time = 200.0)
    # 200 is the last index
    t_expected_data = data[end, :]
    @test t_sliced.dims == OrderedDict(["z" => z])
    @test t_sliced.dim_attributes == OrderedDict(["z" => Dict("b" => 2)])
    @test t_sliced.data == t_expected_data

    @test t_sliced.attributes["long_name"] == "hi time = 1m 50.0s"

    # Test with the general slice

    t_sliced = ClimaAnalysis.slice(var, time = 200.0)
    # 200 is the last index
    t_expected_data = data[end, :]
    @test t_sliced.dims == OrderedDict(["z" => z])
    @test t_sliced.dim_attributes == OrderedDict(["z" => Dict("b" => 2)])
    @test t_sliced.data == t_expected_data

    @test t_sliced.attributes["long_name"] == "hi time = 1m 50.0s"

    @test t_sliced.attributes["slice_time"] == "110.0"
    @test t_sliced.attributes["slice_time_units"] == "s"
end

@testset "Windowing" begin
    z = 0.0:20.0 |> collect
    time = 0.0:10.0 |> collect

    data = reshape(1.0:(11 * 21), (11, 21))

    dims = OrderedDict(["time" => time, "z" => z])
    dim_attributes =
        OrderedDict(["time" => Dict("units" => "s"), "z" => Dict("b" => 2)])
    attribs = Dict("long_name" => "hi")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)

    # Dimension not existing
    @test_throws ErrorException ClimaAnalysis.window(var, "lat")

    # Left right not ordered
    @test_throws ErrorException ClimaAnalysis.window(
        var,
        "time",
        left = 10,
        right = 1,
    )

    var_windowed = ClimaAnalysis.window(var, "time", left = 2.5, right = 5.1)
    expected_data = data[3:6, :]

    @test var_windowed.data == expected_data

    @test var_windowed.dims["time"] == time[3:6]
end
